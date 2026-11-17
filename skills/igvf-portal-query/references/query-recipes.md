# IGVF Portal query recipes

Runnable starting points for the three question shapes this skill covers. All of these
assume the setup from `../SKILL.md`:

```python
import sys
REPO = "/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks-igvf-portal-submission"
sys.path.insert(0, f"{REPO}/workflow/scripts")
from igvf_metadata import portal_client

reader = portal_client.PortalReader(igvf_mode="prod")
```

Run with `/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python`,
inside a Slurm allocation. None of these write anything — every call below is a GET.

## 1. Single-record lookup: "what is this alias/accession?"

Use when you have one alias or accession and want to know what it is, resolve it to a
real `@id`, or check a few fields on it.

```python
record = reader.get_by_alias("jesse-engreitz:some_alias")  # or a bare accession, e.g. "IGVFFI5251UJRB"
if record is None:
    print("not found")
else:
    print(record["@type"], record["@id"])
```

**Probe before you build logic on a field's shape.** The same field name can come back as
either a bare path string or an embedded dict depending on whether it's a plain `linkTo` or
an embedded link, and you cannot tell from the field name alone. Before writing matching
logic, fetch one record and print the field:

```python
import json
print(json.dumps(record, indent=2))  # or just record.get("the_field_you_care_about")
```

This is cheap (one record) and saves you from writing filter logic against a guessed shape
that turns out wrong on row 2.

## 2. Full paged search: "how many / which objects match X?"

Use `get_search_paged` for anything beyond a single record. It already paces itself
(200 rows/page, ~10s between pages, retried on transient errors and 5xx) — this is
deliberate, not a default to tune down. Don't wrap it in your own retry loop or shrink the
pause; the pacing exists because a Portal under submission-deadline load has broken from
tighter polling before.

```python
params = [
    ("type", "PredictionSet"),
    ("software.title", "scE2G"),
    ("field", "@id"),
    ("field", "status"),
    ("field", "submitted_by"),
    ("field", "cell_type"),
]
rows = reader.get_search_paged(params)
print(f"{len(rows)} matching PredictionSets")  # verified against the Portal's own total
```

If the Portal is slow (it can be, especially under load), run this in the background rather
than blocking on it — a few hundred rows across 2-3 pages, each with a courtesy pause, can
take a couple of minutes. That's expected, not a hang.

### Count only, no rows needed

If you genuinely don't need the rows — just "how many" — a single `limit=0` page is lighter:

```python
import urllib.parse, requests
import igvf_utils.utils as iuu

conn = reader._connection()
url = iuu.url_join([conn.igvf_mode.url, "search/?"])
query = urllib.parse.urlencode(params + [("limit", "0")])
resp = requests.get(url + query, auth=conn.auth, timeout=60, headers=iuu.REQUEST_HEADERS_JSON)
resp.raise_for_status()
print(resp.json()["total"])
```

But if there's any chance you'll also want a breakdown (by status, by submitter, by cell
type...), just call `get_search_paged` once and derive everything — including the count —
from that one result. Re-querying per question you ask yourself is the thing to avoid.

## 3. Resolve an alias, then filter a bulk search by it

The common pattern for "which objects reference thing Y, and how does that split by Z":
resolve Y to its real `@id` once (a single `get_by_alias` call), then test membership
against the array field client-side, because the alias itself doesn't live in the field
you're filtering.

```python
import collections

target_id = reader.get_by_alias("jesse-engreitz:some_document_alias")["@id"]

params = [
    ("type", "PredictionSet"),
    ("software.title", "scE2G"),
    ("field", "@id"),
    ("field", "documents"),      # bare path strings, e.g. ["/documents/<uuid>/", ...]
    ("field", "status"),
    ("field", "submitted_by"),   # embedded dict: {"@id": "...", "title": "..."}
]
rows = reader.get_search_paged(params)

matching = [r for r in rows if target_id in (r.get("documents") or [])]
print(f"{len(matching)} / {len(rows)} reference it")

by_status = collections.Counter(r.get("status") for r in matching)
by_submitter = collections.Counter(
    (r.get("submitted_by") or {}).get("@id") for r in matching
)
for label, counts in [("status", by_status), ("submitted_by", by_submitter)]:
    print(f"\n{label} breakdown:")
    for k, n in counts.most_common():
        print(f"  {k}: {n}")
```

## 4. Unique-value tally / group-by, including "mine vs everyone's"

One fetch with the minimal fields, then `collections.Counter` does all the grouping —
including comparing your own subset against the full set, which is just a filter on the
same in-memory rows, not a second query.

```python
params = [
    ("type", "PredictionSet"),
    ("software.title", "scE2G"),
    ("field", "@id"),
    ("field", "cell_type"),
    ("field", "submitted_by"),
]
rows = reader.get_search_paged(params)

def cell_type_tally(subset):
    return collections.Counter((r.get("cell_type") or {}).get("@id") for r in subset)

all_tally = cell_type_tally(rows)
print(f"{len(all_tally)} unique cell types across all {len(rows)} PredictionSets")

MY_USER_ID = "/users/<your-uuid>/"  # get this once from any record you know you submitted
mine = [r for r in rows if (r.get("submitted_by") or {}).get("@id") == MY_USER_ID]
mine_tally = cell_type_tally(mine)
print(f"{len(mine_tally)} unique cell types across your {len(mine)} PredictionSets")
```

Don't know your own `submitted_by.@id` yet? Pull it from the first row of a query you know
is yours, or from `get_by_alias` on anything you submitted — it's in `submitted_by["@id"]`
on any record. No need to search for it separately.

## Credential safety, as it applies to all of the above

None of these recipes ever touch `IGVF_API_KEY` / `IGVF_SECRET_KEY` directly — they only
ever use `conn.auth` (built internally by `igvf_utils.connection.Connection` from the
environment) or `reader.auth`, passed straight into `requests`/`Connection.get()`. Keep it
that way in anything you write from these patterns:

- Never print, log, or write `conn.auth`, `reader.auth`, or the env vars themselves.
- Never run `env` / `printenv` / dump the environment to inspect whether the keys are set —
  if you need to check they're present, do it inside an sbatch script, not interactively.
- Never read `~/.bash_profile`, `~/.bashrc`, `~/.netrc`, or any keypair file.
- If a key ever ends up somewhere it shouldn't (a print, a log file, a committed file), say
  so immediately and loudly — this repo is public.
