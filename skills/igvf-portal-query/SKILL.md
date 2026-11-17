---
name: igvf-portal-query
description: Query the IGVF Data Portal read-only to look things up, count how many objects match a description, and cross-tabulate metadata — e.g. "how many PredictionSets use software X", "count unique cell types across PredictionSets, mine vs everyone's", "which object type does alias Y show up in", "resolve this accession and tell me what it is". Use whenever asked to check the Portal, look up an alias/accession, count or group-by a field across Portal objects, or otherwise answer "what does the Portal say about X" — purely as discovery, never to upload or change anything. For generating upload manifests or running the QC pipeline, use qc-and-predictions instead; this skill only covers the read side.
allowed-tools: Read, Bash, Grep, Glob
---

# IGVF Portal Query: read-only discovery, counting, and cross-tabs

This is the read-only sibling of `qc-and-predictions`: that skill builds and uploads
manifests, this one only looks things up on the Portal — resolving an alias, counting
how many objects match a filter, or tallying a field across a set of objects. It shares
the same hard rules (never upload, never expose credentials) but never touches
`iu_register.py` at all, because it has no write path to guard.

The actual client code already exists and is proven — `workflow/scripts/igvf_metadata/portal_client.py`
in the `QC_pseudobulks-igvf-portal-submission` repo (`$REPO` below). Don't re-implement
`requests.get` loops by hand; that file already handles auth, retry/backoff, and
paging correctly, and re-deriving those from scratch is how credentials get logged by
accident or the Portal gets hammered during a submission deadline. Import it and build
small scripts on top of `PortalReader`.

```
REPO=/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks-igvf-portal-submission
```

## Hard rules (non-negotiable)

- **Read-only.** Never call `iu_register.py`, never pass `--mode upload` or `--mode validate`,
  never POST/PATCH anything. Every query in this skill is a GET. If a question seems to need
  a write to answer (e.g. "did this upload succeed"), that's `get_by_alias` with
  `database=True` re-reading the Portal's own current state, not a write.
- **Never expose `IGVF_API_KEY` / `IGVF_SECRET_KEY`.** Never print, echo, log, or write them
  anywhere; never run `env`/`printenv`; never read `~/.bash_profile`, `~/.bashrc`, `~/.netrc`,
  or any keypair file. The only thing calling code ever touches is `conn.auth` — an opaque
  tuple `igvf_utils.connection.Connection` builds internally from the environment — passed
  straight into `requests.get(..., auth=conn.auth)` or used implicitly by `Connection.get()`.
  If a key is ever accidentally exposed (printed, logged, committed), say so loudly and
  immediately; this repo is public.
- **Sherlock: no network Python on a login node.** Run inside a Slurm allocation
  (`sh_dev` or `sbatch`) — check with `hostname` / `$SLURM_JOB_ID` first. Use the
  `igvf_utils_env` conda environment, which has `igvf_utils` installed:
  `/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python`.
  There's no installed package for `portal_client`, so `sys.path.insert(0, f"{REPO}/workflow/scripts")`
  before importing it.
- **Be nice to the Portal.** `get_search_paged` already paces itself (200 rows/page, a pause
  between pages, generous per-page timeouts, retry-with-backoff on transient errors and on
  5xx). Reuse it rather than inventing a tighter loop — the existing pacing exists because an
  impatient retry loop has caused real multi-hour incidents against a Portal the DACC is
  often actively absorbing a submission deadline on.

## Setup

```python
import sys
sys.path.insert(0, f"{REPO}/workflow/scripts")
from igvf_metadata import portal_client

reader = portal_client.PortalReader(igvf_mode="prod")  # or "staging" / "sandbox"
```

Run this with the `igvf_utils_env` Python, inside a Slurm allocation:

```
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python your_script.py
```

Prefer running it in the background (if your harness supports it) rather than blocking —
a paged search over a few hundred rows can take a couple of minutes once you include the
courtesy pauses between pages, and that's by design, not a hang.

## The three read primitives on `PortalReader`

Pick based on what you're asking:

| Question shape | Use |
|---|---|
| "What is this one alias/accession?" | `reader.get_by_alias(alias_or_accession)` |
| "How many / which objects match type + filters?" | `reader.get_search_paged(params)` |
| "I need fields `/search/` can't give me together" | `reader.get_multireport(query_string)` |

**`get_by_alias(alias, database=False)`** — one object, by alias or accession. Returns the
record dict, or `None` if it doesn't exist. Reads the search index by default; pass
`database=True` only if you need to see a write that just happened (the index can lag a
live PATCH by several seconds).

**`get_search_paged(params, page_size=200, pause_seconds=10, timeout=180)`** — every row
matching a filter set, fetched a page at a time. `params` is a list of `(key, value)`
pairs, same shape as the Portal's own query string, *without* `limit`/`from`/`sort` (those
are added for you). It raises `RuntimeError` if the rows it collected don't add up to the
Portal's own reported `total` — so a result from this function is never a silent partial
page; trust it or get an exception, not a quiet undercount.

**`get_multireport(query_string)`** — a raw GET against `/multireport/`, for the rare case
where you need fields `Connection.search()`/`get_search_paged` can't reach because
`Connection.search()` hardcodes the `search/` path. You'll rarely need this; reach for
`get_search_paged` first.

## Building a filter set

`params` for `get_search_paged` (or the query string for `get_multireport`) is just Portal
search syntax as `(key, value)` pairs:

```python
params = [
    ("type", "PredictionSet"),
    ("software.title", "scE2G"),
    ("field", "@id"),           # repeat "field" per field you want back — keeps pages light
    ("field", "cell_type"),
    ("field", "submitted_by"),
]
```

A few things that aren't obvious until you've hit them:

- **Ask for only the fields you need.** Every `("field", "...")` pair narrows the response.
  Fetching full records for a few hundred rows is slower and heavier than it needs to be for
  a counting or tallying question.
- **Dotted paths filter through a link**, e.g. `software.title=scE2G` or
  `submitted_by.@id=/users/<uuid>/`. The `.@id` suffix matters — a plain
  `submitted_by=/users/<uuid>/` 404s. If a filter 404s, suspect the field needs `.@id` (or
  `.uuid`, or another sub-property) rather than the bare link value.
- **Embedded fields come back as dicts with their own `@id`; plain `linkTo` fields come back
  as bare path strings.** For example `submitted_by` and `cell_type` typically arrive as
  `{"@id": "/users/...", "title": "..."}`, but `documents` arrives as a flat list of path
  strings: `["/documents/<uuid>/", ...]`. You can't tell which shape a field has without
  checking — **probe one record first** (`("limit", "1")` on a raw query, or just look at
  the first row `get_search_paged` returns) before writing matching logic against it.
- **To match a bare-path-string field against something you only have an alias for**,
  resolve the alias first with `get_by_alias`, then compare against the resolved `@id` —
  don't try to filter server-side on an alias that lives inside an array field.

## Counting

Two ways, depending on whether you need the rows:

```python
# Just the count: one light page, read its "total" directly.
import urllib.parse, requests
import igvf_utils.utils as iuu

conn = reader._connection()
url = iuu.url_join([conn.igvf_mode.url, "search/?"])
query = urllib.parse.urlencode(params + [("limit", "0")])
resp = requests.get(url + query, auth=conn.auth, timeout=60, headers=iuu.REQUEST_HEADERS_JSON)
resp.raise_for_status()
total = resp.json()["total"]
```

```python
# Count AND need the rows for further analysis: just use get_search_paged and len() it.
# It already reconciled its row count against the Portal's own total, so len() here
# is the verified count, not a possibly-partial page.
rows = reader.get_search_paged(params)
count = len(rows)
```

Default to the second form whenever you're also going to tally or filter the rows — don't
fetch once for a count and again for the breakdown.

## Group-by / tally pattern

Fetch once with the minimal fields, then tally client-side with `collections.Counter`.
Never re-query per group — the whole point of fetching `@id` + the grouping field(s) in one
pass is that grouping is then free.

```python
import collections

params = [
    ("type", "PredictionSet"),
    ("software.title", "scE2G"),
    ("field", "@id"),
    ("field", "cell_type"),
    ("field", "submitted_by"),
]
rows = reader.get_search_paged(params)

# Unique cell types across ALL submitters:
cell_types_all = collections.Counter(
    (r.get("cell_type") or {}).get("@id") for r in rows
)
print(f"{len(cell_types_all)} unique cell types across {len(rows)} PredictionSets")

# Same, but only for your own submissions:
MY_USER_ID = "/users/326a5b8c-d3f2-4be5-8cc7-77f67398e812/"  # resolve your own via get_by_alias or a probe query first
mine = [r for r in rows if (r.get("submitted_by") or {}).get("@id") == MY_USER_ID]
cell_types_mine = collections.Counter(
    (r.get("cell_type") or {}).get("@id") for r in mine
)
print(f"{len(cell_types_mine)} unique cell types across your {len(mine)} PredictionSets")

# Full breakdown, everyone:
by_submitter = collections.Counter((r.get("submitted_by") or {}).get("@id") for r in rows)
for sid, n in by_submitter.most_common():
    print(sid, n)
```

If you don't already know your own `submitted_by.@id`, get it for free from the first row
`get_search_paged` returns for a query you know is yours (or do a one-off `get_by_alias` on
anything you know you submitted and read its `submitted_by` field) — don't guess the uuid.

## Worked example: does PredictionSet X carry file_format_specification doc Y, and by whom

This is the shape of "which objects reference alias Y, and how does that break down by
field Z" — resolve the alias to its real `@id` once, then filter and tally in memory:

```python
doc_id = reader.get_by_alias("jesse-engreitz:some_file_format_specs_md")["@id"]

params = [
    ("type", "PredictionSet"),
    ("software.title", "scE2G"),
    ("field", "@id"), ("field", "documents"),
    ("field", "status"), ("field", "submitted_by"),
]
rows = reader.get_search_paged(params)

matching = [r for r in rows if doc_id in (r.get("documents") or [])]
print(f"{len(matching)} / {len(rows)} use this doc")

by_status = collections.Counter(r.get("status") for r in matching)
by_submitter = collections.Counter((r.get("submitted_by") or {}).get("@id") for r in matching)
```

## References

- `references/query-recipes.md` — longer runnable scripts for common question shapes
  (resolve-then-filter, full paged tally, single-record inspection), copy-pasteable as a
  starting point for a new one-off script.
