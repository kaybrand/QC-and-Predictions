"""Cell Annotation / Cell Type metadata for Principal Pseudobulk Sets and the
(separately being built) reformatted E2G prediction tabular files -- both
need the same lookup, resolved from one live GET against data.igvf.org's
PseudobulkSet multireport (portal_client.PortalReader.get_multireport),
cached in state.db's cell_annotations table on a 24h TTL so a multi-cluster
run only ever issues the GET once, not once per cluster.

Grouping IS alias-based (reversed 2026-07-22, see below) -- a subsample
(an In-Vitro-System, MULTI-seq-tagged) is NOT a unique key for a pseudobulk:
the final cell-type/cluster annotation is made by downstream human analysis,
so one subsample routinely has MANY distinct primary pseudobulks, one per
resulting cluster (confirmed against real production data: IGVFSM8373MXSW
alone has 18 distinct primary pseudobulks, one per pancreatic differentiation
cluster). The real unique key is (subsample, cluster); today the alias
("{lab}:{dataset}-{cluster}-{subsample}") is the only field that encodes
which cluster a pseudobulk represents, so joining on it (matching the suffix
after the first ":", not assuming a specific lab prefix) is the only correct
mechanism for the ~500 primaries that follow it. The several-thousand
primaries coming soon are unlikely to follow this exact string format (a
separate, ongoing effort is on a calculable-alias mechanism for those) --
until that lands, primaries whose alias doesn't match any local
(dataset, cluster, subsample) candidate just don't join, logged, not fatal.

An EARLIER version of this module deliberately avoided alias parsing,
reasoning that a {subsample: (dataset, cluster)} reverse index built from
local subsample sets alone would be more robust. That reverse index is
WRONG: it assumed subsample uniquely identifies a pseudobulk, which the
multiplexed-cluster finding above disproves -- it would have silently kept
only the last-seen pseudobulk per subsample among the many that share it.

"Exactly 1 contributing sample" is not a data-quality heuristic -- it's the
actual filter for pipeline membership (distinguishes unified-processing
primaries, which QC-and-Predictions is part of, from every other
PseudobulkSet on the portal). One known uniformly-processed dataset doesn't
follow this yet -- a data-side fix, out of scope here.

A group's Cell Annotation is "locked" (shareable) if EITHER every
contributing primary's status == "released", OR some *principal*
pseudobulk row in the same multireport response already carries this exact
cell_annotation value (the portal auto-populates that field on a principal
upload) -- no local upload-ledger lookup needed for that second condition.

Primary/principal classification moved to pseudobulk_sets.classify (2026-08-17)
-- shared with portal_files.py's download discovery, so the one rule that
decides what a pseudobulk set IS cannot drift between the two consumers. See
that module for why file_set_type is NOT a usable discriminator.

Two different kinds of per-scope anomaly get two different treatments
(2026-07-24). A local subsample missing from its portal group, or
disagreement on (cl_id, term_id, term_name) across a scope's contributing
primaries, still skips just that one (dataset, cluster) scope entirely (same
try/except-and-log-and-continue discipline orchestrator.plan_table already
uses for enabled() failures) -- a skipped scope leaves its cache row
stale/absent, and get_metadata_for raises for it, same as any other
still-unresolved stub in this package. term_id/term_name/cl_id disagreement
specifically should never happen (they're all deterministic functions of the
same embedded cell_type object) and signals two genuinely different cell
types merged under one cluster name -- a real data problem to investigate,
not to paper over.

cell_annotation/cell_qualifier disagreement, in contrast, is resolved rather
than skipped: real production data confirms this is routine, benign
messiness (e.g. one contributing subsample's cell_annotation carries a
"derived from {other cell line}" cross-contamination label that the others
don't) -- picking the value from whichever contributing subsample has the
most cells in that cluster's own filtered barcode QC guide (the same
most-to-least-contributive subsample ordering already used for Prediction
Set's/Principal Pseudobulk Set's own `samples` field, see
subsamples.subsamples_by_frequency) resolves this cleanly instead of leaving
the scope permanently uncached.

2026-07-22 correction: EVERY primary/principal-pseudobulk row the
multireport GET returns is saved unconditionally, to
state.cell_metadata_primary_pseudobulks/cell_metadata_principal_pseudobulks
-- NOT gated on whether it happens to match a currently-configured local
cluster, or whether that cluster's whole group validates. The earlier
version only ever persisted the derived per-(dataset, cluster) view
(cell_annotations), so any primary pseudobulk not part of a fully-valid
local group was fetched into memory and then silently discarded --
defeating the entire point of a cache meant to hold what the portal
actually has. cell_annotations remains a stricter, derived, per-cluster
view (still useful for Principal Pseudobulk Set's consistency-checked
cell_type/cell_qualifier), built FROM the raw tables, not instead of them.
"""

import csv
import os
import sys
from datetime import datetime, timezone

from . import pseudobulk_sets, state, subsamples
from .context import Context

CACHE_TTL_HOURS = 24
SYNAPSE_PARENT_ID = "syn53469844"  # E2G Pillar Project's Synapse collaborative space
SHAREABLE_TSV_NAME = "cell_annotation_table.tsv"

# Confirmed against a real production call (2026-07-22): input_file_sets/samples/cell_type
# all come back as fully-embedded objects regardless of which sub-fields are requested here
# (no "@type" sub-field ever appears on input_file_sets entries -- classification instead
# uses their "@id" path prefix, see pseudobulk_sets.classify). limit=all is required: the multireport
# endpoint otherwise silently caps at a default page size (25, confirmed) instead of
# returning every PseudobulkSet.
_MULTIREPORT_QUERY = (
    "type=PseudobulkSet&status%21=deleted&limit=all"
    "&field=%40id&field=cell_annotation&field=aliases&field=cell_type"
    "&field=cell_type.term_name&field=cell_type.term_id&field=summary&field=cell_qualifier"
    "&field=input_file_sets&field=lab&field=samples&field=status"
)

_SHAREABLE_COLUMNS = [
    "dataset",
    "cluster",
    "cell_annotation",
    "cl_id",
    "term_id",
    "term_name",
    "cell_qualifier",
    "all_primary_released",
    "principal_uploaded",
    "principal_alias",
]


def log(msg):
    print(f"[cell_metadata] {msg}", file=sys.stderr)


def _now():
    return datetime.now(timezone.utc).isoformat()


def _is_stale(last_fetch):
    """last_fetch: state.latest_cell_annotation_fetch(conn)'s return value
    (None if the cache has never been populated)."""
    if not last_fetch:
        return True
    age = datetime.now(timezone.utc) - datetime.fromisoformat(last_fetch)
    return age.total_seconds() > CACHE_TTL_HOURS * 3600


def _cl_id_from_cell_type(cell_type):
    """cell_type comes back from the portal as a fully-embedded SampleTerm
    object (confirmed against a real production call, 2026-07-22), e.g.
    {"term_name": "macrophage", "term_id": "CL:0000235",
    "@id": "/sample-terms/CL_0000235/", ...} -- NOT a bare reference path/string.
    Uses "@id" (portal-style "CL_0000235", underscore) rather than "term_id"
    ("CL:0000235", colon/CURIE-style): principal_pseudobulk_set.py's _row()
    re-wraps this as f"/sample-terms/{cl_id}/", which only round-trips
    correctly with the "@id" form. If the E2G reformatted tabular files'
    "CL Term ID" column specifically wants the colon/CURIE form instead,
    that's cell_type["term_id"] on the raw cached row -- flag if so."""
    if not isinstance(cell_type, dict):
        return None
    ref = cell_type.get("@id") or ""
    return ref.strip("/").rsplit("/", 1)[-1] or None


def _term_name_from_cell_type(cell_type):
    """The human-readable name on the same embedded SampleTerm object
    _cl_id_from_cell_type reads (e.g. "macrophage" for CL:0000235) --
    requires field=cell_type.term_name in the multireport query string
    (confirmed against a real production call, 2026-07-23)."""
    if not isinstance(cell_type, dict):
        return None
    return cell_type.get("term_name") or None


def _term_id_from_cell_type(cell_type):
    """The portal's own CURIE-style "term_id" (e.g. "CL:0000235") -- NOT the
    same string as cl_id/_cl_id_from_cell_type's "@id"-derived form
    ("CL_0000235", underscore): the E2G tabular-file reformatting step wants
    this exact CURIE form for its "CL Term ID" column, so it's cached
    verbatim rather than derived by string-swapping cl_id."""
    if not isinstance(cell_type, dict):
        return None
    return cell_type.get("term_id") or None


def _local_subsamples(cluster_keys, cluster_configs):
    """{(dataset, cluster): [local subsample ids, most-to-least contributing]} --
    one entry per local scope. Ordered by descending cell count in that
    cluster's own filtered barcode QC guide (subsamples.subsamples_by_frequency
    -- same ordering already used for Prediction Set's/Principal Pseudobulk
    Set's own `samples` field) rather than an unordered set, so
    refresh_if_stale can resolve cell_annotation/cell_qualifier disagreement
    by picking the most-contributing subsample's own values."""
    by_scope = {}
    for dataset, cluster in cluster_keys:
        ctx = Context(dataset, cluster, None, cluster_configs[(dataset, cluster)], None, None, None)
        by_scope[(dataset, cluster)] = subsamples.subsamples_by_frequency(ctx)
    return by_scope


def _alias_suffix(alias):
    """The part after the first ":" -- "{dataset}-{cluster}-{subsample}" for
    the ~500 primaries that follow that convention, whatever it is for
    primaries that don't (in which case it simply won't match any local
    candidate suffix built the same way, below)."""
    return alias.split(":", 1)[-1] if alias and ":" in alias else alias


def refresh_if_stale(conn, reader, cluster_keys, cluster_configs):
    """One multireport GET per stale cache, shared across every (dataset,
    cluster) in this run -- not one GET per cluster. Callers invoked once
    per cluster (e.g. the pipeline-integrated hook, one cluster at a time --
    see orchestrator.py's own module docstring) still only hit the network
    once per 24h, as long as they share the same --state-db: the 2nd/3rd/...
    invocation's staleness check finds the 1st invocation's fetch still
    fresh and returns here without ever calling reader.get_multireport."""
    last_fetch = state.latest_cell_annotation_fetch(conn)
    if not _is_stale(last_fetch):
        log(f"cell_annotations cache still fresh (last fetched {last_fetch}) -- skipping multireport GET")
        return
    log(f"cell_annotations cache {'empty' if not last_fetch else 'stale'} -- issuing multireport GET")

    rows = reader.get_multireport(_MULTIREPORT_QUERY)
    now = _now()
    # Recorded unconditionally, before any per-scope validation below: the TTL is about
    # "did we hit the network recently," not "did every scope's data validate cleanly."
    # A round where every scope fails local-subset-of-portal validation must still count
    # as fetched, or the next invocation re-GETs immediately instead of waiting out the TTL.
    state.record_cell_annotations_fetch(conn, now)
    local_by_scope = _local_subsamples(cluster_keys, cluster_configs)

    principal_by_annotation = {}
    primary_by_alias_suffix = {}  # "{dataset}-{cluster}-{subsample}" -> its portal row
    saved_primary_count = 0
    skipped_no_alias = 0
    for row in rows:
        kind = pseudobulk_sets.classify(row)
        if kind == "principal":
            annotation = row.get("cell_annotation")
            aliases = row.get("aliases") or []
            alias = aliases[0] if aliases else None
            if annotation and alias:
                principal_by_annotation[annotation] = alias
                # Saved unconditionally -- this is the "already locked in" evidence,
                # independent of whether any local cluster currently references it.
                state.upsert_principal_pseudobulk(conn, alias, annotation, now)
            continue
        if kind != "primary":
            continue
        # samples entries are fully-embedded Sample/In-Vitro-System objects (confirmed
        # against a real production call, 2026-07-22), not bare accession strings.
        sample_accessions = [s.get("accession") for s in row.get("samples") or [] if isinstance(s, dict)]
        if len(sample_accessions) != 1:
            continue  # not part of the unified-processing pipeline this run
        subsample = sample_accessions[0]
        aliases = row.get("aliases") or []
        alias = aliases[0] if aliases else None
        if not alias:
            skipped_no_alias += 1
            continue  # can't identify this pseudobulk uniquely without an alias -- see module docstring
        # Saved unconditionally, keyed by alias, for EVERY primary pseudobulk the portal
        # returned -- not gated on whether it matches a currently-configured local
        # cluster. This is the actual point of the cache: a downstream consumer (e.g.
        # the E2G tabular file reformatting step) can look up a pseudobulk's raw Cell
        # Annotation/CL Term ID/Cell Qualifier here directly.
        state.upsert_primary_pseudobulk(
            conn,
            alias,
            subsample,
            row.get("cell_annotation"),
            _cl_id_from_cell_type(row.get("cell_type")),
            _term_id_from_cell_type(row.get("cell_type")),
            _term_name_from_cell_type(row.get("cell_type")),
            row.get("cell_qualifier"),
            row.get("status"),
            now,
        )
        saved_primary_count += 1
        primary_by_alias_suffix[_alias_suffix(alias)] = row

    log(
        f"saved {saved_primary_count} primary pseudobulk row(s) ({skipped_no_alias} skipped for lacking an alias) "
        f"and {len(principal_by_annotation)} principal pseudobulk row(s) to the raw cache"
    )

    matched_suffixes = set()
    for (dataset, cluster), local_subsamples in local_by_scope.items():
        # (subsample, cluster) -- not subsample alone -- is the real unique key (2026-07-22
        # finding: one MULTI-seq-tagged subsample routinely has many distinct pseudobulks,
        # one per downstream-annotated cluster). The alias is the only field that currently
        # encodes cluster identity, so match candidate suffixes built from OUR OWN known
        # (dataset, cluster, subsample) triples against it, rather than parsing the alias's
        # ambiguous hyphen-separated segments blind.
        candidate_suffixes = {f"{dataset}-{cluster}-{s}": s for s in local_subsamples}
        matched_suffixes |= set(candidate_suffixes)
        missing_locally = [s for suffix, s in candidate_suffixes.items() if suffix not in primary_by_alias_suffix]
        if missing_locally:
            log(
                f"{dataset}/{cluster}: local subsample(s) {sorted(missing_locally)} have no corresponding "
                f"primary pseudobulk alias matching \"{{lab}}:{dataset}-{{cluster}}-{{subsample}}\" on the "
                f"portal -- skipping this scope's cache until resolved (expected for primaries uploaded "
                f"under a different alias convention -- see module docstring)"
            )
            continue
        scope_rows = [primary_by_alias_suffix[suffix] for suffix in candidate_suffixes]

        # term_id/term_name/cl_id all come off the SAME embedded cell_type object and
        # should always agree across every contributing subsample -- unlike
        # cell_annotation (which routinely carries legitimate "derived from X"
        # cross-contamination text), disagreement here means two genuinely different
        # cell types got merged under one cluster name. That's a real data problem to
        # flag and investigate, not to silently paper over -- skip this scope's cache,
        # same as a missing-subsample scope.
        term_triples = {
            (
                _cl_id_from_cell_type(r.get("cell_type")),
                _term_id_from_cell_type(r.get("cell_type")),
                _term_name_from_cell_type(r.get("cell_type")),
            )
            for r in scope_rows
        }
        if len(term_triples) != 1:
            log(
                f"WARNING {dataset}/{cluster}: {len(term_triples)} distinct (cl_id, term_id, term_name) "
                f"triples across its primary pseudobulks -- this should never happen (these are supposed "
                f"to always agree, unlike cell_annotation) -- skipping this scope's cache until investigated"
            )
            continue
        cl_id, term_id, term_name = next(iter(term_triples))

        # cell_annotation/cell_qualifier legitimately vary -- resolve via the
        # most-contributing subsample (local_subsamples is ordered most-to-least
        # contributing, see _local_subsamples) rather than requiring unanimous
        # agreement across every contributing subsample.
        annotation_qualifier_pairs = {(r.get("cell_annotation"), r.get("cell_qualifier")) for r in scope_rows}
        winning_subsample = local_subsamples[0]
        winning_row = primary_by_alias_suffix[f"{dataset}-{cluster}-{winning_subsample}"]
        if len(annotation_qualifier_pairs) != 1:
            log(
                f"{dataset}/{cluster}: {len(annotation_qualifier_pairs)} distinct Cell Annotation/Cell "
                f"Qualifier pairs across its primary pseudobulks -- resolved using the most-contributing "
                f"subsample ({winning_subsample})"
            )
        cell_annotation = winning_row.get("cell_annotation")
        cell_qualifier = winning_row.get("cell_qualifier")

        all_primary_released = all(r.get("status") == "released" for r in scope_rows)
        principal_alias = principal_by_annotation.get(cell_annotation)

        state.upsert_cell_annotation(
            conn,
            dataset,
            cluster,
            cell_annotation,
            cl_id,
            term_id,
            term_name,
            cell_qualifier,
            ",".join(sorted(local_subsamples)),
            all_primary_released,
            principal_alias is not None,
            principal_alias,
            now,
        )

    unmatched = set(primary_by_alias_suffix) - matched_suffixes
    if unmatched:
        log(
            f"{len(unmatched)} portal primary pseudobulk(s) whose alias didn't match any local "
            f"(dataset, cluster, subsample) candidate this run (still saved to the raw cache above -- "
            f"just not part of any currently-configured local cluster's derived group)"
        )


def get_metadata_for(ctx):
    """Should return {"cl_id": ..., "term_id": ..., "term_name": ...,
    "cell_qualifier": ...} for (ctx.dataset, ctx.cluster) -- raises if that
    scope's cache row is missing (never successfully cached yet -- see
    refresh_if_stale's per-scope skip conditions above)."""
    row = state.get_cell_annotation(ctx.conn, ctx.dataset, ctx.cluster)
    if row is None:
        raise ValueError(
            f"{ctx.dataset}/{ctx.cluster}: no cached Cell Annotation metadata -- either the multireport "
            "cache hasn't been refreshed yet, or this scope failed its consistency/subset validation "
            "(see cell_metadata.refresh_if_stale logs)"
        )
    return {
        "cl_id": row["cl_id"],
        "term_id": row["term_id"],
        "term_name": row["term_name"],
        "cell_qualifier": row["cell_qualifier"],
    }


def build_shareable_rows(conn):
    """Only rows that won't change further: every contributing primary
    already released, or a principal pseudobulk already uploaded using this
    exact Cell Annotation."""
    return [r for r in state.all_cell_annotations(conn) if r["all_primary_released"] or r["principal_uploaded"]]


def _write_shareable_tsv(rows, path):
    dirname = os.path.dirname(path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(_SHAREABLE_COLUMNS)
        for row in rows:
            writer.writerow([row.get(c, "") for c in _SHAREABLE_COLUMNS])
    return path


def push_to_synapse(rows, manifest_dir=None):
    """Writes the shareable-rows TSV and stores it as a child of
    SYNAPSE_PARENT_ID -- synapseclient auto-versions a repeated store() to
    the same (name, parent) File, so "keep it updated" needs no extra
    bookkeeping here."""
    if not rows:
        return None
    import synapseclient

    path = os.path.join(manifest_dir or ".", SHAREABLE_TSV_NAME)
    _write_shareable_tsv(rows, path)
    syn = synapseclient.login()
    entity = synapseclient.File(path, parent=SYNAPSE_PARENT_ID, name=SHAREABLE_TSV_NAME)
    return syn.store(entity)
