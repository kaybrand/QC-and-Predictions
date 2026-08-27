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

2026-08-24 split: refresh_if_stale is now a thin wrapper over fetch_if_stale
(network, global-TTL-gated, wholesale) and derive_scopes (pure local, no TTL
gate, runs every time). Conflating them was a real bug -- the TTL is one GLOBAL
row, so warming for dataset A made it fresh and dataset B's later invocation
returned early BEFORE ever reaching the derivation loop, leaving B's clusters
permanently uncached even though the raw rows they needed were already cached.
Split apart, a dataset whose config is written weeks after the fetch still
derives its own cell_annotations offline, and no dataset depends on another
running first. See those two functions for details.
"""

import csv
import hashlib
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
    """({(dataset, cluster): [local subsample ids, most-to-least contributing]},
    {(dataset, cluster): reason}) -- the first dict has one entry per readable
    local scope, the second names the scopes whose QC guide couldn't be read at
    all. Ordered by descending cell count in that cluster's own filtered barcode
    QC guide (subsamples.subsamples_by_frequency -- same ordering already used
    for Prediction Set's/Principal Pseudobulk Set's own `samples` field) rather
    than an unordered set, so derive_scopes can resolve cell_annotation/
    cell_qualifier disagreement by picking the most-contributing subsample's own
    values.

    An unreadable guide is REPORTED, not raised: this runs across every cluster
    in a pipeline config, and one missing/corrupt guide must not stop every
    other scope from resolving -- the same try/except-and-log-and-continue
    discipline the per-scope anomaly handling in derive_scopes already uses."""
    by_scope = {}
    unreadable = {}
    for dataset, cluster in cluster_keys:
        ctx = Context(dataset, cluster, None, cluster_configs[(dataset, cluster)], None, None, None)
        try:
            by_scope[(dataset, cluster)] = subsamples.subsamples_by_frequency(ctx)
        except (OSError, KeyError) as e:
            # OSError: guide absent/unreadable. KeyError: present but missing the
            # "subsample" column, i.e. not actually a QC guide.
            unreadable[(dataset, cluster)] = f"unreadable_qc_guide: {type(e).__name__}: {e}"
    return by_scope, unreadable


def _merge_sources(cluster_cfg, cluster):
    """The annotation name(s) a cluster's primary pseudobulks are aliased under,
    most clusters having exactly one: their own name.

    `pseudobulk_annotation` is this repo's established way of saying "this cluster
    is built from these source annotations", comma-separated when merged (igvf18's
    mcf7 = "mcf7_1,mcf7_2").

    ONLY the comma case is honoured here, and that restriction is load-bearing.
    `pseudobulk_annotation` is NOT always the cluster name for unmerged clusters --
    the ATAC-only variants are exactly the counterexample (igvf18's
    hudep_d7_ATAC_only carries `pseudobulk_annotation: hudep_d7`). Honouring it for
    them would let their `{dataset}-hudep_d7-{subsample}` primaries suddenly match,
    flipping three long-unresolved clusters to resolved -- which makes them
    reformat-eligible and generates scATAC prediction files that are deliberately
    not being shared this round. Those clusters have their own override for
    annotation lookup (`cell_annotation_key`, see cell_annotations.py); alias
    matching is not the place to reinterpret them.

    So: merged -> the constituent names; everything else -> [cluster], byte-identical
    to the behaviour before merged-cluster support existed.
    """
    raw = str(cluster_cfg.get("pseudobulk_annotation") or "")
    if "," not in raw:
        return [cluster]
    names = [part.strip() for part in raw.split(",") if part.strip()]
    return names or [cluster]


def _alias_suffix(alias):
    """The part after the first ":" -- "{dataset}-{cluster}-{subsample}" for
    the ~500 primaries that follow that convention, whatever it is for
    primaries that don't (in which case it simply won't match any local
    candidate suffix built the same way, below)."""
    return alias.split(":", 1)[-1] if alias and ":" in alias else alias


def fetch_if_stale(conn, reader):
    """The NETWORK half. One wholesale multireport GET per stale cache, saving
    EVERY primary/principal pseudobulk row the portal returns to the raw cache
    tables. Returns True if a GET actually happened, False if the cache was
    still fresh.

    Deliberately knows nothing about which clusters anyone cares about: the
    GET's scope is "every PseudobulkSet on the portal", and all of it is saved
    unconditionally (see the module docstring's 2026-07-22 correction). Turning
    that raw cache into the per-(dataset, cluster) `cell_annotations` view is
    derive_scopes' job, and needs no network at all.

    Split out of the former refresh_if_stale (2026-08-24) because the two halves
    have genuinely different staleness semantics, and conflating them was a real
    bug. The TTL is a single GLOBAL row, so warming for one dataset made it
    fresh, and every LATER dataset's invocation returned early -- before ever
    reaching the derivation loop -- leaving those clusters permanently uncached
    even though the raw data they needed was already sitting in the cache.
    Confirmed against the live DB (2026-08-24): 1039 primary rows cached, 62 of
    them igvf0's, while igvf0 had zero `cell_annotations` rows.
    """
    last_fetch = state.latest_cell_annotation_fetch(conn)
    if not _is_stale(last_fetch):
        log(f"cell_annotations cache still fresh (last fetched {last_fetch}) -- skipping multireport GET")
        return False
    log(f"cell_annotations cache {'empty' if not last_fetch else 'stale'} -- issuing multireport GET")

    rows = reader.get_multireport(_MULTIREPORT_QUERY)
    now = _now()
    # Recorded unconditionally, before any per-scope validation: the TTL is about
    # "did we hit the network recently," not "did every scope's data validate cleanly."
    # A round where every scope fails local-subset-of-portal validation must still count
    # as fetched, or the next invocation re-GETs immediately instead of waiting out the TTL.
    #
    # Recorded AFTER get_multireport returns, deliberately: a failed GET must not
    # poison the TTL, or a concurrent second driver would see a phantom "fresh"
    # cache and skip a fetch that never actually succeeded.
    state.record_cell_annotations_fetch(conn, now)

    saved_primary_count = 0
    saved_principal_count = 0
    skipped_no_alias = 0
    for row in rows:
        kind = pseudobulk_sets.classify(row)
        if kind == "principal":
            annotation = row.get("cell_annotation")
            aliases = row.get("aliases") or []
            alias = aliases[0] if aliases else None
            if annotation and alias:
                # Saved unconditionally -- this is the "already locked in" evidence,
                # independent of whether any local cluster currently references it.
                state.upsert_principal_pseudobulk(conn, alias, annotation, now)
                saved_principal_count += 1
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

    log(
        f"saved {saved_primary_count} primary pseudobulk row(s) ({skipped_no_alias} skipped for lacking an alias) "
        f"and {saved_principal_count} principal pseudobulk row(s) to the raw cache"
    )
    return True


def derive_scopes(conn, cluster_keys, cluster_configs):
    """The LOCAL half. Builds the per-(dataset, cluster) `cell_annotations` view
    from the already-cached raw tables plus each cluster's own QC guide. Makes
    NO network call and has NO TTL gate -- it runs every time, so a dataset
    whose config was written long after the fetch still gets its own rows from
    whatever the raw cache already holds.

    Returns [{dataset, cluster, resolved, reason, cell_annotation}] -- one entry
    per requested scope, so a caller can report exactly which clusters resolved
    and why the rest didn't, instead of a resolved scope and an unresolved one
    being indistinguishable.

    Reads the FLATTENED raw rows (state.all_primary_pseudobulks /
    all_principal_pseudobulks) rather than live portal response objects, so
    cl_id/term_id/term_name come straight off the cached columns -- the
    _*_from_cell_type() unpacking happens once, in fetch_if_stale, on the way in.

    Every per-scope validation rule from the former refresh_if_stale is preserved
    exactly: a local subsample with no matching portal alias, or disagreement on
    (cl_id, term_id, term_name), skips the scope; cell_annotation/cell_qualifier
    disagreement is resolved via the most-contributing subsample. See the module
    docstring for why those two anomalies get different treatment.
    """
    now = _now()
    local_by_scope, unreadable = _local_subsamples(cluster_keys, cluster_configs)
    statuses = []
    for (dataset, cluster), reason in sorted(unreadable.items()):
        log(f"{dataset}/{cluster}: {reason} -- cannot derive this scope's Cell Annotation")
        statuses.append(
            {"dataset": dataset, "cluster": cluster, "resolved": False, "reason": reason, "cell_annotation": None}
        )

    primary_by_alias_suffix = {
        _alias_suffix(raw["alias"]): raw for raw in state.all_primary_pseudobulks(conn)
    }
    principal_by_annotation = {
        raw["cell_annotation"]: raw["alias"]
        for raw in state.all_principal_pseudobulks(conn)
        if raw["cell_annotation"] and raw["alias"]
    }
    if not primary_by_alias_suffix:
        log("raw primary-pseudobulk cache is empty -- nothing to derive from (fetch_if_stale first)")

    matched_suffixes = set()
    for (dataset, cluster), local_subsamples in local_by_scope.items():
        # (subsample, cluster) -- not subsample alone -- is the real unique key (2026-07-22
        # finding: one MULTI-seq-tagged subsample routinely has many distinct pseudobulks,
        # one per downstream-annotated cluster). The alias is the only field that currently
        # encodes cluster identity, so match candidate suffixes built from OUR OWN known
        # (dataset, cluster, subsample) triples against it, rather than parsing the alias's
        # ambiguous hyphen-separated segments blind.
        #
        # A MERGED cluster has no primary of its own. igvf18's mcf7 is
        # `pseudobulk_annotation: mcf7_1,mcf7_2` and the portal only ever saw
        # `igvf18-mcf7_1-{subsample}` / `igvf18-mcf7_2-{subsample}` -- nothing is
        # aliased `igvf18-mcf7-...`, so matching on the cluster name alone reports
        # every subsample as no_matching_primary_alias. Match on the CONSTITUENT
        # annotation names instead, and accept a subsample that resolves under any
        # one of them (mcf7_1 carries all 7; mcf7_2 only 4).
        #
        # Gated on the cluster actually being merged (see _merge_sources) so that
        # non-merged clusters keep the exact previous code path. That gate is not
        # mere caution: `pseudobulk_annotation` is NOT reliably equal to `cluster`
        # for unmerged clusters, so matching on it unconditionally would silently
        # resolve the ATAC-only variants too. See _merge_sources' docstring.
        source_names = _merge_sources(cluster_configs[(dataset, cluster)], cluster)
        candidate_suffixes = {}   # suffix -> subsample, only for suffixes that exist
        rows_by_subsample = {}    # subsample -> its matched row(s), in source_names order
        missing_locally = []
        for s in local_subsamples:
            hits = [f"{dataset}-{name}-{s}" for name in source_names
                    if f"{dataset}-{name}-{s}" in primary_by_alias_suffix]
            if hits:
                candidate_suffixes.update({suffix: s for suffix in hits})
                rows_by_subsample[s] = [primary_by_alias_suffix[suffix] for suffix in hits]
            else:
                missing_locally.append(s)
            # Report every suffix we CONSIDERED as matched, not just the ones that
            # hit, so the trailing unmatched-primaries audit isn't polluted by the
            # constituent names we deliberately probed.
            matched_suffixes |= {f"{dataset}-{name}-{s}" for name in source_names}
        if missing_locally:
            log(
                f"{dataset}/{cluster}: local subsample(s) {sorted(missing_locally)} have no corresponding "
                f"primary pseudobulk alias matching \"{{lab}}:{dataset}-{{{'|'.join(source_names)}}}-{{subsample}}\" "
                f"on the portal -- skipping this scope's cache until resolved (expected for primaries uploaded "
                f"under a different alias convention -- see module docstring)"
            )
            statuses.append(
                {
                    "dataset": dataset,
                    "cluster": cluster,
                    "resolved": False,
                    "reason": f"no_matching_primary_alias: {','.join(sorted(missing_locally))}",
                    "cell_annotation": None,
                }
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
        term_triples = {(r["cl_id"], r["term_id"], r["term_name"]) for r in scope_rows}
        if len(term_triples) != 1:
            log(
                f"WARNING {dataset}/{cluster}: {len(term_triples)} distinct (cl_id, term_id, term_name) "
                f"triples across its primary pseudobulks -- this should never happen (these are supposed "
                f"to always agree, unlike cell_annotation) -- skipping this scope's cache until investigated"
            )
            statuses.append(
                {
                    "dataset": dataset,
                    "cluster": cluster,
                    "resolved": False,
                    "reason": f"cell_type_disagreement: {len(term_triples)} distinct (cl_id, term_id, term_name)",
                    "cell_annotation": None,
                }
            )
            continue
        cl_id, term_id, term_name = next(iter(term_triples))

        # cell_annotation/cell_qualifier legitimately vary -- resolve via the
        # most-contributing subsample (local_subsamples is ordered most-to-least
        # contributing, see _local_subsamples) rather than requiring unanimous
        # agreement across every contributing subsample.
        annotation_qualifier_pairs = {(r["cell_annotation"], r["cell_qualifier"]) for r in scope_rows}
        if len(source_names) > 1:
            # MERGED cluster: the per-constituent qualifier describes the CONSTITUENT,
            # not the merge, so it cannot be carried onto the merged object. igvf18's
            # mcf7 merges a MED13L+ half and a VMP1+ half; the most-contributing-
            # subsample tie-break below would stamp the whole thing "MED13L+", which is
            # affirmatively wrong for the half that isn't. Both halves agree on the
            # cell type itself (the term_triples check above already proved that), so
            # the merge's honest annotation is the bare term_name with NO qualifier --
            # which is also exactly the shape its unqualified siblings take
            # (igvf18/hct116 -> 'HCT116', qualifier NULL).
            cell_annotation = term_name
            cell_qualifier = None
            log(
                f"{dataset}/{cluster}: merged cluster over {source_names} with "
                f"{len(annotation_qualifier_pairs)} distinct Cell Annotation/Cell Qualifier pair(s) -- "
                f"using the shared cell type {term_name!r} and dropping the per-constituent qualifier "
                f"({sorted({r['cell_qualifier'] for r in scope_rows if r['cell_qualifier']})})"
            )
        else:
            winning_subsample = local_subsamples[0]
            # Indexed off the rows we actually matched, never re-derived from the
            # cluster name: those two agree only while `pseudobulk_annotation` ==
            # `cluster`, which the ATAC-only variants disprove, and reconstructing
            # the key here raised KeyError mid-loop and took down every remaining
            # scope's derivation with it.
            winning_row = rows_by_subsample[winning_subsample][0]
            if len(annotation_qualifier_pairs) != 1:
                log(
                    f"{dataset}/{cluster}: {len(annotation_qualifier_pairs)} distinct Cell Annotation/Cell "
                    f"Qualifier pairs across its primary pseudobulks -- resolved using the most-contributing "
                    f"subsample ({winning_subsample})"
                )
            cell_annotation = winning_row["cell_annotation"]
            cell_qualifier = winning_row["cell_qualifier"]

        all_primary_released = all(r["status"] == "released" for r in scope_rows)
        principal_alias = principal_by_annotation.get(cell_annotation)

        # cell_annotations declares cell_annotation/cl_id/term_id/term_name NOT NULL,
        # so a portal row missing any of them would raise IntegrityError mid-loop and
        # take down every remaining scope's derivation with it. Check first and report
        # this scope instead -- same per-scope-skip discipline as the anomalies above.
        blank = [
            name
            for name, value in (
                ("cell_annotation", cell_annotation), ("cl_id", cl_id),
                ("term_id", term_id), ("term_name", term_name),
            )
            if not value
        ]
        if blank:
            log(f"{dataset}/{cluster}: portal primaries have no {'/'.join(blank)} -- skipping this scope's cache")
            statuses.append(
                {
                    "dataset": dataset, "cluster": cluster, "resolved": False,
                    "reason": f"blank_required_field: {','.join(blank)}", "cell_annotation": None,
                }
            )
            continue

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
        statuses.append(
            {"dataset": dataset, "cluster": cluster, "resolved": True, "reason": "ok", "cell_annotation": cell_annotation}
        )

    unmatched = set(primary_by_alias_suffix) - matched_suffixes
    if unmatched:
        log(
            f"{len(unmatched)} cached primary pseudobulk(s) whose alias didn't match any local "
            f"(dataset, cluster, subsample) candidate among the {len(cluster_keys)} scope(s) requested "
            f"here -- still in the raw cache, just not part of any of THESE clusters' derived groups "
            f"(another dataset's own derive_scopes call will pick up the ones that belong to it)"
        )
    resolved = sum(1 for s in statuses if s["resolved"])
    log(f"derived {resolved}/{len(statuses)} requested scope(s) into cell_annotations")
    return statuses


def refresh_if_stale(conn, reader, cluster_keys, cluster_configs):
    """Back-compat entry point, unchanged signature: fetch (network, TTL-gated)
    then derive (local, always). Callers invoked once per cluster still only hit
    the network once per 24h, as long as they share the same --state-db -- but
    unlike before, each one now derives its OWN scopes from the raw cache
    regardless of whose fetch populated it. See fetch_if_stale for the bug this
    split fixes."""
    fetch_if_stale(conn, reader)
    return derive_scopes(conn, cluster_keys, cluster_configs)


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


SYNAPSE_DIGEST_KEY = "synapse_cell_annotation_digest"


def _content_digest(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def push_to_synapse(rows, manifest_dir=None, conn=None, force=False):
    """Write the shareable-rows TSV, and store it as a child of
    SYNAPSE_PARENT_ID **only if its content actually changed**.

    Returns a status dict the caller can report:
      {"status": "pushed"|"unchanged"|"empty", "rows": n, "digest": ..., "path": ...}
    plus "entity_id"/"version" when a store() really happened.

    Change detection is by sha256 of the written file, compared against
    state.sync_state[SYNAPSE_DIGEST_KEY]. Synapse's own store() auto-versions a
    repeated upload to the same (name, parent), so without this gate every pass
    of a multi-pass upload session added a fresh, identical version -- and since
    a session is six passes for one cluster, that is five versions of noise for
    no content change.

    The TSV is written locally either way: it is cheap, and it keeps
    manifest_dir's copy current so the digest reflects the file on disk.

    `conn` is required for the gate. Without it (or with force=True) this
    unconditionally pushes, which is the old behaviour -- keep that available
    deliberately, e.g. to repair a Synapse-side edit the digest cannot see.
    """
    if not rows:
        return {"status": "empty", "rows": 0}

    path = os.path.join(manifest_dir or ".", SHAREABLE_TSV_NAME)
    _write_shareable_tsv(rows, path)
    digest = _content_digest(path)

    if conn is not None and not force:
        if state.get_sync_state(conn, SYNAPSE_DIGEST_KEY) == digest:
            return {"status": "unchanged", "rows": len(rows), "digest": digest, "path": path}

    import synapseclient

    syn = synapseclient.login()
    entity = synapseclient.File(path, parent=SYNAPSE_PARENT_ID, name=SHAREABLE_TSV_NAME)
    stored = syn.store(entity)

    if conn is not None:
        state.set_sync_state(conn, SYNAPSE_DIGEST_KEY, digest)

    return {
        "status": "pushed",
        "rows": len(rows),
        "digest": digest,
        "path": path,
        "entity_id": stored.get("id") if hasattr(stored, "get") else None,
        "version": stored.get("versionNumber") if hasattr(stored, "get") else None,
    }
