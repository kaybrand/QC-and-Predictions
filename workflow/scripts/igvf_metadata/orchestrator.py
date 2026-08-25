"""Core upload loop, shared by the pipeline-integrated hook and the
standalone scanner -- both just assemble a different list of (dataset,
cluster) scope-keys and hand it to run(). See registry.py and state.py for
the data model, and portal_client.py for why real submission goes through
iu_register.py rather than direct API calls.

Two phases per table, every run:
  1. plan_table(): for every enabled row not already uploaded-and-unchanged,
     build+validate its payload, hash it, and classify it as "post" (no
     existing portal record for this alias) or "patch" (exists, and the
     hash changed) via a live, read-only alias lookup. This is the ONLY
     network contact in "preview" mode.
  2. Write one post.tsv/patch.tsv per (dataset, table, variant, round) group
     under manifest_dir/<dataset>/, filename-prefixed with that round number
     -- see "Rounds" below. These are deliberately EPHEMERAL: only rows not
     already uploaded-and-unchanged ever appear, so a table's post/patch
     file shrinks down to just what's still outstanding as pieces get
     uploaded (e.g. an already-uploaded cluster's rows don't clutter review
     of newly-added ones) -- confirmed 2026-08-05 as a feature to keep, not
     a gap to fix.

Per-(object_type, dataset) accumulator (added 2026-08-05, manifest_dir/
<dataset>/<object_type>.tsv): a SEPARATE, durable record of every alias ever
confirmed live, one row per record_id, upserted every run (see
portal_client.merge_write_tsv) -- never shrinks, never dropped, values
updated in place if they change. Grouped by the portal's own object_type
(iu_register.py's --profile_id), not our internal table_name, since several
of our table modules share one object_type (e.g. filtered_barcode_list,
filtered_atac_fragment_file, and every prediction_tabular_files variant are
all "tabular_file") -- this is what actually makes a future bulk field edit
across a whole object type a single edit-and-resubmit against one file,
matching iu_register.py's own one-profile-per-invocation constraint. Fed
from two places: plan_table's "unchanged" branch (already live, no live GET
needed -- this is also how the accumulator backfills anything uploaded
before this feature existed, automatically, on the next ordinary run) and
_verify_and_record's return value (rows a real --mode upload pass just
confirmed). Never itself passed to invoke_register -- it's a reference, not
a submission.

depends_on / rounds (redesigned 2026-08-05 -- see git history for the prior
design this replaced): a row's dependencies name other (table_name,
variant_name) pairs. Two DIFFERENT things used to be conflated under one
gate ("is the dependency status='uploaded' in the local ledger"):
  (a) external dependencies we don't control the timing of (the
      Kundaje-lab primary pseudobulks) -- waiting for those to actually
      exist made sense, since we can't predict when they land.
  (b) our OWN dependency chain (QC_documents -> principal_pseudobulk_set ->
      filtered_barcode_list/prediction_set -> ...) -- every alias in this
      chain is a deterministic formula we can compute right now, regardless
      of what's actually live on the portal yet, since WE fully control the
      order things get uploaded in.
Blocking (b) on the same live-status gate as (a) meant nothing past round 1
could even be *built* to look at, forcing a slow "post one layer, wait,
regenerate, next layer unblocks" loop -- fine for testing the pipeline,
bad for actually running a real upload session where you want to review
everything in one pass and then upload round by round.

Now: depends_on's live-status gate still applies to "upload" mode (skip
calling iu_register.py for a row whose dependency isn't confirmed live yet
-- never submit a dangling cross-reference for real). In "preview"/
"validate" mode, that gate is informational only (logged as PENDING, not
skipped) -- every row's payload gets built and written regardless, tagged
with its `round` (1 + max(round of each dependency), via _compute_round,
memoized per (dataset, cluster, model-or-none, table, variant) since the
same dependency is reused across many rows). Output filenames are
round-prefixed (round{N}_{table}[_{variant}]_{post,patch}.tsv) so sorting
the manifest_dir by filename IS the upload order -- ties (same round) get
the same number, by construction.

Row ordering across tables/variants otherwise falls out of just
re-running: a row still deferred on a dependency in "upload" mode is
picked up automatically next time run() is called, whether that's the
pipeline hook firing for a later cluster or the scanner's periodic
reconciliation.

cell_metadata.refresh_if_stale() (2026-07-21) runs once at the top, before
any table -- a single PseudobulkSet-multireport GET (24h TTL) that
Principal Pseudobulk Set's payload-building depends on. Grouped/cached
against EVERY cluster in the pipeline config (all_cluster_configs), not just
this invocation's cluster_keys -- otherwise a later invocation for a
different cluster subset, within the same TTL window, would find nothing
cached for its own clusters. On a real upload pass (mode="upload"),
cell_metadata.push_to_synapse() runs once at the end, re-publishing the
"locked" subset of that cache to the shared Synapse
space -- see cell_metadata.py for both.
"""

import csv
import os
import sys
from datetime import datetime, timezone

from . import cell_metadata, portal_client, registry, state
from .context import Context

# Per-dataset coverage artifact (manifest_dir/<dataset>/). Deliberately NOT one of
# the round{N}_* submission files and never handed to iu_register.py -- it is the
# audit trail for "what did this run decide about every row it could have built",
# which is what makes a gap actionable instead of invisible.
MANIFEST_COVERAGE_NAME = "manifest_coverage.tsv"
MANIFEST_COVERAGE_HEADER = ["dataset", "cluster", "model", "table", "variant", "outcome", "reason"]

# Outcomes that mean a row we EXPECTED is not in the manifest. Everything else is
# either a real row (planned-post/planned-patch/unchanged) or a deliberate,
# expected exclusion (skipped-family-gated) or a normal ordering wait (deferred).
MANIFEST_GAP_OUTCOMES = frozenset({"skipped-missing-file", "invalid", "enabled-check-failed"})


def log(msg):
    print(f"[igvf_metadata] {msg}", file=sys.stderr)


def _group_by_dataset(rows):
    grouped = {}
    for row in rows:
        grouped.setdefault(row["dataset"], []).append(row)
    return grouped


def _retire_superseded_sibling(manifest_dir, dataset, table_name, suffix, round_num, kind):
    """Keep at most one of {post, patch} per (dataset, table, variant, round),
    without ever leaving zero.

    Called immediately AFTER the replacement has been written, so there is no
    code path here that deletes a file whose successor doesn't already exist on
    disk. That matters because in a preview-only workflow these round*.tsv files
    are the ONLY on-disk record of a row's field values and submitted_file_name --
    the durable {object_type}.tsv accumulator only ever receives rows already
    confirmed live (plan_table's "unchanged" branch and _verify_and_record), so it
    stays empty until something is really uploaded through this tool. Deleting the
    last file would destroy the record.

    Deliberately ONE-DIRECTIONAL. Aliases are stable per object type, so the real
    progression is: absent -> post -> live -> patch, terminally.
      - Writing a PATCH retires a stale POST: those objects are live, so the POST
        can never be the right instruction again.
      - Writing a POST does NOT retire a PATCH. That direction is backwards, and
        plan_table's own lookup can produce it spuriously -- it calls
        get_by_alias() WITHOUT database=True (orchestrator.py's plan_table), i.e.
        the Elasticsearch-backed read that _verify_and_record's own comment
        documents as able to "falsely report 'not found'". Silently dropping a
        good patch file on a transient false negative would lose real work, so
        this warns and leaves both for a human instead.
    """
    sibling_kind = "post" if kind == "patch" else "patch"
    sibling = os.path.join(manifest_dir, dataset, f"round{round_num}_{table_name}{suffix}_{sibling_kind}.tsv")
    if not os.path.exists(sibling):
        return None

    if kind == "patch":
        os.remove(sibling)
        log(f"  retired superseded {dataset}/{os.path.basename(sibling)} (its rows are live; patch written)")
        return sibling

    log(
        f"  WARNING {dataset}/{os.path.basename(sibling)} exists but this round planned a POST -- "
        "either those objects were deleted portal-side, or the alias lookup lagged and reported them "
        "absent. Left BOTH files in place; confirm which is right before submitting either."
    )
    return None


def write_coverage_tsv(path, rows):
    """Full deterministic overwrite, sorted -- unlike the round{N}_* submission
    files this is a complete snapshot of one run's decisions, so it should never
    accumulate rows from an earlier run with different inputs."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COVERAGE_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["cluster"], r["table"], r["variant"], r["model"])):
            writer.writerow(row)
    return path


def _now():
    return datetime.now(timezone.utc).isoformat()


def build_payload(table, variant, ctx, item_alias):
    payload = {"aliases": [item_alias], "award": ctx.igvf.award, "lab": ctx.igvf.lab}
    payload.update(table.constant_fields)
    payload.update(table.scope_fields(ctx))
    payload.update(variant.build_row(ctx))
    return payload


def _is_missing(value):
    """Distinguishes "absent" (None, "", [], {}) from "falsy but a real
    value" (False, 0) -- a bare `not value` check would wrongly flag a
    legitimate `False` (e.g. controlled_access) as missing."""
    if value is None:
        return True
    if isinstance(value, (str, list, tuple, dict)) and len(value) == 0:
        return True
    return False


def validate_row(table, variant_name, payload):
    missing = [c for c in table.required_columns if _is_missing(payload.get(c))]
    if missing:
        raise ValueError(f"{table.name}/{variant_name}: missing required column(s) {missing}")


def _iter_scopes(table, cluster_keys, cluster_configs, igvf_cfg):
    for dataset, cluster in cluster_keys:
        cluster_cfg = cluster_configs[(dataset, cluster)]
        if table.scope == "cluster_model":
            # Local import, not top-level: avoids a circular import (tables/*
            # import registry, which orchestrator also imports), same pattern
            # refs.py already uses. cluster_cfg["models"] reflects what scE2G
            # actually ran for this cluster -- igvf_cfg.enabled_families is the
            # separate, IGVF-specific gate on which of those families should
            # actually be uploaded this run (Multiome-only in 2026).
            from .tables.prediction_tabular_files import family

            models = [m for m in cluster_cfg["models"] if family(m) in igvf_cfg.enabled_families]
        else:
            models = [None]
        for model in models:
            yield dataset, cluster, model, cluster_cfg


def _dependency_model(dep_table_name, model_key):
    """The dependency's OWN scope decides its model_key, not the dependent
    variant's -- a cluster_model-scoped table (e.g. prediction_set)
    depending on a cluster-scoped one (e.g. principal_pseudobulk_set) must
    look/compute that row under model="" (how it's actually scoped), never
    under the dependent's real model name, or the lookup silently never
    matches (confirmed 2026-08-05: this blocked
    prediction_set/signal_files/prediction_tabular_files forever, even
    after their cluster-scoped dependency was genuinely uploaded)."""
    return model_key if registry.get(dep_table_name).scope == "cluster_model" else ""


def _compute_round(cache, conn, dataset, cluster, model, cluster_cfg, igvf_cfg, scE2G_dir, data_dir, output_dir, table_name, variant_name):
    """1 + max(round of each dependency), memoized per (dataset, cluster,
    model-or-"", table, variant) -- shared across the whole run() call so a
    dependency reused by many rows (e.g. principal_pseudobulk_set) is only
    ever walked once. Round 1 = no dependencies. Purely for REPORTING/
    output-ordering -- never gates whether a row's payload gets built (see
    module docstring)."""
    key = (dataset, cluster, model or "", table_name, variant_name)
    if key in cache:
        return cache[key]
    table = registry.get(table_name)
    variant = next(v for v in table.variants if v.name == variant_name)
    ctx = Context(dataset, cluster, model, cluster_cfg, igvf_cfg, scE2G_dir, data_dir, output_dir, conn=conn)
    deps = variant.depends_on(ctx)
    if not deps:
        cache[key] = 1
        return 1
    dep_rounds = []
    for dep_table_name, dep_variant in deps:
        # Same scope-aware resolution as _dependency_model, but returning the
        # real model value (or None) for building the dependency's own Context
        # -- not the ledger lookup's "" placeholder.
        dep_model = model if registry.get(dep_table_name).scope == "cluster_model" else None
        dep_rounds.append(
            _compute_round(
                cache, conn, dataset, cluster, dep_model, cluster_cfg, igvf_cfg, scE2G_dir, data_dir, output_dir,
                dep_table_name, dep_variant,
            )
        )
    cache[key] = 1 + max(dep_rounds)
    return cache[key]


def plan_table(conn, reader, table, cluster_keys, cluster_configs, igvf_cfg, scE2G_dir, data_dir, output_dir, mode, round_cache):
    """Returns (to_post, to_patch, to_record, counts).
    to_post/to_patch: list of dicts with keys row_id, alias, payload,
    variant, round, dataset (to_patch additionally has record_id) -- run()
    groups by (dataset, variant, round) to decide output filenames. These
    are the ephemeral, this-run-only working files handed to iu_register.py
    -- ONLY rows not already uploaded-and-unchanged ever appear here, by
    design (see module docstring's note on why this is left alone).

    to_record: list of dicts with keys object_type, record_id, payload --
    every row found already `status=uploaded` with a matching hash (i.e.
    the "unchanged" branch below). Feeds orchestrator.run()'s separate,
    durable per-(object_type, dataset) accumulator -- never handed to
    iu_register.py.

    coverage: one row per (dataset, cluster, model, table, variant) this call
    considered, with its `outcome` and a `reason`. EVERY branch below records
    one, including the ones that skip -- that's the point. Previously a row
    that didn't make it into to_post/to_patch left no trace beyond an
    aggregate counter, so "the manifest is missing this cluster because an
    upstream rule never produced its file" and "this row correctly doesn't
    apply here" were indistinguishable in the output. See
    MANIFEST_GAP_OUTCOMES for which outcomes are real gaps.
    """
    to_post, to_patch, to_record, coverage = [], [], [], []
    counts = {}

    def bump(key):
        counts[key] = counts.get(key, 0) + 1

    def record_coverage(dataset, cluster, model_key, variant_name, outcome, reason=""):
        coverage.append(
            {
                "dataset": dataset, "cluster": cluster, "model": model_key,
                "table": table.name, "variant": variant_name,
                "outcome": outcome, "reason": reason,
            }
        )

    for dataset, cluster, model, cluster_cfg in _iter_scopes(table, cluster_keys, cluster_configs, igvf_cfg):
        ctx = Context(dataset, cluster, model, cluster_cfg, igvf_cfg, scE2G_dir, data_dir, output_dir, conn=conn)
        model_key = model or ""
        for variant in table.variants:
            # enabled() and required_paths() can both do real I/O (e.g.
            # prediction_tabular_files' score-threshold glob/parse, which raises on a
            # missing/duplicate marker file) -- must not let that crash the whole
            # multi-cluster run. The message is preserved rather than swallowed: a
            # duplicate score_threshold_* marker is a real misconfiguration to fix,
            # not a row to quietly drop.
            try:
                is_enabled = variant.enabled(ctx)
            except Exception as e:
                log(f"ENABLED-CHECK FAILED {table.name}/{variant.name} for {dataset}/{cluster}/{model_key}: {e}")
                bump("enabled-check-failed")
                record_coverage(dataset, cluster, model_key, variant.name, "enabled-check-failed", str(e))
                continue
            if not is_enabled:
                # A deliberate semantic exclusion (e.g. "genes" under scATAC). Expected,
                # distinct from a missing file, and never a manifest gap.
                bump("skipped-family-gated")
                record_coverage(dataset, cluster, model_key, variant.name, "skipped-family-gated")
                continue

            try:
                missing_paths = [p for p in variant.required_paths(ctx) if not os.path.exists(p)]
            except Exception as e:
                log(f"PATH-CHECK FAILED {table.name}/{variant.name} for {dataset}/{cluster}/{model_key}: {e}")
                bump("enabled-check-failed")
                record_coverage(dataset, cluster, model_key, variant.name, "enabled-check-failed", str(e))
                continue
            if missing_paths:
                # The row this manifest is supposed to describe has no file behind it.
                # For a manifest-eligible cluster this always means an upstream rule
                # didn't produce its output -- name the path so it's actionable.
                log(f"MISSING FILE {table.name}/{variant.name} for {dataset}/{cluster}/{model_key}: {missing_paths[0]}")
                bump("skipped-missing-file")
                record_coverage(
                    dataset, cluster, model_key, variant.name, "skipped-missing-file", ";".join(missing_paths)
                )
                continue

            item_alias = table.build_alias(ctx, variant.name)

            deferred_on = None
            for dep_table_name, dep_variant in variant.depends_on(ctx):
                dep_model_key = _dependency_model(dep_table_name, model_key)
                dep = state.get_upload(conn, dataset, cluster, dep_model_key, dep_table_name, dep_variant)
                if not dep or dep["status"] != "uploaded":
                    deferred_on = f"{dep_table_name}/{dep_variant or '(default)'}"
                    break

            round_num = _compute_round(
                round_cache, conn, dataset, cluster, model, cluster_cfg, igvf_cfg, scE2G_dir, data_dir, output_dir,
                table.name, variant.name,
            )

            if deferred_on:
                if mode == "upload":
                    # Real submission: never send a payload that cross-references
                    # something not actually live yet -- wait for a later run.
                    log(f"DEFERRED {item_alias}: waiting on {deferred_on}")
                    bump("deferred")
                    record_coverage(dataset, cluster, model_key, variant.name, "deferred", deferred_on)
                    continue
                # preview/validate: build+write it anyway, tagged with its round,
                # so the whole chain can be reviewed in one pass -- see module
                # docstring. Just note it's not actually postable yet.
                log(f"PENDING (round {round_num}) {item_alias}: not yet live -- waiting on {deferred_on}")

            try:
                payload = build_payload(table, variant, ctx, item_alias)
                validate_row(table, variant.name, payload)
            except Exception as e:
                log(f"VALIDATION FAILED {item_alias}: {e}")
                bump("invalid")
                record_coverage(dataset, cluster, model_key, variant.name, "invalid", str(e))
                continue

            new_hash = state.payload_hash(payload)
            existing = state.get_upload(conn, dataset, cluster, model_key, table.name, variant.name)
            if existing and existing["status"] == "uploaded" and existing["payload_hash"] == new_hash:
                bump("unchanged")
                record_coverage(dataset, cluster, model_key, variant.name, "unchanged", item_alias)
                to_record.append(
                    {"dataset": dataset, "object_type": table.object_type, "record_id": existing["portal_id"], "payload": payload}
                )
                continue

            row = state.claim_pending(
                conn, dataset, cluster, model_key, table.name, variant.name, item_alias, new_hash, _now()
            )

            portal_record = reader.get_by_alias(item_alias)
            if portal_record is None:
                to_post.append(
                    {
                        "row_id": row["id"], "alias": item_alias, "payload": payload, "variant": variant.name,
                        "round": round_num, "dataset": dataset,
                    }
                )
                bump("planned-post")
                record_coverage(dataset, cluster, model_key, variant.name, "planned-post", item_alias)
            else:
                record_id = portal_record.get("uuid") or portal_record.get("accession") or portal_record.get("@id")
                to_patch.append(
                    {
                        "row_id": row["id"], "alias": item_alias, "payload": payload, "record_id": record_id,
                        "variant": variant.name, "round": round_num, "dataset": dataset,
                    }
                )
                bump("planned-patch")
                record_coverage(dataset, cluster, model_key, variant.name, "planned-patch", item_alias)

    return to_post, to_patch, to_record, counts, coverage


def _verify_and_record(conn, reader, rows):
    """rows: iterable of dicts with at least row_id/alias (to_post or
    to_patch entries both qualify). Re-checks the portal directly rather
    than trusting iu_register.py's exit code.

    Returns the entries just confirmed live, each augmented with its
    resolved record_id -- run() feeds these into the per-(object_type,
    dataset) accumulator alongside plan_table's "unchanged" rows, so a
    freshly-succeeded real POST/PATCH is reflected there immediately
    rather than waiting for next run's "unchanged" pass to pick it up."""
    verified = []
    for entry in rows:
        # database=True: this runs seconds after invoke_register's real write, in the
        # same process -- the default Elasticsearch-backed read can lag that write and
        # falsely report "not found," which record_result below would then persist as
        # status="failed" (blocking any same-run dependent table via the depends_on
        # gate above) for a row that's actually live. Confirmed against real production
        # data 2026-08-11 while backfilling Prediction Set submitter_comment.
        record = reader.get_by_alias(entry["alias"], database=True)
        if record:
            portal_id = record.get("uuid") or record.get("accession") or record.get("@id")
            state.record_result(conn, entry["row_id"], "uploaded", portal_id=portal_id, now=_now())
            verified.append({**entry, "record_id": portal_id})
        else:
            state.record_result(
                conn, entry["row_id"], "failed", error="not found on portal after upload attempt", now=_now()
            )
    return verified


def run(
    cluster_keys,
    cluster_configs,
    igvf_cfg,
    scE2G_dir,
    data_dir,
    state_db_path,
    manifest_dir,
    mode="preview",
    table_names=None,
    excluded=None,
    igvf_mode=None,
    iu_register_path=portal_client.IU_REGISTER_DEFAULT_PATH,
    all_cluster_configs=None,
    output_dir=None,
):
    """mode:
      "preview"  (default) -- plan + write TSVs only. No call to
                 iu_register.py; the only network contact is the read-only
                 alias existence check.
      "validate" -- also runs iu_register.py --dry-run per TSV (real
                 schema validation/type-casting, still zero portal writes).
      "upload"   -- runs iu_register.py for real, then re-verifies and
                 records outcomes. Only reachable when the caller explicitly
                 asks for it -- see manage_igvf_metadata.py's --mode flag.

    all_cluster_configs: every (dataset, cluster) the pipeline knows about
    (all datasets in the pipeline config, not just this invocation's
    cluster_keys) -- defaults to cluster_configs for callers (e.g. the test
    harness) that only ever know about one scope anyway. Passed to
    cell_metadata.refresh_if_stale so the wholesale multireport GET's data
    gets grouped/validated/cached for every cluster we know about, not
    discarded for every cluster except whichever narrow subset this one
    invocation happens to be uploading -- see cell_metadata.py.
    """
    if mode not in ("preview", "validate", "upload"):
        raise ValueError(f"mode must be one of preview/validate/upload, got {mode!r}")

    conn = state.connect(state_db_path)
    reader = portal_client.PortalReader(igvf_mode=igvf_mode)
    now = _now()
    for dataset, cluster in excluded or []:
        state.mark_excluded(conn, dataset, cluster, "resolve_exclusions.py", now)

    # One multireport GET per stale (24h TTL) cache, covering every cluster we know
    # about -- not just this run's cluster_keys, so a later invocation for a
    # different subset can still be served from cache within the same TTL window.
    # Populates Principal Pseudobulk Set's cell_type/cell_qualifier before that
    # table's payload gets built below.
    cell_metadata.refresh_if_stale(
        conn, reader, set(all_cluster_configs or cluster_configs), all_cluster_configs or cluster_configs
    )

    tables = [t for t in registry.all_specs() if table_names is None or t.name in table_names]
    round_cache = {}  # shared across every table this run -- see _compute_round
    accumulator_entries = []  # collected across every table -- see per-(object_type, dataset) write-out below
    coverage_rows = []  # every (cluster, table, variant) and its outcome -- see write-out below
    report = {}
    for table in tables:
        to_post, to_patch, to_record, counts, coverage = plan_table(
            conn, reader, table, cluster_keys, cluster_configs, igvf_cfg, scE2G_dir, data_dir, output_dir,
            mode, round_cache,
        )
        log(f"{table.name}: {counts}")
        accumulator_entries.extend(to_record)
        coverage_rows.extend(coverage)

        # Group by (dataset, variant, round), not just table: a table whose variants
        # sit at different dependency layers (e.g. prediction_tabular_files' elements/
        # genes vs full vs thresholded vs bedpe) must split into separately
        # round-numbered files -- one round per file, so sorting a dataset's
        # subfolder by filename gives the actual upload order (see module
        # docstring). Splitting by dataset too keeps concurrent datasets' files
        # from ever mixing -- each dataset gets its own manifest_dir subfolder.
        groups = {}
        for entry in to_post:
            groups.setdefault((entry["dataset"], entry["variant"], entry["round"], "post"), []).append(entry)
        for entry in to_patch:
            groups.setdefault((entry["dataset"], entry["variant"], entry["round"], "patch"), []).append(entry)

        table_report = {"counts": counts, "files": []}
        for (dataset, variant_name, round_num, kind), entries in sorted(
            groups.items(), key=lambda kv: (kv[0][0], kv[0][2], kv[0][1])
        ):
            suffix = f"_{variant_name}" if variant_name else ""
            fname = f"round{round_num}_{table.name}{suffix}_{kind}.tsv"
            path = os.path.join(manifest_dir, dataset, fname)
            if kind == "post":
                written = portal_client.write_tsv(path, [e["payload"] for e in entries])
            else:
                written = portal_client.write_tsv(
                    path, [e["payload"] for e in entries], record_ids=[e["record_id"] for e in entries]
                )
            table_report["files"].append(written)
            log(f"  -> {dataset}/{fname}: {len(entries)} row(s)")
            _retire_superseded_sibling(manifest_dir, dataset, table.name, suffix, round_num, kind)

            if mode != "preview" and written:
                result = portal_client.invoke_register(
                    written,
                    table.object_type,
                    patch=(kind == "patch"),
                    dry_run=(mode == "validate"),
                    igvf_mode=igvf_mode,
                    iu_register_path=iu_register_path,
                )
                table_report.setdefault("iu_register_results", []).append(
                    {"file": written, "returncode": result.returncode, "stderr": result.stderr[-2000:]}
                )
                if result.returncode != 0:
                    log(f"iu_register.py FAILED on {written} (exit {result.returncode}): {result.stderr[-500:]}")
                if mode == "upload":
                    verified = _verify_and_record(conn, reader, entries)
                    accumulator_entries.extend(
                        {"dataset": v["dataset"], "object_type": table.object_type, "record_id": v["record_id"], "payload": v["payload"]}
                        for v in verified
                    )

        report[table.name] = table_report

    # Per-(object_type, dataset) accumulator: every alias ever confirmed live, one
    # row per record_id, upserted (never dropped) -- a durable reference shaped for
    # bulk field patches later, deliberately separate from the ephemeral POST/PATCH
    # working files above. Written once per group here (not per table) so an
    # object_type shared by multiple table modules (e.g. "tabular_file") gets one
    # merge-write covering all of them, not several partial ones.
    accumulator_groups = {}
    for entry in accumulator_entries:
        accumulator_groups.setdefault((entry["dataset"], entry["object_type"]), []).append(entry)
    for (dataset, object_type), entries in accumulator_groups.items():
        path = os.path.join(manifest_dir, dataset, f"{object_type}.tsv")
        portal_client.merge_write_tsv(
            path, [e["payload"] for e in entries], record_ids=[e["record_id"] for e in entries]
        )
        log(f"  -> {dataset}/{object_type}.tsv: {len(entries)} row(s) merged")

    # Per-dataset manifest coverage: one row per (cluster, model, table, variant)
    # this run considered, with its outcome and reason. This is the artifact that
    # makes an incomplete manifest legible without reading every manifest TSV --
    # "which clusters are in, and for the ones that aren't, exactly why". Written
    # in every mode, including preview.
    for dataset, rows in _group_by_dataset(coverage_rows).items():
        path = os.path.join(manifest_dir, dataset, MANIFEST_COVERAGE_NAME)
        write_coverage_tsv(path, rows)
        outcomes = {}
        for row in rows:
            outcomes[row["outcome"]] = outcomes.get(row["outcome"], 0) + 1
        log(f"  -> {dataset}/{MANIFEST_COVERAGE_NAME}: {len(rows)} row(s) {outcomes}")

    if mode == "upload":
        # Re-derive and re-push the shareable Cell Annotation table whenever a real
        # upload pass runs -- cheap and idempotent (Synapse auto-versions a repeated
        # store() to the same File), simpler than tracking whether principal_pseudobulk_set
        # rows specifically changed this run.
        cell_metadata.push_to_synapse(cell_metadata.build_shareable_rows(conn), manifest_dir=manifest_dir)

    conn.close()
    return report
