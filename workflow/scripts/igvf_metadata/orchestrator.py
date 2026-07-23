"""Core upload loop, shared by the pipeline-integrated hook and the
standalone scanner -- both just assemble a different list of (dataset,
cluster) scope-keys and hand it to run(). See registry.py and state.py for
the data model, and portal_client.py for why real submission goes through
iu_register.py rather than direct API calls.

Two phases per table, every run:
  1. plan_table(): for every enabled row not already uploaded-and-unchanged,
     check depends_on (against the local ledger), build+validate its
     payload, hash it, and classify it as "post" (no existing portal
     record for this alias) or "patch" (exists, and the hash changed) via
     a live, read-only alias lookup. This is the ONLY network contact in
     "preview" mode.
  2. Write one post.tsv/patch.tsv per table under manifest_dir -- always,
     in every mode, so the assembled tables are there "for the user to
     peruse" regardless of what happens next. Then, depending on `mode`:
       - "preview" (default): stop here. No call to iu_register.py.
       - "validate": call iu_register.py --dry-run against each TSV --
         real igvf_utils schema validation/type-casting, still zero writes.
       - "upload": call iu_register.py for real, then re-verify every
         intended row via a fresh alias lookup (iu_register.py doesn't
         isolate row failures within one file, so its exit code alone
         isn't trustworthy) and update the ledger accordingly.

Row ordering across tables/variants falls out of just re-running: a row
deferred on a dependency is picked up automatically next time run() is
called, whether that's the pipeline hook firing for a later cluster or the
scanner's periodic reconciliation.

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

import os
import sys
from datetime import datetime, timezone

from . import cell_metadata, portal_client, registry, state
from .context import Context


def log(msg):
    print(f"[igvf_metadata] {msg}", file=sys.stderr)


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


def plan_table(conn, reader, table, cluster_keys, cluster_configs, igvf_cfg, scE2G_dir):
    """Returns (to_post, to_patch, counts).
    to_post:  list of (state_row_id, item_alias, payload)
    to_patch: list of (state_row_id, item_alias, payload, record_id)
    """
    to_post, to_patch = [], []
    counts = {}

    def bump(key):
        counts[key] = counts.get(key, 0) + 1

    for dataset, cluster, model, cluster_cfg in _iter_scopes(table, cluster_keys, cluster_configs, igvf_cfg):
        ctx = Context(dataset, cluster, model, cluster_cfg, igvf_cfg, scE2G_dir, conn=conn)
        model_key = model or ""
        for variant in table.variants:
            # enabled() can do real I/O (e.g. prediction_tabular_files' score-threshold
            # glob/parse) that can raise on a genuinely bad input (missing/duplicate
            # marker file) -- must not let that crash the whole multi-cluster run.
            try:
                is_enabled = variant.enabled(ctx)
            except Exception as e:
                log(f"ENABLED-CHECK FAILED {table.name}/{variant.name} for {dataset}/{cluster}/{model_key}: {e}")
                bump("enabled-check-failed")
                continue
            if not is_enabled:
                bump("skipped-disabled")
                continue

            item_alias = table.build_alias(ctx, variant.name)

            deferred_on = None
            for dep_table_name, dep_variant in variant.depends_on(ctx):
                dep = state.get_upload(conn, dataset, cluster, model_key, dep_table_name, dep_variant)
                if not dep or dep["status"] != "uploaded":
                    deferred_on = f"{dep_table_name}/{dep_variant or '(default)'}"
                    break
            if deferred_on:
                log(f"DEFERRED {item_alias}: waiting on {deferred_on}")
                bump("deferred")
                continue

            try:
                payload = build_payload(table, variant, ctx, item_alias)
                validate_row(table, variant.name, payload)
            except Exception as e:
                log(f"VALIDATION FAILED {item_alias}: {e}")
                bump("invalid")
                continue

            new_hash = state.payload_hash(payload)
            existing = state.get_upload(conn, dataset, cluster, model_key, table.name, variant.name)
            if existing and existing["status"] == "uploaded" and existing["payload_hash"] == new_hash:
                bump("unchanged")
                continue

            row = state.claim_pending(
                conn, dataset, cluster, model_key, table.name, variant.name, item_alias, new_hash, _now()
            )

            portal_record = reader.get_by_alias(item_alias)
            if portal_record is None:
                to_post.append((row["id"], item_alias, payload))
                bump("planned-post")
            else:
                record_id = portal_record.get("uuid") or portal_record.get("accession") or portal_record.get("@id")
                to_patch.append((row["id"], item_alias, payload, record_id))
                bump("planned-patch")

    return to_post, to_patch, counts


def _verify_and_record(conn, reader, rows):
    """rows: iterable of tuples whose first two elements are (state_row_id,
    item_alias) -- works for both to_post and to_patch entries. Re-checks
    the portal directly rather than trusting iu_register.py's exit code."""
    for entry in rows:
        row_id, item_alias = entry[0], entry[1]
        record = reader.get_by_alias(item_alias)
        if record:
            portal_id = record.get("uuid") or record.get("accession") or record.get("@id")
            state.record_result(conn, row_id, "uploaded", portal_id=portal_id, now=_now())
        else:
            state.record_result(conn, row_id, "failed", error="not found on portal after upload attempt", now=_now())


def run(
    cluster_keys,
    cluster_configs,
    igvf_cfg,
    scE2G_dir,
    state_db_path,
    manifest_dir,
    mode="preview",
    table_names=None,
    excluded=None,
    igvf_mode=None,
    iu_register_path=portal_client.IU_REGISTER_DEFAULT_PATH,
    all_cluster_configs=None,
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
    report = {}
    for table in tables:
        to_post, to_patch, counts = plan_table(conn, reader, table, cluster_keys, cluster_configs, igvf_cfg, scE2G_dir)

        post_path = os.path.join(manifest_dir, f"{table.name}_post.tsv")
        patch_path = os.path.join(manifest_dir, f"{table.name}_patch.tsv")
        written_post = portal_client.write_tsv(post_path, [p for _, _, p in to_post]) if to_post else None
        written_patch = (
            portal_client.write_tsv(
                patch_path, [p for _, _, p, _ in to_patch], record_ids=[r for _, _, _, r in to_patch]
            )
            if to_patch
            else None
        )

        table_report = {"counts": counts, "post_tsv": written_post, "patch_tsv": written_patch}
        log(f"{table.name}: {counts} (post_tsv={written_post}, patch_tsv={written_patch})")

        if mode != "preview":
            iu_results = []
            for path, patch_flag, rows in ((written_post, False, to_post), (written_patch, True, to_patch)):
                if not path:
                    continue
                result = portal_client.invoke_register(
                    path,
                    table.object_type,
                    patch=patch_flag,
                    dry_run=(mode == "validate"),
                    igvf_mode=igvf_mode,
                    iu_register_path=iu_register_path,
                )
                iu_results.append({"file": path, "returncode": result.returncode, "stderr": result.stderr[-2000:]})
                if result.returncode != 0:
                    log(f"iu_register.py FAILED on {path} (exit {result.returncode}): {result.stderr[-500:]}")
                if mode == "upload":
                    _verify_and_record(conn, reader, rows)
            table_report["iu_register_results"] = iu_results

        report[table.name] = table_report

    if mode == "upload":
        # Re-derive and re-push the shareable Cell Annotation table whenever a real
        # upload pass runs -- cheap and idempotent (Synapse auto-versions a repeated
        # store() to the same File), simpler than tracking whether principal_pseudobulk_set
        # rows specifically changed this run.
        cell_metadata.push_to_synapse(cell_metadata.build_shareable_rows(conn), manifest_dir=manifest_dir)

    conn.close()
    return report
