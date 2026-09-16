#!/usr/bin/env python3
"""
Seed this CATlas-local igvf_metadata state.db's cell_annotations cache via a
LIVE IGVF Data Portal fetch (lab.@id=/labs/yang-li/), instead of the live
portal-refresh path (igvf_metadata.cell_metadata.refresh_if_stale) that
reformat.smk's portal_cell_metadata() normally expects: that refresh matches
portal rows to local (dataset, cluster, subsample) triples by the alias
suffix "{dataset}-{cluster}-{subsample}", which requires each cluster to
carry a QC guide (for subsample derivation). CATlas clusters have no QC
guide and no subsample concept at all -- each cluster is already a single,
uniquely-aliased, curated PseudobulkSet on the Portal (no merging of raw
primary pseudobulks happens for CATlas), so the QC-guide/subsample machinery
is unnecessary and, for these clusters, unreadable
(workflow/rules/prefiltered_fragments.smk; confirmed 2026-09-16: every
CATlas cluster's derive_scopes call fails with `unreadable_qc_guide:
KeyError: 'qc_guide'`, 0/258 resolved).

Fetches every non-deleted yang-li PseudobulkSet regardless of species (627
as of 2026-09-16: 258 human + 369 mouse, confirmed zero un-parseable
aliases) and reads cl_id/term_id/term_name off the real embedded cell_type
object, via the exact same extraction helpers igvf_metadata.cell_metadata
uses for every other dataset's Cell Annotation cache -- so CATlas gets real
ontology data (e.g. "PCL:0051300"), not a stub.

Dataset key is `f"catlas-{species}"` (e.g. "catlas-human", "catlas-mouse"),
matching every real pipeline config's own dataset key -- NOT the flat
"catlas" this script used previously, which orphaned 10 rows under a key
nothing else in this repo ever reads (catlas-human_pipeline_config.yaml,
results/igvf_metadata/catlas-human_*, portal_files.dataset all use
"catlas-human"). Species and cluster are both parsed from the alias by
splitting on the fixed "-pseudobulk-" marker, not by rsplit("-", 1) as
before -- several human cluster names contain internal hyphens of their own
(e.g. "supercluster-Deep-layer_corticothalamic_and_6b",
"supercluster-Upper-layer_intratelencephalic"), which rsplit("-", 1) would
mis-parse (e.g. down to just "6b").

This writes to config["igvf"]["state_db_path"] as set in
catlas-human_pipeline_config.yaml -- a database file local to this worktree,
never the shared one at
QC_pseudobulks-igvf-portal-submission/resources/igvf_metadata_state.db (the
three worktrees are separate forks with separate submission targets).

Idempotent: upsert_cell_annotation is ON CONFLICT(dataset, cluster) DO
UPDATE, so re-running this whenever the portal data changes is safe. Writes
are wrapped in the same state.db.warmlock flock every other stage that
touches state.db uses (run_pipeline.py's state_db_lock), so a concurrent
driver invocation can never interleave writes with this one.

Usage:
    python seed_catlas_cell_annotations.py \
        --state-db-path /path/to/CATlas_predictions/resources/igvf_metadata_state.db
"""

import argparse
import fcntl
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(__file__))
from get_washu_pseudobulk_report import fetch_multireport  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))
from igvf_metadata import state as igvf_state  # noqa: E402
from igvf_metadata.cell_metadata import (  # noqa: E402
    _cl_id_from_cell_type,
    _term_id_from_cell_type,
    _term_name_from_cell_type,
)

ALIAS_MARKER = "-pseudobulk-"


def parse_args():
    p = argparse.ArgumentParser(description="Seed this CATlas worktree's state.db cell_annotations cache via a live IGVF Portal fetch.")
    p.add_argument("--state-db-path", required=True, help="Path to this worktree's own igvf_metadata_state.db (created if missing).")
    return p.parse_args()


def dataset_and_cluster_from_alias(alias):
    """"yang-li:catlas-human-pseudobulk-supercluster-Deep-layer_near-projecting"
    -> ("catlas-human", "supercluster-Deep-layer_near-projecting").

    Splits on the fixed ALIAS_MARKER, not the last hyphen: several cluster
    names contain their own internal hyphens (see module docstring), which a
    rsplit("-", 1) approach would truncate incorrectly.

    Raises ValueError (loud, not skipped) if the marker is missing -- a
    future alias convention change must be seen and handled deliberately,
    not silently mis-parsed.
    """
    suffix = alias.split(":", 1)[-1]
    if ALIAS_MARKER not in suffix:
        raise ValueError(f"alias {alias!r} has no {ALIAS_MARKER!r} marker -- convention may have changed")
    prefix, cluster = suffix.split(ALIAS_MARKER, 1)
    if not prefix.startswith("catlas-"):
        raise ValueError(f"alias {alias!r}'s prefix {prefix!r} doesn't start with 'catlas-'")
    species = prefix[len("catlas-"):]
    return f"catlas-{species}", cluster


def main():
    args = parse_args()
    now = datetime.now(timezone.utc).isoformat()

    rows = fetch_multireport()
    print(f"[seed_catlas_cell_annotations] fetched {len(rows)} PseudobulkSet(s)")

    parsed = []
    errors = []
    for row in rows:
        alias = row["aliases"][0]
        try:
            dataset, cluster = dataset_and_cluster_from_alias(alias)
        except ValueError as e:
            errors.append((alias, str(e)))
            continue
        parsed.append((dataset, cluster, row))

    if errors:
        print(f"[seed_catlas_cell_annotations] {len(errors)} alias(es) could not be parsed -- NOT seeded:")
        for alias, why in errors:
            print(f"  {alias}: {why}")

    by_dataset = {}
    for dataset, cluster, row in parsed:
        by_dataset.setdefault(dataset, []).append((cluster, row))

    lock_path = f"{args.state_db_path}.warmlock"
    fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o644)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX)
        conn = igvf_state.connect(args.state_db_path)
        try:
            for dataset, entries in sorted(by_dataset.items()):
                for cluster, row in sorted(entries):
                    cell_type = row.get("cell_type")
                    igvf_state.upsert_cell_annotation(
                        conn, dataset, cluster,
                        cell_annotation=row.get("cell_annotation", ""),
                        cl_id=_cl_id_from_cell_type(cell_type) or "TODO: ontology_id (not returned by portal)",
                        term_id=_term_id_from_cell_type(cell_type) or "TODO: ontology_id (not returned by portal)",
                        term_name=_term_name_from_cell_type(cell_type) or row.get("cell_annotation", ""),
                        cell_qualifier=row.get("cell_qualifier") or None,
                        portal_samples=row["@id"],
                        all_primary_released=(row.get("status") == "released"),
                        principal_uploaded=False,
                        principal_alias=None,
                        now=now,
                    )
                print(f"[seed_catlas_cell_annotations] {dataset}: seeded {len(entries)} cluster(s)")
        finally:
            conn.close()
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)


if __name__ == "__main__":
    main()
