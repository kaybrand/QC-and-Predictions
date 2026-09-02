#!/usr/bin/env python3
"""
Seed a brand-new, CATlas-local igvf_metadata state.db cell_annotations cache
via a LIVE IGVF Data Portal fetch (lab.@id=/labs/yang-li/), instead of the
live portal-refresh path (igvf_metadata.cell_metadata.refresh_if_stale) that
reformat.smk's portal_cell_metadata() normally expects: that refresh matches
portal rows to local (dataset, cluster, subsample) triples by the alias
suffix "{dataset}-{cluster}-{subsample}", which catlas's WashU cluster IDs
(e.g. EMSN_234, Bgl_51) don't follow -- there's no local subsample/QC-guide
concept for these clusters at all (see workflow/rules/prefiltered_fragments.smk).

An earlier version of this script read washu_pseudobulk_report.tsv instead
of fetching live -- that TSV never had a real ontology term ID column (lost
somewhere in an ad hoc, unreproducible flattening step), so cl_id/term_id
got written as a "TODO: ontology_id..." placeholder for every cluster. This
version fetches the portal directly (reusing get_washu_pseudobulk_report.py's
own query) and reads cl_id/term_id/term_name off the real embedded
cell_type object, via the exact same extraction helpers
igvf_metadata.cell_metadata uses for every other dataset's Cell Annotation
cache -- so catlas gets real ontology data (e.g. "PCL:0051300"), not a stub.

This writes to config["igvf"]["state_db_path"] as set in
catlas_pipeline_config.yaml -- a NEW database file local to this worktree,
never the shared one at
QC_pseudobulks-igvf-portal-submission/resources/igvf_metadata_state.db (the
three worktrees are separate forks with separate submission targets).

Idempotent: upsert_cell_annotation is ON CONFLICT(dataset, cluster) DO UPDATE,
so re-running this whenever the portal data changes is safe.

Usage:
    python seed_catlas_cell_annotations.py \
        --state-db-path /path/to/CATlas_predictions/resources/igvf_metadata_state.db
"""

import argparse
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

DATASET = "catlas"


def parse_args():
    p = argparse.ArgumentParser(description="Seed a CATlas-local state.db cell_annotations cache via a live IGVF Portal fetch.")
    p.add_argument("--state-db-path", required=True, help="Path to this worktree's own igvf_metadata_state.db (created if missing).")
    return p.parse_args()


def cluster_from_alias(alias):
    # e.g. "yang-li:catlas-human-pseudobulk-EMSN_234" -> "EMSN_234"
    return alias.rsplit("-", 1)[-1]


def main():
    args = parse_args()
    conn = igvf_state.connect(args.state_db_path)
    now = datetime.now(timezone.utc).isoformat()

    rows = fetch_multireport()
    print(f"[seed_catlas_cell_annotations] fetched {len(rows)} PseudobulkSet(s)")

    for row in rows:
        cluster = cluster_from_alias(row["aliases"][0])
        cell_type = row.get("cell_type")
        igvf_state.upsert_cell_annotation(
            conn, DATASET, cluster,
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
        print(f"[seed_catlas_cell_annotations] seeded {DATASET}/{cluster}: term_id={_term_id_from_cell_type(cell_type)}")

    conn.close()


if __name__ == "__main__":
    main()
