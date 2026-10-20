#!/usr/bin/env python3
"""Set, list, or clear hand-curated Cell Annotation overrides in state.db.

An override replaces the cell_annotation/cell_qualifier that
cell_metadata.derive_scopes would otherwise derive for one (dataset, cluster).
It lives in its own table (cell_annotation_overrides), which no Portal refresh
writes, so it survives every fetch and is re-applied on every derive. The
driver's warm stage then carries it into the CellAnnotation snapshot
(reformatted file headers) and the IGVF manifests (cell_qualifier). Cell type
terms (cl_id/term_id/term_name) always stay the Portal's.

    python set_cell_annotation_override.py --state-db resources/igvf_metadata_state.db \\
        --dataset IGVFDS3445MTOM --cluster hcasm_calcification \\
        --cell-annotation "calcification HCASMC ..." --cell-qualifier "calcification HCASMC ..." \\
        --reason "merged d6+d10: no day in the annotation"
    python set_cell_annotation_override.py --state-db ... --list
    python set_cell_annotation_override.py --state-db ... --dataset D --cluster C --clear

Takes the same flock as the driver ({state_db}.warmlock): state.db is WAL-mode
SQLite on Lustre, so it must have one writer at a time across hosts. Use an env
interpreter (SQLite >= 3.24), not the login node's system python.
"""

import argparse
import fcntl
import os
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from igvf_metadata import state  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--state-db", required=True)
    p.add_argument("--dataset")
    p.add_argument("--cluster")
    p.add_argument("--cell-annotation")
    p.add_argument("--cell-qualifier", help="omit for no qualifier")
    p.add_argument("--reason", help="why this override exists (required when setting)")
    p.add_argument("--list", action="store_true", help="print every override and exit")
    p.add_argument("--clear", action="store_true", help="remove the override for --dataset/--cluster")
    args = p.parse_args(argv)

    if not args.list and not (args.dataset and args.cluster):
        p.error("--dataset and --cluster are required unless --list")
    if not (args.list or args.clear) and not (args.cell_annotation and args.reason):
        p.error("setting an override needs --cell-annotation and --reason")

    with open(f"{args.state_db}.warmlock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        conn = state.connect(args.state_db)
        if args.clear:
            state.delete_cell_annotation_override(conn, args.dataset, args.cluster)
        elif not args.list:
            state.upsert_cell_annotation_override(
                conn, args.dataset, args.cluster, args.cell_annotation, args.cell_qualifier,
                args.reason, datetime.now(timezone.utc).isoformat(),
            )
        for row in sorted(state.all_cell_annotation_overrides(conn), key=lambda r: (r["dataset"], r["cluster"])):
            print(f"{row['dataset']}/{row['cluster']}\n  cell_annotation: {row['cell_annotation']}\n"
                  f"  cell_qualifier:  {row['cell_qualifier']}\n  reason: {row['reason']} (set {row['set_at']})")
        conn.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
