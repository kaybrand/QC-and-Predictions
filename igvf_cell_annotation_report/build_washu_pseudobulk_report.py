#!/usr/bin/env python3
"""
Flattens get_washu_pseudobulk_report.py's raw multireport JSON into a
human-readable TSV -- one row per Yang Li (WashU) PseudobulkSet.

Replaces an earlier ad hoc, unreproducible version of this same flattening
(never checked in as a script) that dropped SampleTermID entirely -- that's
what left the state.db seeding step (seed_catlas_cell_annotations.py, which
now does its own live fetch instead of reading this TSV) with no real
ontology term ID to cache, only the placeholder "TODO: ontology_id..."
string. This script exists mainly for download_and_verify_fragments.py
(reads FragmentsFile_Href/FragmentsFile_MD5/FragmentsFile_Aliases) and for
human inspection -- keep those three column names stable if this ever needs
extending again.

Usage:
    python build_washu_pseudobulk_report.py --out washu_pseudobulk_report.tsv
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from get_washu_pseudobulk_report import fetch_multireport  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts"))
from igvf_metadata.cell_metadata import (  # noqa: E402
    _cl_id_from_cell_type,
    _term_id_from_cell_type,
    _term_name_from_cell_type,
)

HEADER = [
    "PseudobulkSet_ID", "Alias", "CellAnnotation", "CellQualifier",
    "SampleTermName", "SampleTermID", "SampleTermCLID", "Status",
    "FragmentsFile_Aliases", "FragmentsFile_Status", "FragmentsFile_Href",
    "FragmentsFile_MD5", "N_Fragments_Files",
]


def fragments_files(row):
    return [f for f in row.get("files", []) if f.get("content_type") == "fragments"]


def build_row(row):
    cell_type = row.get("cell_type")
    frags = fragments_files(row)
    frag = frags[0] if frags else {}
    return {
        "PseudobulkSet_ID": row["@id"],
        "Alias": row["aliases"][0] if row.get("aliases") else "",
        "CellAnnotation": row.get("cell_annotation", ""),
        "CellQualifier": row.get("cell_qualifier") or "",
        "SampleTermName": _term_name_from_cell_type(cell_type) or "",
        "SampleTermID": _term_id_from_cell_type(cell_type) or "",
        "SampleTermCLID": _cl_id_from_cell_type(cell_type) or "",
        "Status": row.get("status", ""),
        "FragmentsFile_Aliases": ",".join(frag.get("aliases", [])),
        "FragmentsFile_Status": frag.get("status", ""),
        "FragmentsFile_Href": frag.get("href", ""),
        "FragmentsFile_MD5": frag.get("md5sum", ""),
        "N_Fragments_Files": len(frags),
    }


def main():
    p = argparse.ArgumentParser(description="Flatten the WashU PseudobulkSet multireport into a TSV.")
    p.add_argument("--out", required=True, help="Output TSV path.")
    args = p.parse_args()

    rows = fetch_multireport()
    print(f"# {len(rows)} PseudobulkSet(s) returned", file=sys.stderr)

    with open(args.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(build_row(row))


if __name__ == "__main__":
    main()
