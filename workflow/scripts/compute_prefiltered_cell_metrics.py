#!/usr/bin/env python3
"""
Derive a filtered_cell_subsample_metrics.tsv-compatible row from an ATAC
fragments file that has already been filtered/QC'd upstream (e.g. catlas's
per-cluster fragments pulled straight from the IGVF Data Portal) -- no QC
guide, no per-cell QC table exists for these clusters.

Written in the exact column schema QC_pseudobulks/scripts/plotting_scripts/
plot_per_cell_qc.R produces (subsample, n_cells, total_fragments,
total_RNA_reads, mean_frag_per_cell, mean_RNA_per_cell, mean_frip, mean_tss),
one row representing the whole cluster (there's no real subsampling here),
so resolve_exclusions.py/aggregate_qc_stats.py consume it unmodified.

Usage:
    python compute_prefiltered_cell_metrics.py \
        --frag-file /path/to/fragments.tsv.gz \
        --subsample-name my_cluster \
        --out /path/to/filtered_cell_subsample_metrics.tsv
"""

import argparse
import csv
import gzip
import os

HEADER = [
    "subsample", "n_cells", "total_fragments", "total_RNA_reads",
    "mean_frag_per_cell", "mean_RNA_per_cell", "mean_frip", "mean_tss",
]


def parse_args():
    p = argparse.ArgumentParser(description="Derive cell/fragment counts from an already-filtered ATAC fragments file.")
    p.add_argument("--frag-file", required=True, help="Gzipped fragments TSV (chrom, start, end, barcode, count).")
    p.add_argument("--subsample-name", required=True, help="Value for the 'subsample' column (this cluster's name).")
    p.add_argument("--out", required=True, help="Output path for the metrics TSV.")
    return p.parse_args()


def count_cells_and_fragments(frag_file):
    """Streams the fragments file rather than `sort -u`-ing it: memory stays
    bounded to the (small) number of distinct barcodes instead of scaling
    with the (potentially multi-million) line count."""
    barcodes = set()
    total_fragments = 0
    with gzip.open(frag_file, "rt") as f:
        for line in f:
            total_fragments += 1
            barcodes.add(line.split("\t", 4)[3])
    return len(barcodes), total_fragments


def main():
    args = parse_args()
    n_cells, total_fragments = count_cells_and_fragments(args.frag_file)

    row = {
        "subsample": args.subsample_name,
        "n_cells": n_cells,
        "total_fragments": total_fragments,
        "total_RNA_reads": 0,  # no RNA for these clusters
        "mean_frag_per_cell": total_fragments / n_cells if n_cells else "NA",
        "mean_RNA_per_cell": 0,
        "mean_frip": "NA",  # not computable without a per-cell QC table
        "mean_tss": "NA",
    }

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    tmp_path = f"{args.out}.tmp.{os.getpid()}"
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)
    os.replace(tmp_path, args.out)  # atomic: concurrent readers never see a partial file


if __name__ == "__main__":
    main()
