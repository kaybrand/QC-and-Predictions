#!/usr/bin/env python3
"""Merge a merged cluster's component filtered_cell_subsample_metrics.tsv files.

plot_per_cell_qc.R writes filtered_cell_subsample_metrics.tsv per cluster it is
run on, so a MERGED cluster -- one whose pipeline config gives a comma-separated
`pseudobulk_annotation`, e.g. igvf18's mcf7 = "mcf7_1,mcf7_2" -- never gets one
under its merged name. That is the whole reason igvf18/mcf7 is the only real
cluster falling through resolve_exclusions.py's quality gate to the datatable
branch (resolve_exclusions.py:73-91): the AND at :77 fails on the missing file,
not on the guide name.

The datatable branch is the wrong basis for a depth gate, because
{cluster}_per_cell_qc.tsv is UNFILTERED -- it holds every cell, pre-QC -- whereas
filtered_cell_subsample_metrics.tsv is the post-filter aggregate. The branch
compensates by joining against the guide's barcodes, but the honest fix is to
give the merged cluster the same post-filter metrics file every other cluster
has. This script does that, so mcf7 uses the primary branch with real filtered
numbers and the gate needs no logic change at all.

    python merge_cluster_metrics.py \
        --metrics .../plots/igvf18/mcf7_1/filtered_cell_subsample_metrics.tsv \
        --metrics .../plots/igvf18/mcf7_2/filtered_cell_subsample_metrics.tsv \
        --out     .../plots/igvf18/mcf7/filtered_cell_subsample_metrics.tsv

GROUPING BY SUBSAMPLE IS REQUIRED, NOT COSMETIC. A subsample routinely
contributes to several clusters, so the components genuinely overlap: FOUR of
mcf7_1's seven subsamples (IGVFSM1269ANTZ, IGVFSM3772AARV, IGVFSM7026LUVK,
IGVFSM8404EKXH) also appear in mcf7_2. A plain concatenation would emit two rows
for each. The gate sums across rows, so it would still total correctly today --
but the file would be malformed for every other reader, and one that expects one
row per subsample would silently double-count.

Column semantics, so the merge stays arithmetically honest:
  - n_cells / total_fragments / total_RNA_reads are ADDITIVE -> summed. These
    are the only three resolve_exclusions.py and aggregate_qc_stats.py read.
  - mean_frag_per_cell / mean_RNA_per_cell are RATIOS -> recomputed from the
    summed numerator and denominator, never averaged.
  - mean_frip / mean_tss are per-cell means -> recombined as n_cells-WEIGHTED
    means, which is the only way to get the true pooled mean. Averaging two
    means would silently misweight a 3059-cell group against a 1-cell one.

Writes via a temp file plus os.replace, matching
resolve_exclusions.write_cluster_stats_table:197-215.

Note --out is explicit and has no default: the canonical plots/ tree lives under
config["data_dir"], which this pipeline treats as READ-ONLY (common.smk:42-48),
so nothing is written there implicitly.
"""

import argparse
import csv
import os
import sys

SUBSAMPLE = "subsample"
ADDITIVE = ("n_cells", "total_fragments", "total_RNA_reads")
# ratio column -> (numerator, denominator)
RATIOS = {
    "mean_frag_per_cell": ("total_fragments", "n_cells"),
    "mean_RNA_per_cell": ("total_RNA_reads", "n_cells"),
}
# per-cell means, recombined weighted by n_cells
WEIGHTED = ("mean_frip", "mean_tss")
COLUMNS = [SUBSAMPLE, *ADDITIVE, "mean_frag_per_cell", "mean_RNA_per_cell", *WEIGHTED]


def log(msg):
    print(f"[merge_cluster_metrics] {msg}", file=sys.stderr)


def _fmt(value):
    """15 significant digits, matching what R's data.table::fwrite produced in
    the existing files (e.g. 14030.2670807453) so a merged file is visually
    consistent with a natively-written one."""
    if isinstance(value, int):
        return str(value)
    return f"{value:.15g}"


def read_metrics(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path}: missing column(s) {missing}; found {reader.fieldnames}")
        return list(reader)


def merge(paths):
    """Returns (rows, overlaps) -- rows keyed and sorted by subsample, overlaps
    the subsamples that appeared in more than one component."""
    acc = {}
    seen_in = {}
    for path in paths:
        for row in read_metrics(path):
            sub = row[SUBSAMPLE]
            seen_in.setdefault(sub, []).append(os.path.basename(os.path.dirname(path)))
            cur = acc.setdefault(sub, {c: 0 for c in ADDITIVE} | {f"_w_{c}": 0.0 for c in WEIGHTED})
            n_cells = int(row["n_cells"])
            for col in ADDITIVE:
                cur[col] += int(row[col])
            for col in WEIGHTED:
                # Weight by this component's own n_cells before pooling.
                cur[f"_w_{col}"] += float(row[col]) * n_cells

    rows = []
    for sub in sorted(acc):
        cur = acc[sub]
        out = {SUBSAMPLE: sub}
        for col in ADDITIVE:
            out[col] = cur[col]
        n_cells = cur["n_cells"]
        for col, (num, den) in RATIOS.items():
            out[col] = (cur[num] / cur[den]) if cur[den] else 0
        for col in WEIGHTED:
            out[col] = (cur[f"_w_{col}"] / n_cells) if n_cells else 0
        rows.append(out)
    overlaps = {s: v for s, v in seen_in.items() if len(v) > 1}
    return rows, overlaps


def write_metrics(rows, out_path):
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        for row in rows:
            writer.writerow([_fmt(row[c]) if c != SUBSAMPLE else row[c] for c in COLUMNS])
    os.replace(tmp, out_path)


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "--metrics",
        action="append",
        required=True,
        help="a component cluster's filtered_cell_subsample_metrics.tsv (repeat for each)",
    )
    p.add_argument("--out", required=True, help="output path for the merged metrics file")
    p.add_argument("--force", action="store_true", help="overwrite --out if it exists")
    args = p.parse_args(argv)

    if len(args.metrics) < 2:
        log("WARNING: only one --metrics given; merging a single file just copies it")
    for path in args.metrics:
        if not os.path.exists(path):
            log(f"ERROR: no such file: {path}")
            return 2
    if os.path.exists(args.out) and not args.force:
        log(f"ERROR: {args.out} exists (use --force to overwrite)")
        return 2

    rows, overlaps = merge(args.metrics)
    write_metrics(rows, args.out)

    totals = {c: sum(r[c] for r in rows) for c in ADDITIVE}
    log(f"merged {len(args.metrics)} component(s) -> {len(rows)} subsample row(s): {args.out}")
    log(f"  totals the quality gate will read: {totals}")
    if overlaps:
        log(f"  {len(overlaps)} subsample(s) appeared in more than one component and were SUMMED:")
        for sub, where in sorted(overlaps.items()):
            log(f"    {sub}: {', '.join(where)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
