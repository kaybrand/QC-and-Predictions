#!/usr/bin/env python3
"""
Filter and concatenate ATAC fragment files by QC-passing barcodes for a given cell type.

Usage:
    python filter_atac_fragments.py \
        --qc-guide /path/to/qc_filter_guide.tsv.gz \
        --pseudobulks /path/to/pseudobulks/ \
        --cell-type MycellType \
        --chrom-sizes /path/to/hg38.chrom.sizes \
        --out /path/to/output/fragments.tsv.gz
"""

import argparse
import gzip
import os
import glob
import subprocess
import sys
import tempfile


def parse_args():
    p = argparse.ArgumentParser(description="Filter ATAC fragments by QC barcodes and cell type.")
    p.add_argument("--qc-guide",     required=True, help="Path to gzipped QC filter guide TSV.")
    p.add_argument("--pseudobulks",  required=True, help="Path to pseudobulks folder.")
    p.add_argument("--cell-type",    required=True, help="Exact cell type identifier used in pseudobulk directory names (the string between 'annotation-' and '-IGVF').")
    p.add_argument("--chrom-sizes",  required=True, help="Chromosome sizes file defining the expected sort order (e.g. hg38.chrom.sizes).")
    p.add_argument("--out",          required=True, help="Output path for the filtered, sorted, bgzipped fragment file.")
    return p.parse_args()


def load_passing_barcodes(qc_guide_path: str) -> set:
    """
    Read the QC filter guide and return a set of full barcode strings
    (e.g. CCATATTTCGATAACC_IGVFSM4662QKFQ). Matching is always on the full
    name -- no trimming, no remapping (mirrors filter_rna_counts.py).
    The 16bp nucleotide prefix alone is NOT unique: the same 10x barcode
    sequence can recur across different multiplexed lanes/analyses within
    a single subsample, so truncating to 16 chars collapses distinct cells
    onto the same key and silently corrupts output (confirmed in production:
    one cell's fragments got misattributed to another, and the first was
    dropped entirely).
    """
    barcodes = set()
    opener = gzip.open if qc_guide_path.endswith(".gz") else open
    with opener(qc_guide_path, "rt") as fh:
        fh.readline()  # skip header
        for line in fh:
            line = line.strip()
            if not line:
                continue
            barcodes.add(line.split("\t")[0])
    print(f"[info] Loaded {len(barcodes)} passing barcodes from QC guide.", file=sys.stderr)
    return barcodes


def find_fragment_dirs(pseudobulks_dir: str, cell_type: str) -> list:
    """Return all direct child directories matching annotation-{cellType}-IGVF*.
    Falls back to annotation-{cellType} (no subsample suffix) if no -IGVF* dirs found."""
    pattern = os.path.join(pseudobulks_dir, f"annotation-{cell_type}-IGVF*")
    dirs = [d for d in glob.glob(pattern) if os.path.isdir(d)]
    if dirs:
        print(f"[info] Found {len(dirs)} matching directories.", file=sys.stderr)
        return sorted(dirs)
    # Fall back to annotation-only pseudobulks (no subsample suffix)
    fallback = os.path.join(pseudobulks_dir, f"annotation-{cell_type}")
    if os.path.isdir(fallback):
        print(f"[info] No -IGVF* directories found; using annotation-level pseudobulk: {fallback}", file=sys.stderr)
        return [fallback]
    print(f"[warning] No directories found matching: {pattern}", file=sys.stderr)
    return []


def filter_fragments(frag_path: str, passing_barcodes: set, out_fh, found_barcodes: set):
    """
    Stream through a gzipped fragment file, writing rows whose 4th-column
    barcode (full string, e.g. CCATATTTCGATAACC_IGVFSM4662QKFQ) is in
    passing_barcodes to out_fh, unchanged. Updates found_barcodes in place.
    Returns (n_pass, n_total).
    """
    n_pass = 0
    n_total = 0
    opener = gzip.open if frag_path.endswith(".gz") else open
    with opener(frag_path, "rt") as fh:
        for line_num, line in enumerate(fh, start=1):
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 5:
                raise ValueError(
                    f"Malformed fragment record in {frag_path}, line {line_num}: "
                    f"expected 5 tab-separated fields (chrom, start, end, barcode, "
                    f"duplicate count), got {len(fields)}: {line.rstrip(chr(10))!r}"
                )
            n_total += 1
            bc = fields[3]
            if bc in passing_barcodes:
                out_fh.write("\t".join(fields) + "\n")
                found_barcodes.add(bc)
                n_pass += 1
    return n_pass, n_total


def check_tool(tool: str):
    if subprocess.run(["which", tool], capture_output=True).returncode != 0:
        print(f"[error] Required tool not found in PATH: {tool}", file=sys.stderr)
        sys.exit(1)


def main():
    args = parse_args()

    for tool in ("sort-bed", "bgzip", "tabix"):
        check_tool(tool)

    passing_barcodes = load_passing_barcodes(args.qc_guide)
    frag_dirs        = find_fragment_dirs(args.pseudobulks, args.cell_type)

    if not frag_dirs:
        print("[error] No fragment directories found. Exiting.", file=sys.stderr)
        sys.exit(1)

    # Ensure output directory exists
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)

    # Use a temp file to collect raw (unsorted) passing rows
    found_barcodes = set()
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tmp:
        tmp_path = tmp.name
        total_pass = 0
        total_rows = 0
        for d in frag_dirs:
            frag_path = os.path.join(d, "fragments.tsv.gz")
            if not os.path.isfile(frag_path):
                print(f"[warning] fragments.tsv.gz not found in {d}, skipping.", file=sys.stderr)
                continue
            n_pass, n_total = filter_fragments(frag_path, passing_barcodes, tmp, found_barcodes)
            pct = 100 * n_pass / n_total if n_total > 0 else 0.0
            print(f"[info]   {os.path.basename(d)}: {n_pass}/{n_total} rows retained ({pct:.1f}%)", file=sys.stderr)
            total_pass += n_pass
            total_rows += n_total

    total_pct = 100 * total_pass / total_rows if total_rows > 0 else 0.0
    print(f"[info] Total passing rows: {total_pass}/{total_rows} ({total_pct:.1f}%)", file=sys.stderr)

    # Sanity check: report any QC-passing barcodes never seen in any fragment file
    missing_barcodes = passing_barcodes - found_barcodes
    if missing_barcodes:
        print(f"[error] {len(missing_barcodes)} barcode(s) from the QC guide were not found in any fragment file:", file=sys.stderr)
        for bc in sorted(missing_barcodes):
            print(f"  {bc}", file=sys.stderr)
        sys.exit(1)
    else:
        print(f"[info] Sanity check passed: all {len(passing_barcodes)} QC-guide barcodes were found in at least one fragment file.", file=sys.stderr)

    # Ensure output ends with .gz for bgzip
    out_path = args.out if args.out.endswith(".gz") else args.out + ".gz"

    # Sort respecting chrom order from sizes file -> bgzip -> tabix index
    print("[info] Sorting by chromosome order, bgzipping, and indexing ...", file=sys.stderr)

    split_dir = tmp_path + ".split"
    os.makedirs(split_dir, exist_ok=True)

    # Split the concatenated temp file by chromosome (col 1)
    with open(tmp_path, "r") as fh:
        for line in fh:
            chrom = line.split("\t")[0]
            with open(os.path.join(split_dir, f"{chrom}.bed"), "a") as chrom_fh:
                chrom_fh.write(line)

    # Load chrom sizes order; warn once per chromosome missing from sizes file
    sizes_chroms = []
    with open(args.chrom_sizes, "r") as fh:
        for line in fh:
            line = line.strip()
            if line:
                sizes_chroms.append(line.split()[0])
    sizes_set = set(sizes_chroms)

    data_chroms = {f[:-4] for f in os.listdir(split_dir) if f.endswith(".bed")}
    chrom_dropped = 0
    for chrom in sorted(data_chroms - sizes_set):
        chrom_bed = os.path.join(split_dir, f"{chrom}.bed")
        n_dropped = sum(1 for _ in open(chrom_bed))
        chrom_dropped += n_dropped
        print(f"[warning] Chrom '{chrom}' not in chrom sizes file — skipping ({n_dropped} rows dropped).", file=sys.stderr)

    # Stream chromosomes in sizes-file order through sort-bed, pipe to bgzip
    with open(out_path, "wb") as out_fh:
        bgzip_proc = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE, stdout=out_fh)

        for chrom in sizes_chroms:
            chrom_bed = os.path.join(split_dir, f"{chrom}.bed")
            if not os.path.isfile(chrom_bed):
                continue
            sort_proc = subprocess.run(
                ["sort-bed", chrom_bed],
                stdout=subprocess.PIPE,
                check=True
            )
            bgzip_proc.stdin.write(sort_proc.stdout)

        bgzip_proc.stdin.close()
        bgzip_proc.wait()

    if bgzip_proc.returncode != 0:
        print(f"[error] bgzip exited with code {bgzip_proc.returncode}", file=sys.stderr)
        sys.exit(1)

    # Clean up split dir
    import shutil
    shutil.rmtree(split_dir)

    # tabix index (fragment files: chr, start, end in cols 1-3)
    tabix_cmd = ["tabix", "-p", "bed", out_path]
    subprocess.run(tabix_cmd, check=True)

    print(f"[info] Done. Output: {out_path}", file=sys.stderr)
    print(f"[info] Index:  {out_path}.tbi", file=sys.stderr)

    os.unlink(tmp_path)


if __name__ == "__main__":
    main()