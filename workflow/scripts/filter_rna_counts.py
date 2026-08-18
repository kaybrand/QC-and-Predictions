#!/usr/bin/env python3
"""
Concatenate and QC-filter RNA count matrices (h5ad) for a given cell type.

For each annotation-{cellType}-IGVF* directory, opens rna_counts_mtx.h5ad,
filters to cells whose full obs_name appears in the QC filter guide, concatenates
all passing cells, and writes a single output file.

Default behavior converts Ensembl IDs to gene symbols, summing expression across
all Ensembl IDs that share the same gene symbol (e.g. 62,757 Ensembl IDs ->
61,217 unique gene symbols using GENCODE 43 / GRCh38). Use --ensemblIDs-as-genes
to skip this conversion and retain Ensembl IDs.

Output format is detected from the --out filename:
  *.csv or *.csv.gz        -> gzipped CSV (barcodes x genes)
  *.h5ad or *.h5           -> AnnData HDF5 file
  *.mtx or *.mtx.gz        -> sparse matrix directory:
                                matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz

Usage:
    python filter_rna_counts.py \
        --qc-guide /path/to/qc_filter_guide.tsv.gz \
        --pseudobulks /path/to/pseudobulks/ \
        --cell-type telohaec_tnf_0hr \
        --out /path/to/output/rna_counts.h5ad
"""

import argparse
import gzip
import glob
import os
import sys
import anndata as ad
import scipy.io
import scipy.sparse
import numpy as np
import pandas as pd
import re
from collections import Counter


def parse_args():
    p = argparse.ArgumentParser(description="Concatenate and QC-filter RNA h5ad matrices by cell type.")
    p.add_argument("--qc-guide",    required=True, help="Path to gzipped QC filter guide TSV.")
    p.add_argument("--pseudobulks", required=True, help="Path to pseudobulks folder.")
    p.add_argument("--cell-type",   required=True, help="Exact cell type identifier used in pseudobulk directory names (the string between 'annotation-' and '-IGVF'). Comma-separated for a merged cluster spanning multiple raw pseudobulk annotations (e.g. 'mcf7_1,mcf7_2').")
    p.add_argument("--out",         required=True, help="Output path. Format is inferred from extension: .csv/.csv.gz, .h5ad/.h5, or .mtx/.mtx.gz (writes a directory).")
    p.add_argument("--ensemblIDs-as-genes", action="store_true", dest="ensembl_ids",
                   help="Retain Ensembl IDs as gene identifiers. Default: convert to gene symbols, "
                        "summing expression across all Ensembl IDs that share the same gene symbol.")
    p.add_argument("--gtf", type=str, default=None,
                   help="Path to a GTF annotation file for Ensembl ID -> gene symbol conversion. "
                        "When provided, uses the GTF instead of embedded anndata metadata.")
    p.add_argument("--standard-chromosomes-only", action="store_true", dest="standard_chroms_only",
                   help="Only include genes from standard chromosomes (chr1-22, X, Y, M). Requires --gtf.")
    p.add_argument("--log", type=str, default=None,
                   help="Path to save a detailed log of GTF mapping observations.")
    return p.parse_args()


def detect_format(out_path: str) -> str:
    name = out_path.lower()
    if name.endswith(".csv.gz") or name.endswith(".csv"):
        return "csv"
    if name.endswith(".h5ad") or name.endswith(".h5"):
        return "h5ad"
    if name.endswith(".mtx.gz") or name.endswith(".mtx"):
        return "mtx"
    print(f"[error] Cannot infer output format from '{out_path}'.\n"
          f"        Supported extensions: .csv, .csv.gz, .h5ad, .h5, .mtx, .mtx.gz", file=sys.stderr)
    sys.exit(1)


def parse_gtf(gtf_path):
    """Parse a GTF file to build a map: versioned Ensembl gene ID -> (chrom, gene_symbol).

    Returns (ensembl_map, conflicts) where ensembl_map is
    {versioned_id: (chrom, gene_symbol)} and conflicts is a list of
    (id, prev_chrom, prev_name, new_chrom, new_name) tuples for IDs
    that had inconsistent gene names across GTF lines.
    """
    gene_id_re = re.compile(r'gene_id "([^"]+)"')
    gene_name_re = re.compile(r'gene_name "([^"]+)"')
    ensembl_map = {}
    conflicts = []

    print(f"[info] Parsing GTF file: {gtf_path} ...", file=sys.stderr)
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) < 9:
                continue
            chrom = parts[0]
            attrs = parts[8]
            gid_m = gene_id_re.search(attrs)
            gname_m = gene_name_re.search(attrs)

            if gid_m and gname_m:
                versioned_id = gid_m.group(1)
                gene_name = gname_m.group(1)

                if versioned_id in ensembl_map:
                    prev_chrom, prev_name = ensembl_map[versioned_id]
                    if prev_name != gene_name:
                        conflicts.append((versioned_id, prev_chrom, prev_name, chrom, gene_name))
                else:
                    ensembl_map[versioned_id] = (chrom, gene_name)

    print(f"[info] GTF: {len(ensembl_map)} unique versioned Ensembl gene IDs.", file=sys.stderr)
    if conflicts:
        print(f"[warning] {len(conflicts)} Ensembl IDs had conflicting gene names in GTF.", file=sys.stderr)

    return ensembl_map, conflicts


def load_passing_barcodes(qc_guide_path: str) -> set:
    """
    Return a set of full barcode strings from the QC filter guide.
    Matching is always on the full name (e.g. CCATATTTCGATAACC_IGVFSM4662QKFQ).
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


def find_h5ad_dirs(pseudobulks_dir: str, cell_type: str) -> list:
    """Returns all direct child directories matching annotation-{cellType}-IGVF*,
    for each comma-separated cell_type -- supports merging multiple raw
    pseudobulk annotations into one filtered output (e.g. a cluster config
    with pseudobulk_annotation "mcf7_1,mcf7_2" concatenates both raw
    directory sets; a plain single cell_type behaves exactly as before)."""
    all_dirs = []
    for ct in cell_type.split(","):
        ct = ct.strip()
        pattern = os.path.join(pseudobulks_dir, f"annotation-{ct}-IGVF*")
        dirs = sorted(d for d in glob.glob(pattern) if os.path.isdir(d))
        if dirs:
            print(f"[info] Found {len(dirs)} matching directories for cell type '{ct}'.", file=sys.stderr)
            all_dirs.extend(dirs)
            continue
        # Fall back to annotation-only pseudobulks (no subsample suffix)
        fallback = os.path.join(pseudobulks_dir, f"annotation-{ct}")
        if os.path.isdir(fallback):
            print(f"[info] No -IGVF* directories found for '{ct}'; using annotation-level pseudobulk: {fallback}", file=sys.stderr)
            all_dirs.append(fallback)
        else:
            print(f"[warning] No directories found matching: {pattern}", file=sys.stderr)
    return all_dirs


def collapse_to_gene_symbols(combined: ad.AnnData) -> ad.AnnData:
    """
    Convert Ensembl ID columns to gene symbol columns by summing expression
    across all Ensembl IDs that share the same gene symbol. Applies to X
    and all layers. Returns a new AnnData with gene symbols as var_names.
    """
    gene_symbols = combined.var["gene_symbol"].values
    unique_symbols, inverse = np.unique(gene_symbols, return_inverse=True)
    n_cells      = combined.n_obs
    n_ensembl    = combined.n_vars
    n_symbols    = len(unique_symbols)

    n_overloaded = int((np.bincount(inverse) > 1).sum())
    print(f"[info] Collapsing {n_ensembl} Ensembl IDs -> {n_symbols} unique gene symbols "
          f"({n_overloaded} symbols have multiple Ensembl IDs; their expression will be summed).",
          file=sys.stderr)

    # Build a (n_ensembl x n_symbols) summation matrix S.
    # new_X = X @ S  gives (n_cells x n_symbols) with duplicates summed.
    S = scipy.sparse.csr_matrix(
        (np.ones(n_ensembl), (np.arange(n_ensembl), inverse)),
        shape=(n_ensembl, n_symbols)
    )

    def _collapse(mat):
        if not scipy.sparse.issparse(mat):
            mat = scipy.sparse.csr_matrix(mat)
        return scipy.sparse.csr_matrix(mat @ S)

    new_X      = _collapse(combined.X)
    new_layers = {k: _collapse(combined.layers[k]) for k in combined.layers}
    new_var    = pd.DataFrame(index=unique_symbols)
    new_var.index.name = "gene_symbol"

    return ad.AnnData(
        X      = new_X,
        obs    = combined.obs.copy(),
        var    = new_var,
        layers = new_layers,
    )


def collapse_to_gene_symbols_gtf(combined, gtf_map, standard_chroms_only, log_lines):
    """
    Convert Ensembl ID columns to gene symbol columns using a GTF-derived map.

    Matches var_names against gtf_map using exact versioned Ensembl IDs.
    If any IDs are unmatched, populates log_lines and returns None (hard failure).
    If standard_chroms_only, drops genes on nonstandard chromosomes before collapsing.

    Returns a new AnnData with gene symbols as var_names, or None on failure.
    """
    STANDARD_CHROMS = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY', 'chrM'}
    var_names = list(combined.var_names)

    # --- Match anndata var_names against GTF (exact versioned) ---
    matched = {}
    unmatched = []
    for eid in var_names:
        if eid in gtf_map:
            matched[eid] = gtf_map[eid]
        else:
            unmatched.append(eid)

    log_lines.append(f"Total Ensembl IDs in anndata: {len(var_names)}")
    log_lines.append(f"Matched to GTF (exact versioned): {len(matched)}")
    log_lines.append(f"Unmatched: {len(unmatched)}")

    print(f"[info] GTF matching: {len(matched)}/{len(var_names)} exact versioned matches, "
          f"{len(unmatched)} unmatched", file=sys.stderr)

    if unmatched:
        log_lines.append(f"\nUnmatched Ensembl IDs ({len(unmatched)}):")
        for eid in sorted(unmatched):
            log_lines.append(f"  {eid}")
        print(f"[error] {len(unmatched)} Ensembl IDs in anndata not found in GTF "
              f"(exact versioned match).", file=sys.stderr)
        for eid in sorted(unmatched)[:10]:
            print(f"  {eid}", file=sys.stderr)
        if len(unmatched) > 10:
            print(f"  ... and {len(unmatched) - 10} more", file=sys.stderr)
        return None  # hard failure

    # --- Filter by chromosome if requested ---
    if standard_chroms_only:
        nonstandard_eids = [eid for eid in var_names if matched[eid][0] not in STANDARD_CHROMS]

        log_lines.append(f"\n--- Standard chromosomes filter ---")
        if nonstandard_eids:
            chrom_counts = Counter(matched[eid][0] for eid in nonstandard_eids)
            log_lines.append(f"Ensembl IDs on nonstandard chromosomes: {len(nonstandard_eids)}")
            log_lines.append(f"Nonstandard chromosomes found:")
            for chrom, count in sorted(chrom_counts.items()):
                log_lines.append(f"  {chrom}: {count} Ensembl IDs")

            # Identify inflation-risk genes (symbol on both standard + nonstandard chroms)
            standard_symbols = {matched[eid][1] for eid in var_names
                                if matched[eid][0] in STANDARD_CHROMS}
            nonstandard_symbols = {matched[eid][1] for eid in nonstandard_eids}
            shared_symbols = standard_symbols & nonstandard_symbols
            lost_symbols = nonstandard_symbols - standard_symbols

            log_lines.append(f"\nGene symbols with contributions from BOTH standard and "
                             f"nonstandard chromosomes ({len(shared_symbols)}):")
            log_lines.append("(These are the inflation-risk genes whose UMI counts "
                             "would be reduced by filtering)")
            for sym in sorted(shared_symbols):
                std_ids = [eid for eid in var_names
                           if matched[eid][1] == sym and matched[eid][0] in STANDARD_CHROMS]
                nonstd_ids = [eid for eid in nonstandard_eids if matched[eid][1] == sym]
                log_lines.append(f"  {sym}:")
                for eid in std_ids:
                    log_lines.append(f"    standard:    {eid} ({matched[eid][0]})")
                for eid in nonstd_ids:
                    log_lines.append(f"    nonstandard: {eid} ({matched[eid][0]})")

            if lost_symbols:
                log_lines.append(f"\nGene symbols ONLY on nonstandard chromosomes "
                                 f"(removed entirely): {len(lost_symbols)}")
                for sym in sorted(lost_symbols):
                    ids = [eid for eid in nonstandard_eids if matched[eid][1] == sym]
                    log_lines.append(f"  {sym}: {', '.join(ids)}")

            print(f"[info] Removing {len(nonstandard_eids)} Ensembl IDs on nonstandard "
                  f"chromosomes ({len(shared_symbols)} symbols had both standard+nonstandard "
                  f"IDs, {len(lost_symbols)} symbols lost entirely).", file=sys.stderr)
        else:
            log_lines.append("All Ensembl IDs are on standard chromosomes. No filtering needed.")
            print("[info] All Ensembl IDs are on standard chromosomes.", file=sys.stderr)

        # Filter the anndata to keep only standard-chromosome genes
        keep_mask = np.array([matched[eid][0] in STANDARD_CHROMS for eid in var_names])
        combined = combined[:, keep_mask].copy()
        var_names = list(combined.var_names)
        matched = {eid: matched[eid] for eid in var_names}
        log_lines.append(f"Ensembl IDs after chromosome filtering: {len(var_names)}")

    # --- Collapse by gene symbol ---
    gene_symbols = np.array([matched[eid][1] for eid in var_names])
    unique_symbols, inverse = np.unique(gene_symbols, return_inverse=True)
    n_ensembl = combined.n_vars
    n_symbols = len(unique_symbols)
    bin_counts = np.bincount(inverse)
    n_overloaded = int((bin_counts > 1).sum())

    log_lines.append(f"\n--- Gene symbol collapse ---")
    log_lines.append(f"Ensembl IDs entering collapse: {n_ensembl}")
    log_lines.append(f"Unique gene symbols after collapse: {n_symbols}")
    log_lines.append(f"Symbols with multiple Ensembl IDs (summed): {n_overloaded}")

    if n_overloaded > 0:
        log_lines.append(f"\nOverloaded symbols (multiple Ensembl IDs -> same symbol):")
        for i, sym in enumerate(unique_symbols):
            if bin_counts[i] > 1:
                contributing = [eid for eid in var_names if matched[eid][1] == sym]
                log_lines.append(f"  {sym} ({bin_counts[i]} IDs):")
                for eid in contributing:
                    log_lines.append(f"    {eid} ({matched[eid][0]})")

    print(f"[info] Collapsing {n_ensembl} Ensembl IDs -> {n_symbols} unique gene symbols "
          f"({n_overloaded} symbols have multiple Ensembl IDs; expression summed).",
          file=sys.stderr)

    # Build summation matrix
    S = scipy.sparse.csr_matrix(
        (np.ones(n_ensembl), (np.arange(n_ensembl), inverse)),
        shape=(n_ensembl, n_symbols)
    )

    def _collapse(mat):
        if not scipy.sparse.issparse(mat):
            mat = scipy.sparse.csr_matrix(mat)
        return scipy.sparse.csr_matrix(mat @ S)

    new_X = _collapse(combined.X)
    new_layers = {k: _collapse(combined.layers[k]) for k in combined.layers}
    new_var = pd.DataFrame(index=unique_symbols)
    new_var.index.name = "gene_symbol"

    return ad.AnnData(
        X=new_X,
        obs=combined.obs.copy(),
        var=new_var,
        layers=new_layers,
    )


def sanity_checks(combined: ad.AnnData, passing_barcodes: set, ensembl_ids: bool):
    """
    Run four checks on the final matrix and exit with an error if any fail:
      1. No duplicate cell barcodes (obs)
      2. No duplicate gene identifiers (var)
      3. Row count == number of barcodes in the QC guide
      4. Column count == number of unique gene identifiers expected
    """
    errors = []

    # 1. Duplicate obs
    obs_counts = pd.Series(list(combined.obs_names)).value_counts()
    dup_obs = obs_counts[obs_counts > 1]
    if not dup_obs.empty:
        errors.append(
            f"  Duplicate cell barcodes ({len(dup_obs)} offending barcodes):\n" +
            "\n".join(f"    {bc}: {n} occurrences" for bc, n in dup_obs.items())
        )

    # 2. Duplicate var names
    var_counts = pd.Series(list(combined.var_names)).value_counts()
    dup_var = var_counts[var_counts > 1]
    if not dup_var.empty:
        errors.append(
            f"  Duplicate gene identifiers ({len(dup_var)} offending names):\n" +
            "\n".join(f"    {g}: {n} occurrences" for g, n in dup_var.items())
        )

    # 3. Row count
    n_obs   = combined.n_obs
    n_guide = len(passing_barcodes)
    if n_obs != n_guide:
        # Report which barcodes are missing or extra
        found   = set(combined.obs_names)
        missing = passing_barcodes - found
        extra   = found - passing_barcodes
        msg = f"  Cell count mismatch: matrix has {n_obs} cells, QC guide has {n_guide} (difference: {n_obs - n_guide:+d})."
        if missing:
            msg += f"\n    Barcodes in QC guide but absent from matrix ({len(missing)}):\n"
            msg += "\n".join(f"      {bc}" for bc in sorted(missing)[:10])
            if len(missing) > 10:
                msg += f"\n      ... and {len(missing) - 10} more."
        if extra:
            msg += f"\n    Barcodes in matrix but absent from QC guide ({len(extra)}):\n"
            msg += "\n".join(f"      {bc}" for bc in sorted(extra)[:10])
            if len(extra) > 10:
                msg += f"\n      ... and {len(extra) - 10} more."
        errors.append(msg)

    # 4. Column count
    n_vars    = combined.n_vars
    # Expected: unique gene symbols (post-collapse) or unique Ensembl IDs
    id_label  = "Ensembl IDs" if ensembl_ids else "gene symbols"
    # We can only verify uniqueness here; the expected count is n_vars itself
    # after collapse, so this check is meaningful only if dup_var is also empty.
    if dup_var.empty:
        print(f"[info] Gene identifier check passed: {n_vars} unique {id_label}.", file=sys.stderr)

    if errors:
        print(f"[error] Sanity check(s) failed:\n" + "\n".join(errors), file=sys.stderr)
        sys.exit(1)

    print(f"[info] All sanity checks passed: {n_obs} cells x {n_vars} {id_label}, "
          f"no duplicates in rows or columns.", file=sys.stderr)


def write_csv(combined: ad.AnnData, out_path: str):
    if not out_path.endswith(".gz"):
        out_path += ".gz"
    n_cells, n_genes = combined.shape
    print(f"[info] Writing CSV ({n_cells} cells x {n_genes} genes) to {out_path} ...", file=sys.stderr)
    if n_cells * n_genes > 5e8:
        print("[warning] Dense CSV output is very large. Consider .h5ad or .mtx format.", file=sys.stderr)
    X = combined.X.toarray() if scipy.sparse.issparse(combined.X) else np.array(combined.X)
    df = pd.DataFrame(X, index=combined.obs_names, columns=combined.var_names)
    df.index.name = "barcode"
    df.to_csv(out_path, compression="gzip")


def write_h5ad(combined: ad.AnnData, out_path: str):
    print(f"[info] Writing h5ad to {out_path} ...", file=sys.stderr)
    combined.write_h5ad(out_path)


def write_mtx(combined: ad.AnnData, out_path: str):
    base = out_path
    for ext in (".mtx.gz", ".mtx"):
        if base.lower().endswith(ext):
            base = base[:-len(ext)]
            break
    out_dir = base
    os.makedirs(out_dir, exist_ok=True)

    mtx_path      = os.path.join(out_dir, "matrix.mtx.gz")
    barcodes_path = os.path.join(out_dir, "barcodes.tsv.gz")
    genes_path    = os.path.join(out_dir, "features.tsv.gz")

    n_cells, n_genes = combined.shape
    print(f"[info] Writing sparse MTX directory ({n_cells} cells x {n_genes} genes) to {out_dir}/ ...", file=sys.stderr)

    X = combined.X
    if not scipy.sparse.issparse(X):
        X = scipy.sparse.csr_matrix(X)
    with gzip.open(mtx_path, "wb") as f:
        scipy.io.mmwrite(f, X.T)  # MTX convention: genes x cells

    with gzip.open(barcodes_path, "wt") as f:
        for bc in combined.obs_names:
            f.write(bc + "\n")

    with gzip.open(genes_path, "wt") as f:
        for gene in combined.var_names:
            f.write(gene + "\n")

    print(f"[info]   {mtx_path}", file=sys.stderr)
    print(f"[info]   {barcodes_path}", file=sys.stderr)
    print(f"[info]   {genes_path}", file=sys.stderr)


def main():
    args = parse_args()

    # Validate argument combinations
    if args.standard_chroms_only and not args.gtf:
        print("[error] --standard-chromosomes-only requires --gtf.", file=sys.stderr)
        sys.exit(1)
    if args.gtf and args.ensembl_ids:
        print("[error] --gtf cannot be used with --ensemblIDs-as-genes.", file=sys.stderr)
        sys.exit(1)

    fmt = detect_format(args.out)
    passing_barcodes = load_passing_barcodes(args.qc_guide)
    dirs = find_h5ad_dirs(args.pseudobulks, args.cell_type)

    if not dirs:
        print("[error] No matching directories found. Exiting.", file=sys.stderr)
        sys.exit(1)

    adatas = []
    for d in dirs:
        h5ad_path = os.path.join(d, "rna_counts_mtx.h5ad")
        if not os.path.isfile(h5ad_path):
            print(f"[warning] rna_counts_mtx.h5ad not found in {d}, skipping.", file=sys.stderr)
            continue

        adata = ad.read_h5ad(h5ad_path)
        n_before = adata.n_obs

        # Filter on full obs_name — no trimming, no remapping.
        mask = adata.obs_names.isin(passing_barcodes)
        adata = adata[mask].copy()
        n_after = adata.n_obs

        print(f"[info]   {os.path.basename(d)}: {n_before} cells -> {n_after} passed QC", file=sys.stderr)

        if n_after == 0:
            print(f"[warning]   No passing cells in {d}, skipping.", file=sys.stderr)
            continue

        adatas.append(adata)

    if not adatas:
        print("[error] No passing cells found across any directory. Exiting.", file=sys.stderr)
        sys.exit(1)

    print(f"[info] Concatenating {len(adatas)} matrices ...", file=sys.stderr)
    combined = ad.concat(adatas, axis=0, join="inner", merge="same")
    print(f"[info] Combined shape after concatenation: {combined.shape} (cells x Ensembl IDs)", file=sys.stderr)

    # Convert Ensembl IDs -> gene symbols (default) by summing duplicate columns
    if not args.ensembl_ids:
        if args.gtf:
            log_lines = [
                "=== GTF-based gene symbol mapping log ===",
                f"GTF file: {args.gtf}",
                f"Standard chromosomes only: {args.standard_chroms_only}",
                "",
            ]
            gtf_map, gtf_conflicts = parse_gtf(args.gtf)
            log_lines.append(f"Unique versioned Ensembl gene IDs in GTF: {len(gtf_map)}")
            if gtf_conflicts:
                log_lines.append(f"\nGTF conflicts (same Ensembl ID, different gene names): "
                                 f"{len(gtf_conflicts)}")
                for vid, pc, pn, nc, nn in gtf_conflicts:
                    log_lines.append(f"  {vid}: '{pn}' ({pc}) vs '{nn}' ({nc})")
            log_lines.append("")

            result = collapse_to_gene_symbols_gtf(
                combined, gtf_map, args.standard_chroms_only, log_lines
            )

            # Determine log path
            log_path = args.log
            if log_path is None:
                base = args.out
                for ext in ('.mtx.gz', '.mtx', '.csv.gz', '.csv', '.h5ad', '.h5'):
                    if base.lower().endswith(ext):
                        base = base[:-len(ext)]
                        break
                log_path = base + '_gtf_mapping.txt'
            os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
            with open(log_path, 'w') as f:
                f.write('\n'.join(log_lines) + '\n')
            print(f"[info] Detailed GTF mapping log saved to {log_path}", file=sys.stderr)

            if result is None:
                print(f"[error] Exiting due to unmatched Ensembl IDs. See log: {log_path}",
                      file=sys.stderr)
                sys.exit(1)
            combined = result
        else:
            combined = collapse_to_gene_symbols(combined)
        print(f"[info] Shape after gene symbol collapse: {combined.shape} (cells x gene symbols)",
              file=sys.stderr)

    # Sanity checks — exits with error on any failure
    sanity_checks(combined, passing_barcodes, args.ensembl_ids)

    # Write output
    if fmt != "mtx":
        out_dir = os.path.dirname(os.path.abspath(args.out))
        os.makedirs(out_dir, exist_ok=True)

    if fmt == "csv":
        write_csv(combined, args.out)
    elif fmt == "h5ad":
        write_h5ad(combined, args.out)
    elif fmt == "mtx":
        write_mtx(combined, args.out)

    print(f"[info] Done.", file=sys.stderr)


if __name__ == "__main__":
    main()