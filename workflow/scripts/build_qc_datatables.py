#!/usr/bin/env python3
"""Build per-cluster per-cell QC datatables from downloaded pseudobulks.

Concatenates each pseudobulk's per-subsample per_cell_qc file into one table per
cluster, reproducing the layout the QC tooling already expects:

    {out-dir}/QC_datatables/{dataset}_data/{cluster}_per_cell_qc.tsv

That flat convention -- no {cluster}/ subdirectory, the underscore is part of
the filename -- matches the 179 files already under
QC_pseudobulks/datatables/ and is documented in
QC_pseudobulks/scripts/plotting_scripts/CELL_UMAP_QC_README.md:22-25. It is what
the ~13 QC_pseudobulks scripts glob and what resolve_exclusions.py builds a path
to, so keeping it byte-compatible means nothing downstream needs to change.

This step had no implementation anywhere before now. The transformation was
recovered by comparing the existing outputs against their inputs byte for byte
(igvf0/wtc11_endo_early_mesoderm: 1677 data rows across 16 subsample
directories; igvf18/mcf7, a merged cluster) and is exactly:

    header of the first source file, then every subsequent file's rows minus its
    header, with subsample directories visited in SORTED order.

Nothing is added or rewritten -- notably NOT a subsample column, because the raw
per-subsample file already carries one, populated with that subsample's own
accession. The 12 columns are:
    analysis_accession, barcode, subsample, rna_read_count, gene_count,
    pct_mito, pct_ribo, num_frags, pct_duplicated_reads, nucleosomal_signal,
    tss_enrichment, frip

The output is UNCOMPRESSED .tsv, matching the existing files, even though the
downloaded inputs are per_cell_qc.tsv.gz -- the portal gzips this file where the
old archive did not, so reading gzip is mandatory here.

Header agreement across a cluster's sources is asserted, not assumed: a silent
upstream schema change should stop the run rather than produce a table whose
rows disagree about what its columns mean.

Merged clusters: a cluster whose pipeline config gives a comma-separated
`pseudobulk_annotation` (e.g. igvf18's mcf7 = "mcf7_1,mcf7_2") produces ONE
table spanning both annotations, keyed by the cluster name -- matching both the
file on disk today and resolve_exclusions.py:82-88's deliberate choice to key on
`cluster` rather than `pseudobulk_annotation`. Annotations are concatenated in
the order the config lists them. Without a config, cluster == annotation, which
is the right default for a newly discovered pseudobulk that has no config yet.
"""

import argparse
import glob
import gzip
import os
import re
import sys

DIRNAME_RE = re.compile(r"^annotation-(?P<annotation>.+)-(?P<subsample>IGVFSM\w+)$")
PER_CELL_QC_NAMES = ("per_cell_qc.tsv.gz", "per_cell_qc.tsv")
EXPECTED_COLUMNS = [
    "analysis_accession", "barcode", "subsample", "rna_read_count", "gene_count",
    "pct_mito", "pct_ribo", "num_frags", "pct_duplicated_reads",
    "nucleosomal_signal", "tss_enrichment", "frip",
]


def log(msg):
    print(f"[build_qc_datatables] {msg}", file=sys.stderr)


def _open(path):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, "r", newline="")


def per_cell_qc_path(directory):
    """The per-cell QC file in a pseudobulk directory, gz or plain. Returns None
    if absent -- an incomplete pseudobulk (some sets genuinely lack this file
    yet) must be reported, not crash the build."""
    for name in PER_CELL_QC_NAMES:
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            return candidate
    return None


def annotation_dirs(pseudobulks_root, dataset, annotation):
    """This annotation's pseudobulk directories, sorted -- the order that
    reproduces the existing datatables exactly."""
    pattern = os.path.join(pseudobulks_root, dataset, "pseudobulks", f"annotation-{annotation}-IGVFSM*")
    return sorted(d for d in glob.glob(pattern) if os.path.isdir(d))


def discover_annotations(pseudobulks_root, dataset):
    base = os.path.join(pseudobulks_root, dataset, "pseudobulks")
    if not os.path.isdir(base):
        return []
    found = set()
    for entry in os.listdir(base):
        m = DIRNAME_RE.match(entry)
        if m and os.path.isdir(os.path.join(base, entry)):
            found.add(m.group("annotation"))
    return sorted(found)


def build_one(pseudobulks_root, dataset, cluster, annotations, out_path):
    """Concatenate one cluster's sources. Returns (n_rows, n_sources, skipped)."""
    sources = []
    skipped = []
    for annotation in annotations:
        dirs = annotation_dirs(pseudobulks_root, dataset, annotation)
        if not dirs:
            skipped.append(f"{annotation}:no_directories")
            continue
        for directory in dirs:
            path = per_cell_qc_path(directory)
            if path is None:
                skipped.append(f"{os.path.basename(directory)}:no_per_cell_qc")
                continue
            sources.append(path)
    if not sources:
        return 0, 0, skipped

    header = None
    total = 0
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "w", newline="") as out:
        for path in sources:
            with _open(path) as fh:
                first = fh.readline()
                if header is None:
                    header = first
                    cols = header.rstrip("\n").split("\t")
                    if cols != EXPECTED_COLUMNS:
                        os.unlink(tmp)
                        raise ValueError(
                            f"{path}: unexpected columns.\n  expected {EXPECTED_COLUMNS}\n  got      {cols}"
                        )
                    out.write(header)
                elif first != header:
                    # Refuse rather than emit a table whose rows disagree about
                    # what its columns mean.
                    os.unlink(tmp)
                    raise ValueError(
                        f"{path}: header differs from the first source in this cluster.\n"
                        f"  first : {header.rstrip()}\n  this  : {first.rstrip()}"
                    )
                for line in fh:
                    out.write(line)
                    total += 1
    os.replace(tmp, out_path)
    return total, len(sources), skipped


def load_cluster_map(config_paths):
    """{dataset: {cluster: [annotation, ...]}} from pipeline configs, so merged
    clusters produce a single table. Parsed with a small reader rather than PyYAML
    because the interpreter that has requests+igvf_utils has no yaml, and the
    generated `clusters:` block is a fixed two-level shape."""
    mapping = {}
    for path in config_paths:
        dataset = cluster = None
        in_clusters = False
        with open(path) as fh:
            for raw in fh:
                line = raw.rstrip("\n")
                if not line.strip() or line.lstrip().startswith("#"):
                    continue
                if re.match(r"^\S", line):
                    in_clusters = line.startswith("clusters:")
                    dataset = cluster = None
                    continue
                if not in_clusters:
                    continue
                m = re.match(r"^  (\S+):\s*$", line)
                if m:
                    dataset, cluster = m.group(1), None
                    mapping.setdefault(dataset, {})
                    continue
                m = re.match(r"^    (\S+):\s*$", line)
                if m and dataset:
                    cluster = m.group(1)
                    mapping[dataset].setdefault(cluster, [cluster])
                    continue
                m = re.match(r"^      pseudobulk_annotation:\s*(.+?)\s*$", line)
                if m and dataset and cluster:
                    mapping[dataset][cluster] = [a.strip() for a in m.group(1).split(",") if a.strip()]
    return mapping


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--pseudobulks-root", required=True, help="root of the downloaded pseudobulk tree")
    p.add_argument("--out-dir", required=True, help="QC_datatables/ is created under here")
    p.add_argument("--datasets", default="", help="comma-separated datasets (default: all found)")
    p.add_argument(
        "--pipeline-config",
        action="append",
        default=[],
        help="a *_pipeline_config.yaml, repeatable. Supplies cluster -> annotation(s) so merged "
        "clusters yield one table. Without any, cluster == annotation.",
    )
    p.add_argument("--force", action="store_true", help="rebuild tables that already exist")
    args = p.parse_args(argv)

    cluster_map = load_cluster_map(args.pipeline_config)
    if cluster_map:
        merged = {
            (d, c): a for d, cs in cluster_map.items() for c, a in cs.items() if len(a) > 1
        }
        log(f"loaded cluster map for {len(cluster_map)} dataset(s); {len(merged)} merged cluster(s): {sorted(merged)}")

    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    if not datasets:
        root = args.pseudobulks_root
        datasets = sorted(
            d for d in os.listdir(root)
            if os.path.isdir(os.path.join(root, d, "pseudobulks"))
        )
    log(f"datasets: {datasets}")

    built = failed = 0
    for dataset in datasets:
        available = set(discover_annotations(args.pseudobulks_root, dataset))
        if not available:
            log(f"{dataset}: no pseudobulk directories, skipping")
            continue
        # Prefer the config's cluster grouping; fall back to one table per
        # annotation for anything the config doesn't mention (a new pseudobulk).
        clusters = {}
        for cluster, annotations in (cluster_map.get(dataset) or {}).items():
            if set(annotations) & available:
                clusters[cluster] = annotations
        covered = {a for anns in clusters.values() for a in anns}
        for annotation in sorted(available - covered):
            clusters[annotation] = [annotation]

        out_dir = os.path.join(args.out_dir, "QC_datatables", f"{dataset}_data")
        for cluster, annotations in sorted(clusters.items()):
            out_path = os.path.join(out_dir, f"{cluster}_per_cell_qc.tsv")
            if os.path.exists(out_path) and not args.force:
                log(f"{dataset}/{cluster}: exists, skipping (use --force)")
                continue
            try:
                rows, n_sources, skipped = build_one(
                    args.pseudobulks_root, dataset, cluster, annotations, out_path
                )
            except ValueError as exc:
                log(f"{dataset}/{cluster}: FAILED -- {exc}")
                failed += 1
                continue
            if not n_sources:
                log(f"{dataset}/{cluster}: no per-cell QC sources ({'; '.join(skipped) or 'none found'})")
                continue
            built += 1
            note = f" (skipped: {'; '.join(skipped)})" if skipped else ""
            log(f"{dataset}/{cluster}: {rows} rows from {n_sources} subsample file(s){note}")

    log(f"built {built} datatable(s), {failed} failure(s)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
