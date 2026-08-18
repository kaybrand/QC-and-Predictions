"""
Parse-time exclusion resolution for the QC-and-Predictions pipeline.

Computes cell_count / fragments_total / umi_count for each configured
(dataset, cluster) pair from files that already exist BEFORE this pipeline's
own rules ever run (the QC guide, and either the pre-existing
`filtered_cell_subsample_metrics.tsv` or a merged per-cell-QC table), so
exclusion decisions can be made before any rule is defined -- no Snakemake
checkpoint required, and excluded clusters never get a single rule
instantiated for them.

A cluster's true unique identity is the (dataset, cluster) pair, not
`cluster` alone -- cluster names are not guaranteed unique across datasets,
only within one dataset's own config/cell_clusters table. Every set/dict
here is keyed by (dataset, cluster) tuples for that reason.

See the "Prefilter stats" / "Multi-dataset support" sections of the design
plan for the reasoning and the telohaec_crispri validation numbers this
approach was checked against.
"""

import csv
import gzip
import os

DEFAULT_QC_GUIDE_NAME = "filtered_barcodes_with_subsamples.tsv.gz"


def _read_qc_guide_barcodes(qc_guide_path):
    """Return the set of barcodes listed in a QC guide TSV."""
    with gzip.open(qc_guide_path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return {row["barcode"] for row in reader}


def _stats_from_subsample_metrics(metrics_path):
    """Sum n_cells/total_fragments/total_RNA_reads across all subsample rows."""
    cell_count = fragments_total = umi_count = 0
    with open(metrics_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cell_count += int(row["n_cells"])
            fragments_total += int(row["total_fragments"])
            umi_count += int(row["total_RNA_reads"])
    return cell_count, fragments_total, umi_count


def _stats_from_per_cell_qc_join(qc_guide_path, per_cell_qc_path):
    """Fallback for a non-default QC guide: join barcodes against the merged
    per-cluster per_cell_qc.tsv table and sum num_frags/rna_read_count."""
    barcodes = _read_qc_guide_barcodes(qc_guide_path)
    cell_count = fragments_total = umi_count = 0
    with open(per_cell_qc_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["barcode"] in barcodes:
                cell_count += 1
                fragments_total += int(row["num_frags"])
                umi_count += int(row["rna_read_count"])
    return cell_count, fragments_total, umi_count


def compute_cluster_stats(plots_dir, datatables_dir, dataset, cluster, pseudobulk_annotation, qc_guide_path, has_rna):
    """
    Returns (stats, reason). stats is a dict {cell_count, fragments_total,
    umi_count}, or None if the required source files don't exist yet.
    reason is "ok" on success, else "missing_qc_guide" (no QC guide file at
    all) or "missing_per_cell_qc_table" (QC guide exists, but the joined
    per-cell-QC table a non-default guide needs doesn't) -- distinct failure
    modes a caller (e.g. the coverage report) needs to tell apart.

    umi_count is None (not zero) for ATAC-only clusters -- the min_umi_count
    threshold is skipped for them, not failed.
    """
    qc_guide_name = os.path.basename(qc_guide_path)
    cluster_plots_dir = os.path.join(plots_dir, dataset, cluster)
    metrics_path = os.path.join(cluster_plots_dir, "filtered_cell_subsample_metrics.tsv")

    if qc_guide_name == DEFAULT_QC_GUIDE_NAME and os.path.exists(metrics_path):
        cell_count, fragments_total, umi_count = _stats_from_subsample_metrics(metrics_path)
    else:
        if not os.path.exists(qc_guide_path):
            return None, "missing_qc_guide"
        # Keyed by `cluster` (the cluster's own identity), not
        # pseudobulk_annotation -- for a merged cluster (comma-separated
        # pseudobulk_annotation, e.g. igvf18's mcf7 = mcf7_1 + mcf7_2) this
        # expects a single pre-built mcf7_per_cell_qc.tsv, not one per raw
        # source name. For every other cluster `cluster` == `pseudobulk_annotation`,
        # so this is unchanged from before.
        per_cell_qc_path = os.path.join(datatables_dir, f"{dataset}_data", f"{cluster}_per_cell_qc.tsv")
        if not os.path.exists(per_cell_qc_path):
            return None, "missing_per_cell_qc_table"
        cell_count, fragments_total, umi_count = _stats_from_per_cell_qc_join(qc_guide_path, per_cell_qc_path)

    stats = {
        "cell_count": cell_count,
        "fragments_total": fragments_total,
        "umi_count": umi_count if has_rna else None,
    }
    return stats, "ok"


def _iter_cluster_configs(clusters_by_dataset):
    """Yields (dataset, cluster, cluster_cfg) for every configured cluster, across all datasets."""
    for dataset, dataset_clusters in clusters_by_dataset.items():
        for cluster, cluster_cfg in dataset_clusters.items():
            yield dataset, cluster, cluster_cfg


def _parse_qualified_names(qualified_names):
    """Parses "dataset/cluster" strings (as used in exclusion.user_specified) into (dataset, cluster) tuples."""
    parsed = set()
    for name in qualified_names:
        dataset, _, cluster = name.partition("/")
        parsed.add((dataset, cluster))
    return parsed


def resolve_exclusions(config, data_dir):
    """
    Returns (included_clusters, upload_eligible_clusters, excluded_clusters, stats_by_cluster),
    all keyed by (dataset, cluster) tuples.

    included_clusters:        (dataset, cluster) pairs this run will actually process (filter/scE2G/etc)
    upload_eligible_clusters: subset of included_clusters whose products may be uploaded
    excluded_clusters:        (dataset, cluster) pairs this run will skip entirely (unless process_excluded_no_upload)

    data_dir: config["data_dir"] -- corrected 2026-08-03, this used to be the
    pipeline's own code-checkout root (WDIR), which only found real
    plots/datatables when the code and the data happened to live in the same
    place.
    """
    plots_dir = os.path.join(data_dir, "plots")
    datatables_dir = os.path.join(data_dir, "datatables")

    exclusion_cfg = config.get("exclusion", {})
    user_specified = _parse_qualified_names(exclusion_cfg.get("user_specified", []))
    process_excluded_no_upload = exclusion_cfg.get("process_excluded_no_upload", False)
    thresholds = exclusion_cfg.get("auto_thresholds", {})
    min_cell_count = thresholds.get("min_cell_count", 0)
    min_fragments_total = thresholds.get("min_fragments_total", 0)
    min_umi_count = thresholds.get("min_umi_count", 0)

    excluded_clusters = set()
    stats_by_cluster = {}
    all_clusters = set()

    for dataset, cluster, cluster_cfg in _iter_cluster_configs(config["clusters"]):
        key = (dataset, cluster)
        all_clusters.add(key)
        has_rna = cluster_cfg["models"] != ["scATAC_powerlaw_v3"]
        stats, stat_reason = compute_cluster_stats(
            plots_dir, datatables_dir, dataset, cluster,
            cluster_cfg["pseudobulk_annotation"], cluster_cfg["qc_guide"], has_rna,
        )

        if key in user_specified:
            excluded_clusters.add(key)
            reason = "user_specified"
        elif stats is None:
            # No stats available yet (brand-new cluster, or missing QC-guide
            # inputs) -- exclude it. A cluster with no QC guide on disk can't
            # mechanically be processed at all (atac_fragment_file/
            # rna_count_matrix both require it as an input Snakemake can't
            # produce), so this is a real exclusion, not just an unresolved
            # quality check -- distinguished from quality-threshold failures
            # via `reason` (missing_qc_guide/missing_per_cell_qc_table vs
            # below_min_*), which the coverage report keys off of.
            excluded_clusters.add(key)
            reason = stat_reason
        elif stats["cell_count"] < min_cell_count:
            excluded_clusters.add(key)
            reason = "below_min_cell_count"
        elif stats["fragments_total"] < min_fragments_total:
            excluded_clusters.add(key)
            reason = "below_min_fragments_total"
        elif stats["umi_count"] is not None and stats["umi_count"] < min_umi_count:
            excluded_clusters.add(key)
            reason = "below_min_umi_count"
        else:
            reason = "pass"

        stats_by_cluster[key] = {
            "cell_count": stats["cell_count"] if stats else None,
            "fragments_total": stats["fragments_total"] if stats else None,
            "umi_count": stats["umi_count"] if stats else None,
            "reason": reason,
        }

    if process_excluded_no_upload:
        included_clusters = all_clusters
    else:
        included_clusters = all_clusters - excluded_clusters
    upload_eligible_clusters = included_clusters - excluded_clusters

    return included_clusters, upload_eligible_clusters, excluded_clusters, stats_by_cluster


CLUSTER_STATS_HEADER = ["dataset", "cluster", "cell_count", "fragments_total", "umi_count", "reason"]


def write_cluster_stats_table(dataset, stats_by_cluster, out_dir):
    """Deterministic full overwrite per dataset -- unlike write_scE2G_config.py's
    tables, there's no cross-invocation accumulation to merge, so a plain
    temp-file + os.replace is enough (no flock needed): every worker node that
    re-parses this Snakefile computes byte-identical content for the same
    dataset, so the only real hazard is a torn read, which os.replace avoids."""
    path = os.path.join(out_dir, f"{dataset}_cluster_stats.tsv")
    os.makedirs(out_dir, exist_ok=True)
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CLUSTER_STATS_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (ds, cluster), row in sorted(stats_by_cluster.items()):
            if ds == dataset:
                writer.writerow({"dataset": ds, "cluster": cluster, **row})
    os.replace(tmp, path)
