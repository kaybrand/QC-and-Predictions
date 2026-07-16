"""
Regenerates the per-dataset QC sanity-check artifacts scE2G itself produces
(all_qc_stats.tsv, its 6 summary PDFs, and predictions_qc_report.html).
These are never requested by this pipeline's own `rule all` unless
explicitly added here, since we don't reuse scE2G's own `rule all` (which is
what normally requests them).

`aggregate_qc_stats` (Python) and `plot_qc_stats` (forked R script) REPLACE
scE2G's own `plot_stats` rule, which is excluded when importing scE2G (see
common.smk) -- that rule rebuilds all_qc_stats.tsv strictly from whichever
clusters are in THIS pipeline's own cell_clusters.tsv for the current run,
so it can't reflect clusters processed some other way or preserve clusters
untouched in a given invocation. Our replacement writes to the exact same
output path (RESULTS_DIR/qc_plots/all_qc_stats.tsv), so scE2G's own
(unexcluded, still-imported) `hover_plots` rule -- whose only input is that
same path -- transparently depends on our replacement instead of scE2G's
original, without any changes needed on its end; we just target its output
file directly, by path, further below.

Coverage is dataset-wide, not run-scoped: all_qc_stats.tsv always covers
every cluster directory that exists on disk for a dataset (see
aggregate_qc_stats.py), including clusters this particular invocation never
touched -- this is a sanity-check view and must show a low-quality cluster
(e.g. one with 15 cells) sitting below every threshold line for comparison,
not silently omit it just because it wasn't part of this run's target set.
"""

QC_PLOTS_DIR = os.path.join(RESULTS_DIR_BASE, "{dataset}", "qc_plots")
SCE2G_SC_E2G_ENV = os.path.join(config["scE2G_dir"], "workflow", "envs", "sc_e2g.yml")


def per_cluster_stats_files(dataset):
    """This run's freshly-(re)generated stats.tsv paths for `dataset` -- declared
    as explicit Snakemake inputs purely so aggregate_qc_stats waits for them to be
    current before it re-scans the dataset directory (see aggregate_qc_stats.py's
    own docstring for why the scan itself needs no such list)."""
    files = []
    for ds, cluster in INCLUDED_CLUSTERS:
        if ds != dataset:
            continue
        for model in config["clusters"][dataset][cluster]["models"]:
            threshold = get_model_threshold(config["scE2G_dir"], model)
            files.append(os.path.join(
                RESULTS_DIR_BASE, dataset, cluster, model,
                f"scE2G_predictions_threshold{threshold}_stats.tsv",
            ))
    return files


def get_qc_stats_targets():
    targets = []
    for dataset in DATASETS:
        plots_dir = QC_PLOTS_DIR.format(dataset=dataset)
        targets.append(os.path.join(plots_dir, "all_qc_stats.tsv"))
        targets.append(os.path.join(plots_dir, "predictions_qc_report.html"))
        targets += [
            os.path.join(plots_dir, name)
            for name in (
                "dataset_metrics.pdf", "enhancer_metrics.pdf", "eg_metrics.pdf",
                "gene_metrics.pdf", "distance_and_size_metrics.pdf", "qc_metric_distributions.pdf",
            )
        ]
    return targets


def get_igv_track_files():
    """ATAC_norm.bw (per cluster) + thresholded .bedpe (per cluster/model),
    produced by scE2G's own (unexcluded) frag_to_norm_bigWig/
    write_sc_e2g_predictions_bedpe rules -- targeted directly by path, no new
    rule needed. Gated on scE2G_options.make_IGV_tracks exactly like scE2G's
    own `rule all` gates the same files, since we don't reuse that rule."""
    if not config.get("scE2G_options", {}).get("make_IGV_tracks", False):
        return []
    files = []
    for dataset, cluster in INCLUDED_CLUSTERS:
        files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, "ATAC_norm.bw"))
        for model in config["clusters"][dataset][cluster]["models"]:
            threshold = get_model_threshold(config["scE2G_dir"], model)
            files.append(os.path.join(
                RESULTS_DIR_BASE, dataset, cluster, model,
                f"scE2G_predictions_threshold{threshold}.bedpe",
            ))
    return files


rule aggregate_qc_stats:
    input:
        this_run_stats=lambda wildcards: per_cluster_stats_files(wildcards.dataset),
    params:
        dataset_dir=os.path.join(RESULTS_DIR_BASE, "{dataset}"),
        plots_dir=os.path.join(WDIR, "plots", "{dataset}"),
    output:
        all_stats=os.path.join(QC_PLOTS_DIR, "all_qc_stats.tsv"),
    conda:
        "../envs/filter_multiome_env.yaml"
    resources:
        mem_mb=8000,
    shell:
        """
        python {workflow.basedir}/scripts/aggregate_qc_stats.py {params.dataset_dir} {params.plots_dir} {output.all_stats}
        """


rule plot_qc_stats:
    input:
        all_stats=os.path.join(QC_PLOTS_DIR, "all_qc_stats.tsv"),
    output:
        ds_stats=os.path.join(QC_PLOTS_DIR, "dataset_metrics.pdf"),
        enh_stats=os.path.join(QC_PLOTS_DIR, "enhancer_metrics.pdf"),
        eg_stats=os.path.join(QC_PLOTS_DIR, "eg_metrics.pdf"),
        gene_stats=os.path.join(QC_PLOTS_DIR, "gene_metrics.pdf"),
        dist_size_stats=os.path.join(QC_PLOTS_DIR, "distance_and_size_metrics.pdf"),
        all_distributions=os.path.join(QC_PLOTS_DIR, "qc_metric_distributions.pdf"),
    conda:
        SCE2G_SC_E2G_ENV
    resources:
        mem_mb=lambda wildcards, input, attempt: determine_mem_mb(wildcards, input, attempt),
    script:
        "../scripts/plot_all_qc_stats_from_merged.R"
