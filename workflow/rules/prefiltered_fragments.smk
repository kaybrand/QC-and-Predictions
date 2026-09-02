"""
Stats for clusters whose ATAC fragments arrive already filtered/QC'd (no QC
guide barcode list, no per-cell QC table -- e.g. catlas's WashU pseudobulks,
pulled pre-filtered straight from the IGVF Data Portal). These clusters are
marked by an `atac_frag_file` key in their cluster config (see
write_scE2G_config.py::write_cell_clusters_table and
resolve_exclusions.py::compute_cluster_stats for the other two places this
same marker is checked) instead of `pseudobulk_annotation`/`qc_guide`.

filter_pseudobulks.smk's atac_fragment_file/rna_count_matrix rules are simply
never instantiated for these clusters -- nothing in the DAG ever requests
their filtered-output paths -- so no changes were needed there. This file
supplies the one thing that flow WOULD have produced as a side effect that
these clusters still need: filtered_cell_subsample_metrics.tsv, so
resolve_exclusions.py/aggregate_qc_stats.py can compute real
cell_count/fragments_total for them exactly as they already do for every
other dataset, with no changes to either script's actual file-parsing logic.
"""

PRE_FILTERED_CLUSTERS = [
    (d, c) for d, c in INCLUDED_CLUSTERS if "atac_frag_file" in config["clusters"][d][c]
]


def cluster_prefiltered_frag_file(wildcards):
    return config["clusters"][wildcards.dataset][wildcards.cluster]["atac_frag_file"]


rule prefiltered_cell_subsample_metrics:
    input:
        frag_file=cluster_prefiltered_frag_file,
    output:
        metrics=os.path.join(config["data_dir"], "plots", "{dataset}", "{cluster}", "filtered_cell_subsample_metrics.tsv"),
    conda:
        "../envs/filter_multiome_env.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/compute_prefiltered_cell_metrics.py \
            --frag-file {input.frag_file} \
            --subsample-name {wildcards.cluster} \
            --out {output.metrics}
        """
