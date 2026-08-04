"""
ATAC/RNA filtering for the (dataset, cluster) pairs this run is processing
(INCLUDED_CLUSTERS, from common.smk). Mirrors the legacy top-level
Snakefile's atac_fragment_file/rna_count_matrix rules, but per (dataset,
cluster) pair -- both are genuine wildcards, jointly constrained -- and
ATAC-only aware: the RNA rule is simply never instantiated for a cluster
whose models == [scATAC_powerlaw_v3].

Outputs land in {data_dir}/multiome_data/{dataset}/{cluster}/, the same
staging convention the legacy Snakefile and the qc-filter-pseudobulks skill
use -- scE2G reads these paths via the cell_clusters table written by
write_scE2G_config.py, not from this pipeline's own results_dir.

OUT_DIR_BASE matches common.smk's own multiome_data_dir(dataset) minus the
{dataset} join (both read config["data_dir"]) -- corrected 2026-08-03, this
used to be WDIR-relative (this pipeline's own code checkout), which only
worked by coincidence when code and data lived in the same place.
"""

OUT_DIR_BASE = os.path.join(config["data_dir"], "multiome_data")

MULTIOME_CLUSTERS = [
    (d, c) for d, c in INCLUDED_CLUSTERS if config["clusters"][d][c]["models"] != ["scATAC_powerlaw_v3"]
]


def cluster_pseudobulk_annotation(wildcards):
    return config["clusters"][wildcards.dataset][wildcards.cluster]["pseudobulk_annotation"]


def cluster_qc_guide(wildcards):
    return config["clusters"][wildcards.dataset][wildcards.cluster]["qc_guide"]


rule atac_fragment_file:
    input:
        barcodes=cluster_qc_guide,
    output:
        filtered_atac=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "atac_fragments_{dataset}_{cluster}.tsv.gz"),
        index=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "atac_fragments_{dataset}_{cluster}.tsv.gz.tbi"),
    params:
        pseudobulks=lambda wildcards: pseudobulks_dir(wildcards.dataset),
        cell_type=cluster_pseudobulk_annotation,
        chrom_sizes=os.path.join(WDIR, "reference", "IGVF.DACC.GRCh38.chrom.sizes.tsv"),
    conda:
        "../envs/filter_multiome_env.yaml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        python {workflow.basedir}/scripts/filter_atac_fragments.py \
            --qc-guide    {input.barcodes} \
            --pseudobulks {params.pseudobulks} \
            --cell-type   {params.cell_type} \
            --chrom-sizes {params.chrom_sizes} \
            --out         {output.filtered_atac}
        """


rule rna_count_matrix:
    input:
        barcodes=cluster_qc_guide,
    output:
        matrix=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}", "matrix.mtx.gz"),
        barcodes_out=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}", "barcodes.tsv.gz"),
        features_out=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}", "features.tsv.gz"),
    params:
        pseudobulks=lambda wildcards: pseudobulks_dir(wildcards.dataset),
        cell_type=cluster_pseudobulk_annotation,
        rna_out=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}.mtx"),
        gtf=config["gtf"],
    log:
        os.path.join(OUT_DIR_BASE, "{dataset}", "logs", "{cluster}_gtf_mapping.txt"),
    conda:
        "../envs/filter_multiome_env.yaml"
    resources:
        mem_mb=determine_mem_mb,
    # No per-rule wildcard_constraints here: `cluster` alone can't express "only
    # multiome (dataset, cluster) pairs" since cluster names aren't unique across
    # datasets. Safety instead comes from never requesting this rule's output for
    # an ATAC-only pair anywhere in this pipeline's own target-building code
    # (get_reformat_output_files, write_cell_clusters_table, get_filtered_data_files).
    shell:
        """
        python {workflow.basedir}/scripts/filter_rna_counts.py \
            --qc-guide    {input.barcodes} \
            --pseudobulks {params.pseudobulks} \
            --cell-type   {params.cell_type} \
            --out         {params.rna_out} \
            --gtf         {params.gtf} \
            --standard-chromosomes-only \
            --log         {log}
        """
