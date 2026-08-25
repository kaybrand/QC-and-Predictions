"""
ATAC/RNA filtering for the (dataset, cluster) pairs this run is processing
(INCLUDED_CLUSTERS, from common.smk). Mirrors the legacy top-level
Snakefile's atac_fragment_file/rna_count_matrix rules, but per (dataset,
cluster) pair -- both are genuine wildcards, jointly constrained -- and
ATAC-only aware: the RNA rule is simply never instantiated for a cluster
whose models == [scATAC_powerlaw_v3].

Outputs land in {output_dir}/multiome_data/{dataset}/{cluster}/, the same
staging convention the legacy Snakefile and the qc-filter-pseudobulks skill
used for data_dir -- scE2G reads these paths via the cell_clusters table
written by write_scE2G_config.py, not from this pipeline's own results_dir.

OUT_DIR_BASE matches common.smk's own multiome_data_dir(dataset) minus the
{dataset} join (both read OUTPUT_DIR) -- 2026-08-15: repointed from
config["data_dir"] to the new consolidated OUTPUT_DIR (config["output_dir"],
default ./results), since data_dir is now read-only (QC-guide plots/
datatables inputs only). Before that, corrected 2026-08-03 from WDIR-relative
(this pipeline's own code checkout), which only worked by coincidence when
code and data lived in the same place.
"""

OUT_DIR_BASE = os.path.join(OUTPUT_DIR, "multiome_data")

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
        # Also declared as its own directory() output: scE2G's own
        # generate_atac_matrix rule (the Kendall-feature branch) takes this
        # directory itself as a literal input path (via the cell_clusters
        # table's rna_matrix_file column), not the three files above. Without
        # this, Snakemake can't find a producing rule for that literal path
        # on a from-scratch dataset and raises MissingInputException.
        rna_dir=directory(os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}")),
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


rule package_rna_count_matrix:
    """IGVF Portal packaging for the RNA count matrix -- Filtered Matrix
    Files' own file format spec (QC_pseudobulks/FILE_SPEC_RNA_COUNT_MATRIX.txt)
    requires a tar.gz of DEcompressed matrix.mtx/barcodes.tsv/features.tsv,
    flat (no subdirectory nesting inside the archive) -- distinct from the
    gzipped-individually directory rule rna_count_matrix produces above,
    which is untouched here and remains scE2G's own input (referenced
    directly by the cell_clusters table write_scE2G_config.py writes;
    scE2G needs the untarred, per-file-gzipped form, never this tarball).

    No Snakemake rule downstream of this one consumes its output -- unlike
    atac_fragment_file/rna_count_matrix, which get pulled into the DAG
    transitively via scE2G's own inputs, this tarball exists only for the
    IGVF metadata uploader (igvf_metadata.tables.filtered_rna_count_matrix)
    to find on disk, so it needs its own explicit target
    (get_rna_matrix_packages(), wired into rule all in the top-level
    Snakefile).
    """
    input:
        matrix=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}", "matrix.mtx.gz"),
        barcodes=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}", "barcodes.tsv.gz"),
        features=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}", "features.tsv.gz"),
    output:
        tarball=os.path.join(OUT_DIR_BASE, "{dataset}", "{cluster}", "rna_count_matrix_{dataset}_{cluster}.tar.gz"),
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        tmp={resources.tmpdir}/rna_count_matrix_{wildcards.dataset}_{wildcards.cluster}
        mkdir -p "$tmp"
        zcat {input.matrix}    > "$tmp/matrix.mtx"
        zcat {input.barcodes}  > "$tmp/barcodes.tsv"
        zcat {input.features}  > "$tmp/features.tsv"
        tar -czf {output.tarball} -C "$tmp" matrix.mtx barcodes.tsv features.tsv
        rm -rf "$tmp"
        """


def get_rna_matrix_packages():
    """Tarball targets for every upload-eligible cluster with RNA data --
    excluded clusters (unless predictions_on_everything_but_do_not_upload) never get
    distribution-format outputs, same reasoning as
    reformat.smk's get_reformat_output_files()."""
    return [
        os.path.join(OUT_DIR_BASE, dataset, cluster, f"rna_count_matrix_{dataset}_{cluster}.tar.gz")
        for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS
        if config["clusters"][dataset][cluster]["models"] != ["scATAC_powerlaw_v3"]
    ]
