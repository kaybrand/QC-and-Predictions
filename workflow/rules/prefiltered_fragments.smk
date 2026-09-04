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

from pipeline_paths import sorted_fragments_dir, sorted_frag_file  # noqa: E402

PRE_FILTERED_CLUSTERS = [
    (d, c) for d, c in INCLUDED_CLUSTERS if "atac_frag_file" in config["clusters"][d][c]
]


def cluster_prefiltered_frag_file(wildcards):
    return config["clusters"][wildcards.dataset][wildcards.cluster]["atac_frag_file"]


rule sort_atac_fragments:
    """Coordinate-sort an `atac_frag_file` cluster's fragments before scE2G sees
    them.

    NEEDED because scE2G does not do this and assumes it has been done. Its
    frag_to_norm_bigWig pipes fragments straight into `bedtools genomecov -bg`,
    which requires chromosome-GROUPED input, and only sorts the resulting
    bedGraph -- too late to repair overlapping intervals genomecov already
    emitted. ~18% of the CATlas files revisit chromosomes later in the file
    (positions ARE sorted within each run), and every one of those fails with
    "Error - overlapping regions in bedGraph". `fragments_preprocessed: False`
    does not help: that branch only filters to chromosomes in chr_sizes.
    See pipeline_paths.sorted_frag_file for the full diagnosis.

    INTERIM. scE2G's own frag_to_tagAlign sorts the identical input via
    bedSplitSort.sh, so the real fix is upstream in frag_to_norm_bigWig /
    frag_to_bigWig; this rule can be deleted once that lands.

    `sort -k1,1 -k2,2n` under LC_ALL=C is exactly what scE2G's config comment
    prescribes ("sorted with sort -k1,1 2n,n"), and matches the lexicographic
    order its own post-genomecov sort already produces -- so the bedGraph that
    reaches bedGraphToBigWig is ordered the same way it is today for the files
    that happen to work.

    temp(): these are full-size copies of the fragment tree (~776 GB for the 257
    human clusters), so they are deleted once nothing needs them. The tradeoff
    is that re-running anything downstream re-sorts; drop the temp() wrapper if
    that churn ever costs more than the disk.

    Written to $SCRATCH (via download_root's sibling, see
    pipeline_paths.sorted_fragments_dir), never $OAK -- job data I/O, and $OAK
    is near quota."""
    input:
        frag_file=cluster_prefiltered_frag_file,
    output:
        sorted_frags=temp(
            os.path.join(
                sorted_fragments_dir(config), "{dataset}", "{cluster}_fragments.sorted.tsv.gz",
            )
        ),
    params:
        temp_dir=lambda wildcards, output: os.path.dirname(output.sorted_frags),
    threads: 4
    resources:
        mem_mb=determine_mem_mb,
    conda:
        "../envs/filter_multiome_env.yaml"
    shell:
        """
        LC_ALL=C
        BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} \
            'BEGIN {{ result = mem_mb/threads/2; print int(result) "M" }}')
        zcat {input.frag_file} \
            | sort -k1,1 -k2,2n --parallel={threads} -S $BUFFER_SIZE -T {params.temp_dir} \
            | bgzip -c > {output.sorted_frags}
        """


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
