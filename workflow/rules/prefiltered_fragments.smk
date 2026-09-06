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
    """Make sure an `atac_frag_file` cluster's fragments are chromosome-grouped
    before scE2G sees them.

    NEEDED because scE2G does not do this and assumes it has been done. Its
    frag_to_norm_bigWig pipes fragments straight into `bedtools genomecov -bg`,
    which requires chromosome-GROUPED input, and only sorts the resulting
    bedGraph -- too late to repair overlapping intervals genomecov already
    emitted. ~18% of the CATlas files revisit chromosomes later in the file and
    every one of those fails with "Error - overlapping regions in bedGraph".
    See pipeline_paths.sorted_frag_file for the full diagnosis.

    INTERIM: scE2G's own frag_to_tagAlign sorts the identical input via
    bedSplitSort.sh, so the real fix is upstream. Delete this rule when that
    lands.

    ONLY SORTS WHAT NEEDS SORTING. A single streaming pass decides; if the file
    is already grouped and position-sorted (about 82% of them) it is HARDLINKED
    instead of copied, which costs no time and no space. That matters a lot:
    sorting everything unconditionally made this rule the dominant cost of the
    whole 258-cluster run, and the largest files (25-52 GiB compressed, hundreds
    of GB decompressed) blew through the profile's 6h default and had to be
    retried -- Bipolar_cells exhausted all 3 attempts and was lost.
    The hardlink is safe: it shares an inode with the download, and temp()
    removal unlinks only our link, never the original.

    runtime is declared HERE rather than inherited. The profile's `runtime="6h"`
    is fine for the pipeline's other rules but not for sorting a 52 GiB
    fragments file to Lustre.

    sort -T points at a DEDICATED directory, never the output directory. When a
    sort is killed mid-flight (walltime, preemption) GNU sort cannot clean up,
    and spill files left beside the outputs are indistinguishable from them --
    the first version of this rule left 102 orphaned spills, the largest 213 GB,
    accounting for most of 2.6 TB on $SCRATCH.
    """
    input:
        frag_file=cluster_prefiltered_frag_file,
    output:
        sorted_frags=temp(
            os.path.join(
                sorted_fragments_dir(config), "{dataset}", "{cluster}_fragments.sorted.tsv.gz",
            )
        ),
    params:
        # Sibling of the outputs, not the outputs' own directory -- see docstring.
        sort_tmp=os.path.join(sorted_fragments_dir(config), "sort_tmp"),
    threads: 8
    resources:
        mem_mb=determine_mem_mb,
        runtime=24 * 60,
    conda:
        "../envs/filter_multiome_env.yaml"
    shell:
        """
        set -euo pipefail
        LC_ALL=C
        mkdir -p {params.sort_tmp}

        # One streaming pass: exit 1 the moment a chromosome is revisited or a
        # position goes backwards within a run. No sort, no spill.
        if zcat {input.frag_file} | awk '
              {{ if ($1 != p) {{ if ($1 in s) exit 1; s[$1]=1; p=$1; q=-1 }}
                 if ($2+0 < q) exit 1; q=$2+0 }}
              END {{ exit 0 }}'; then
            echo "[sort_atac_fragments] {wildcards.cluster}: already grouped and sorted -- hardlinking"
            ln {input.frag_file} {output.sorted_frags} \
              || cp {input.frag_file} {output.sorted_frags}
        else
            echo "[sort_atac_fragments] {wildcards.cluster}: interleaved -- sorting"
            BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} \
                'BEGIN {{ result = mem_mb/threads/2; print int(result) "M" }}')
            zcat {input.frag_file} \
                | sort -k1,1 -k2,2n --parallel={threads} -S $BUFFER_SIZE -T {params.sort_tmp} \
                | bgzip -c > {output.sorted_frags}
        fi
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
