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

from pipeline_paths import (  # noqa: E402
    sorted_fragments_dir,
    sorted_frag_file,
    scE2G_chr_sizes,
    bed_split_sort_script,
)

PRE_FILTERED_CLUSTERS = [
    (d, c) for d, c in INCLUDED_CLUSTERS if "atac_frag_file" in config["clusters"][d][c]
]


def cluster_prefiltered_frag_file(wildcards):
    return config["clusters"][wildcards.dataset][wildcards.cluster]["atac_frag_file"]


rule sort_atac_fragments:
    """Put an `atac_frag_file` cluster's fragments into the chromosome order
    scE2G expects, before scE2G sees them.

    NEEDED because scE2G does not do this and assumes it has been done. Its
    frag_to_norm_bigWig pipes fragments straight into `bedtools genomecov -bg`,
    which requires each chromosome's records to be CONTIGUOUS, and only sorts
    the resulting bedGraph -- too late to repair overlapping intervals genomecov
    already emitted. ~18% of the CATlas files revisit chromosomes later in the
    file and every one of those failed with "Error - overlapping regions in
    bedGraph". See pipeline_paths.sorted_frag_file for the full diagnosis.

    "CORRECTLY SORTED" MEANS chr_sizes ORDER, not merely grouped. An earlier
    version of this rule only checked grouping and then GNU-sorted, which gives
    LEXICOGRAPHIC order (chr1 chr10 chr11 ... chr2) -- while chr_sizes is
    numeric (chr1 chr2 chr3 ...). That left the rule with two different output
    contracts depending on which branch it took, and lexicographic order is the
    exact thing scE2G 299ca9e / PR #94 had to fix for frag_to_tagAlign
    ("`sort -k 1,1V` placed chrM before chrX ... downstream bedtools commands
    expect chromosome order matching the user's chr_sizes file").

    So the check now requires all three, in one streaming pass:
      1. every chromosome's records contiguous (never revisited),
      2. positions non-decreasing within each chromosome's run,
      3. chromosomes appearing in chr_sizes ORDER, and none absent from it.
    Anything else is handed to scE2G's own bedSplitSort.sh, which splits by
    chromosome, sorts each by position and concatenates in the sizes file's
    order -- the same guarantee frag_to_tagAlign relies on, rather than a second
    hand-rolled sort that could disagree with it. It also drops chromosomes
    absent from the sizes file (warning about each), which is consistent with
    the chromosome filtering process_fragment_file applies downstream anyway.

    EXPECT ALMOST EVERY FILE TO BE SORTED. The CATlas files that are already
    grouped are grouped LEXICOGRAPHICALLY, which fails condition 3, so the
    fast path will rarely fire. It is kept because it is correct and free, not
    because it is expected to save work.

    bedSplitSort is memory-efficient by construction (per-chromosome sorts, not
    one monolithic sort), which is why it is the right tool for the 25-45 GiB
    files that timed out under the previous version: its own commit notes
    `bedtools sort -g` OOM-killed on a 797M-line file.

    INTERIM: the real fix is upstream, in scE2G's frag_to_norm_bigWig /
    frag_to_bigWig. Delete this rule when that PR lands. See
    .baton/2026-09-04-sce2g-sort-pr-baton.md.

    runtime is declared HERE rather than inherited; the profile's 6h is not
    enough for the largest files. sort_tmp is a DEDICATED directory, never the
    output directory -- a killed sort cannot clean up after itself, and spill
    files left beside the outputs are indistinguishable from them (the first
    version left 102 orphans, largest 213 GB).
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
        chr_sizes=scE2G_chr_sizes(config),
        bed_split_sort=bed_split_sort_script(config),
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

        # One streaming pass over the fragments, with chr_sizes read first to
        # establish the canonical chromosome ranks.
        #
        # NOTE the `bad` flag rather than `exit 1` directly. In awk, `exit N` in
        # a main block STILL RUNS END, so an `exit 0` there overwrites the status
        # and every file reports "sorted" -- which would pass interleaved files
        # through untouched and silently reinstate the bigWig failures. Caught by
        # validating against files of known order; do not "simplify" this back.
        if awk '
              NR == FNR {{ rank[$1] = FNR; next }}
              {{ if ($1 != prev) {{
                     if ($1 in seen) {{ bad=1; exit }}                 # revisited
                     if (!($1 in rank)) {{ bad=1; exit }}              # not in chr_sizes
                     if (rank[$1] < lastrank) {{ bad=1; exit }}        # wrong order
                     seen[$1]=1; prev=$1; lastrank=rank[$1]; pos=-1
                 }}
                 if ($2+0 < pos) {{ bad=1; exit }}                     # position regressed
                 pos=$2+0 }}
              END {{ exit (bad ? 1 : 0) }}' \
             {params.chr_sizes} <(zcat {input.frag_file}); then
            echo "[sort_atac_fragments] {wildcards.cluster}: already in chr_sizes order -- hardlinking"
            ln {input.frag_file} {output.sorted_frags} \
              || cp {input.frag_file} {output.sorted_frags}
        else
            echo "[sort_atac_fragments] {wildcards.cluster}: not in chr_sizes order -- bedSplitSort"
            BUFFER_SIZE=$(awk -v mem_mb={resources.mem_mb} -v threads={threads} \
                'BEGIN {{ result = mem_mb/threads/2; print int(result) "M" }}')
            SPLIT_DIR=$(mktemp -d -p {params.sort_tmp} bedSplitSort.XXXXXX)
            zcat {input.frag_file} \
                | bash {params.bed_split_sort} \
                    -i /dev/stdin \
                    -g {params.chr_sizes} \
                    -t "$SPLIT_DIR" \
                    -p {threads} \
                    -s "$BUFFER_SIZE" \
                    -T {params.sort_tmp} \
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
