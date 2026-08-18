"""
Portal-facing bgzip+tabix packaging of scE2G's ABC Neighborhoods candidate
enhancer list. EnhancerList.bed (produced by scE2G's own
ENCODE_rE2G/ABC/workflow/scripts/neighborhoods.py, already coordinate-sorted)
is a byproduct of the predictions DAG, not a named scE2G output -- these
rules make it a first-class, indexed target of our own.

Outputs land at the CLUSTER level ({dataset}/{cluster}/EnhancerList.bed.gz[.tbi]),
one directory above Neighborhoods/, not inside it: scE2G's own neighborhoods
rule declares the whole Neighborhoods/ directory as a Snakemake directory()
output, and Snakemake raises ChildIOException for any other rule's output
nested inside another rule's directory() output (confirmed empirically
2026-08-15) -- there's no way to make this pipeline's own outputs literal
siblings of EnhancerList.bed within that directory.

Scoped to UPLOAD_ELIGIBLE_CLUSTERS (quality-passing, independent of
CellAnnotation availability -- same scope as candidates_features.smk, since
core generation is never gated on portal-metadata availability). Not
model-scoped: Neighborhoods/ is shared across every model for a cluster
(RESULTS_DIR/{cluster}/Neighborhoods/..., no {model} component in the path).
"""

ENHANCER_LIST_SRC_DIR = os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "Neighborhoods")
ENHANCER_LIST_OUT_DIR = os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}")


rule bgzip_enhancer_list:
    input:
        # EnhancerList.bed itself is an undeclared side-effect of scE2G's own
        # neighborhoods rule -- only EnhancerList.txt/GeneList.txt and the
        # whole Neighborhoods directory (via directory()) are declared
        # outputs, so Snakemake can only resolve this dependency through the
        # directory, not the individual .bed file (confirmed empirically:
        # requesting the .bed path directly raises MissingInputException).
        neighborhood_dir=ENHANCER_LIST_SRC_DIR,
    output:
        bedgz=os.path.join(ENHANCER_LIST_OUT_DIR, "EnhancerList.bed.gz"),
    conda:
        "../envs/filter_multiome_env.yaml"  # has bioconda::htslib (bgzip, tabix) already
    shell:
        "bgzip -c {input.neighborhood_dir}/EnhancerList.bed > {output.bedgz}"


rule tabix_enhancer_list:
    input:
        bedgz=os.path.join(ENHANCER_LIST_OUT_DIR, "EnhancerList.bed.gz"),
    output:
        tbi=os.path.join(ENHANCER_LIST_OUT_DIR, "EnhancerList.bed.gz.tbi"),
    conda:
        "../envs/filter_multiome_env.yaml"
    shell:
        # EnhancerList.bed is already coordinate-sorted -- no sort step needed.
        "tabix -p bed {input.bedgz}"


def get_enhancer_list_targets():
    files = []
    for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS:
        d = ENHANCER_LIST_OUT_DIR.format(dataset=dataset, cluster=cluster)
        files += [os.path.join(d, "EnhancerList.bed.gz"), os.path.join(d, "EnhancerList.bed.gz.tbi")]
    return files
