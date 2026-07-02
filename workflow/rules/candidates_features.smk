"""
Candidate E2G pairs + scE2G feature tables.

Forked from IGVF/workflow/rules/e2g_candidates_uni_proc.smk and
scE2G_features.smk -- both files define a rule confusingly named
`isolate_e2g_candidate_files` (a real naming collision, not a typo); both are
reimplemented here under distinct names.

Forking (not verbatim copying) was necessary because both original rules
hardcode their input to .../multiome_powerlaw_v3/scE2G_predictions.tsv.gz,
which breaks for ATAC-only clusters. Here the source model is
resolve_primary_model(cluster's models) instead: multiome_powerlaw_v3 when
both models ran, else scATAC_powerlaw_v3 -- exactly one candidates file and
one feature table per cluster, and the ATAC-only feature table's header
correctly reads "scATAC_powerlaw_v3" (see the forked
isolate_scE2G_feature_columns.R) rather than a hardcoded multiome string.

The candidates rule's wildcards are `dataset` + `cluster` (both genuine,
jointly-constrained) instead of the original's single `cell_type`, which
silently assumed cell_type == cluster == globally unique -- true for the old
McGinnis datasets, not guaranteed here where pseudobulk_annotation can
differ from cluster (split-cluster case) and cluster names can recur across
datasets. Candidates output intentionally stays FLAT (no {dataset}/{cluster}
nesting), matching the original -- only the wildcard set was wrong, not the
output layout; the dataset+cluster prefix in the filename is enough to keep
outputs across datasets from colliding.
"""

CANDIDATES_DIR = os.path.join(RESULTS_DIR, "candidate_e2g_pairs")
FEATURES_DIR = os.path.join(RESULTS_DIR, "scE2G_features")


def get_candidates_output_files():
    return [
        os.path.join(CANDIDATES_DIR, f"{dataset}_{cluster}_candidate_e2g_pairs.tsv.gz")
        for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS
    ]


def get_features_output_files():
    files = []
    for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS:
        model = resolve_primary_model(config["clusters"][dataset][cluster]["models"])
        files.append(os.path.join(FEATURES_DIR, dataset, cluster, f"{dataset}_{cluster}_scE2G_{model}_features.tsv.gz"))
    return files


def candidates_input_model_dir(wildcards):
    model = resolve_primary_model(config["clusters"][wildcards.dataset][wildcards.cluster]["models"])
    return os.path.join(RESULTS_DIR, wildcards.dataset, wildcards.cluster, model, "scE2G_predictions.tsv.gz")


rule isolate_candidate_pairs:
    input:
        candidates_input_model_dir,
    output:
        os.path.join(CANDIDATES_DIR, "{dataset}_{cluster}_candidate_e2g_pairs.tsv.gz"),
    params:
        summary=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["summary"],
    conda:
        "../envs/e2g_jamboree_env.yml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        Rscript {workflow.basedir}/scripts/isolate_candidate_e2g_pairs.R \
        -i {input} \
        -o {output} \
        -s "{params.summary}"
        """


rule isolate_feature_table:
    input:
        candidates_input_model_dir,
    output:
        os.path.join(FEATURES_DIR, "{dataset}", "{cluster}", "{dataset}_{cluster}_scE2G_{model}_features.tsv.gz"),
    params:
        model="{model}",
        cell_type=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["cell_type"],
        term_id=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["ontology_id"],
        summary=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["summary"],
    conda:
        "../envs/e2g_jamboree_env.yml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        Rscript {workflow.basedir}/scripts/isolate_scE2G_feature_columns.R \
        -i {input} \
        -o {output} \
        -c "{params.cell_type}" \
        -d "{params.term_id}" \
        -s "{params.summary}" \
        -m {params.model}
        """
