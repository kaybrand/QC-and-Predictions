"""
Reformat scE2G predictions/lists into the standard sharing format.

Adapted from IGVF/workflow/rules/reformat.smk: the three target rules below
(reformat_predictions, reformat_predictions_thresholded, reformat_lists) are
unchanged in logic -- model is already a wildcard and
update_scE2G_pred_formats.R already branches on `grepl("ATAC", opt$method)`,
so no forking was needed here (unlike candidates_features.smk). The
`update_*` Synapse-upload rules from the original file are intentionally not
carried over -- Synapse upload for this pipeline is handled by
synapse_manifest.smk instead. `id`/`cell_type`/`summary` now come from
CLUSTER_METADATA[(dataset, cluster)] (common.smk, joined from
lab_annotations_with_cl.tsv) instead of the pipeline config directly, and
thresholds come from get_model_threshold() (reads
models/{model}/score_threshold_*) instead of a static config['models'] dict.

RESULTS_DIR_BASE here is the dataset-AGNOSTIC parent
(scE2G_dir/results/uniformly_processed) -- {dataset} is a genuine,
jointly-constrained wildcard alongside {cluster}, matching scE2G's own
per-dataset results_dir (RESULTS_DIRS[dataset] in common.smk) one level
down. This also matches the original IGVF/workflow/rules/reformat.smk's own
convention (its results_dir spanned all datasets too).
"""

def get_reformat_output_files():
    files = []
    for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS:
        models = config["clusters"][dataset][cluster]["models"]
        for model in models:
            threshold = get_model_threshold(config["scE2G_dir"], model)
            files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_scE2G_{model}.e2g.tsv.gz"))
            files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_scE2G_{model}_threshold{threshold}.e2g.tsv.gz"))
        if "multiome_powerlaw_v3" in models:
            for meta in ["element", "gene"]:
                files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_scE2G_multiome_v3_{meta}_list.tsv.gz"))
    return files


rule reformat_predictions:
    input:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{model}", "scE2G_predictions.tsv.gz"),
    output:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_scE2G_{model}.e2g.tsv.gz"),
    params:
        model="{model}",
        version=config["scE2G_version"],
        name=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["cell_type"],
        term_id=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["ontology_id"],
        summary=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["summary"],
    conda:
        "../envs/e2g_jamboree_env.yml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        Rscript {workflow.basedir}/scripts/update_scE2G_pred_formats.R \
        -i {input} \
        -o {output} \
        -c "{params.name}" \
        -d "{params.term_id}" \
        -s "{params.summary}" \
        -m scE2G_{params.model} \
        -v {params.version}
        """


rule reformat_predictions_thresholded:
    input:
        lambda wildcards: os.path.join(
            RESULTS_DIR_BASE, wildcards.dataset, wildcards.cluster, wildcards.model,
            f"scE2G_predictions_threshold{get_model_threshold(config['scE2G_dir'], wildcards.model)}.tsv.gz",
        ),
    output:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_scE2G_{model}_threshold{threshold}.e2g.tsv.gz"),
    params:
        model="{model}",
        version=config["scE2G_version"],
        name=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["cell_type"],
        term_id=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["ontology_id"],
        summary=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["summary"],
        threshold="{threshold}",
    conda:
        "../envs/e2g_jamboree_env.yml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        Rscript {workflow.basedir}/scripts/update_scE2G_pred_formats.R \
        -i {input} \
        -o {output} \
        -c "{params.name}" \
        -d "{params.term_id}" \
        -s "{params.summary}" \
        -m scE2G_{params.model} \
        -v {params.version} \
        --threshold 'Score >= {params.threshold}'
        """


rule reformat_lists:
    input:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "multiome_powerlaw_v3", "scE2G_{meta}_list.tsv.gz"),
    output:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_scE2G_multiome_v3_{meta}_list.tsv.gz"),
    params:
        model="multiome_powerlaw_v3",
        version=config["scE2G_version"],
        name=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["cell_type"],
        term_id=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["ontology_id"],
        summary=lambda wildcards: CLUSTER_METADATA[(wildcards.dataset, wildcards.cluster)]["summary"],
    conda:
        "../envs/e2g_jamboree_env.yml"
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        Rscript {workflow.basedir}/scripts/update_scE2G_pred_formats.R \
        -i {input} \
        -o {output} \
        -c "{params.name}" \
        -d "{params.term_id}" \
        -s "{params.summary}" \
        -m scE2G_{params.model} \
        -v {params.version}
        """
