"""
Reformat scE2G predictions/lists into the standard sharing format.

Adapted from IGVF/workflow/rules/reformat.smk: the three target rules below
(reformat_predictions, reformat_predictions_thresholded, reformat_lists) are
unchanged in logic -- model is already a wildcard and
update_scE2G_pred_formats.R already branches on `grepl("ATAC", opt$method)`,
so no forking was needed here (unlike candidates_features.smk). The
`update_*` Synapse-upload rules from the original file are intentionally not
carried over -- Synapse upload for this pipeline is handled by
synapse_manifest.smk instead. Thresholds come from get_model_threshold()
(reads models/{model}/score_threshold_*) instead of a static config['models']
dict.

RESULTS_DIR_BASE here is the dataset-AGNOSTIC parent
(scE2G_dir/results/uniformly_processed) -- {dataset} is a genuine,
jointly-constrained wildcard alongside {cluster}, matching scE2G's own
per-dataset results_dir (RESULTS_DIRS[dataset] in common.smk) one level
down. This also matches the original IGVF/workflow/rules/reformat.smk's own
convention (its results_dir spanned all datasets too).

2026-07-23 header-field update (per the IGVF Data Portal manager feedback
pass): SampleTermName/SampleTermID/CellAnnotation (formerly SampleSummaryShort
-- renamed to match the portal's own field name) now all come from the IGVF
Portal itself, via igvf_metadata.state's cell_annotations cache
(portal_cell_metadata() below), not CLUSTER_METADATA/lab_annotations_with_cl.tsv
-- that local table was a stand-in until the portal-side cache existed.
SampleTermName is cell_type's portal "term_name" (e.g. "macrophage");
SampleTermID is cell_type's portal "term_id", the CURIE form (e.g.
"CL:0000235") -- NOT cl_id, which is the "@id"-derived, underscore form used
elsewhere (e.g. Principal Pseudobulk Set's own sample_terms reference).
CellAnnotation is still cell_annotation, unchanged. A new `--portal-link`
carries the file's own IGVF alias (computed via igvf_metadata's own
build_alias, the single source of truth for that formula -- not reimplemented
here), which the R script spells out into a full https://data.igvf.org/... URL.

STUB, flagged for whoever wires up the "pipeline-integrated Snakemake rule"
manage_igvf_metadata.py's own docstring already anticipates: portal_cell_metadata()
below only READS an already-populated state.db (plain read-only sqlite3, no
network/igvf_utils needed at Snakemake parse time) -- it does not trigger a
live portal refresh itself. Assumes some other invocation (manage_igvf_metadata.py,
run standalone, or that future pipeline-integrated rule) has already refreshed
the cache -- via cell_metadata.refresh_if_stale, gated on cache age > 24h AND
(reformatting requested OR Principal Pseudobulk Set construction) -- for
these clusters. Raises clearly, per cluster, if it hasn't.
"""

import sys

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from igvf_metadata import state as igvf_state
from igvf_metadata.context import Context, IgvfConfig
from igvf_metadata.tables.prediction_tabular_files import build_alias as _ptf_build_alias

IGVF_CFG = IgvfConfig.from_dict(config.get("igvf", {}))
_STATE_DB_PATH = config.get("igvf", {}).get("state_db_path")
_STATE_CONN = igvf_state.connect(_STATE_DB_PATH) if _STATE_DB_PATH else None


def portal_cell_metadata(dataset, cluster):
    """{"cell_annotation":..., "cl_id":..., "term_id":..., "term_name":...,
    "cell_qualifier":..., ...} for (dataset, cluster) from the IGVF Portal
    cell-metadata cache -- see this module's docstring for what this does
    (and doesn't) do, and for which of these fields feeds which header
    field."""
    if _STATE_CONN is None:
        raise ValueError(
            "igvf.state_db_path not set in the pipeline config -- required for "
            "SampleTermID/SampleSummaryShort, sourced from the IGVF Portal "
            "cell-metadata cache (see workflow/scripts/igvf_metadata/cell_metadata.py)"
        )
    row = igvf_state.get_cell_annotation(_STATE_CONN, dataset, cluster)
    if row is None:
        raise ValueError(
            f"{dataset}/{cluster}: no cached Cell Annotation metadata in state.db -- "
            "refresh the cache for this cluster first (e.g. via manage_igvf_metadata.py)"
        )
    return row


def portal_link_alias(dataset, cluster, model, variant):
    """The file's own IGVF alias, via Prediction Tabular Files' own
    build_alias -- single source of truth for that formula."""
    ctx = Context(
        dataset, cluster, model, config["clusters"][dataset][cluster], IGVF_CFG,
        config["scE2G_dir"], config["data_dir"], OUTPUT_DIR,
    )
    return _ptf_build_alias(ctx, variant)


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
        name=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_name"],
        term_id=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_id"],
        summary=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["cell_annotation"],
        portal_link=lambda wildcards: portal_link_alias(wildcards.dataset, wildcards.cluster, wildcards.model, "full"),
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
        -l "{params.portal_link}"
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
        name=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_name"],
        term_id=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_id"],
        summary=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["cell_annotation"],
        portal_link=lambda wildcards: portal_link_alias(
            wildcards.dataset, wildcards.cluster, wildcards.model, "thresholded"
        ),
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
        -l "{params.portal_link}" \
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
        name=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_name"],
        term_id=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_id"],
        summary=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["cell_annotation"],
        # reformat_lists only ever runs against multiome_powerlaw_v3 (see input: above) --
        # "elements"/"genes" is Prediction Tabular Files' own variant naming (plural),
        # not the {meta} wildcard's singular "element"/"gene".
        portal_link=lambda wildcards: portal_link_alias(
            wildcards.dataset, wildcards.cluster, "multiome_powerlaw_v3",
            "elements" if wildcards.meta == "element" else "genes",
        ),
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
        -l "{params.portal_link}"
        """
