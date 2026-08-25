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

2026-08-24: portal_cell_metadata() now reads common.smk's CELL_ANNOTATIONS -- the
driver-written, timestamp- and digest-checked snapshot of state.db's
cell_annotations (workflow/scripts/cell_annotation_snapshot.py) -- rather than
opening state.db itself. Nothing under workflow/rules/ connects to that database
any more: every Slurm worker re-parses this file, state.db lives on Lustre in WAL
mode, and WAL's shared-memory index is not supported across hosts.

This file still never contacts the network, at parse time or at execution time.
The cache is warmed once, before Snakemake starts, by
workflow/scripts/run_pipeline.py's warm stage (cell_metadata.fetch_if_stale +
derive_scopes). In pipeline_mode=local_only none of this is touched at all and
get_reformat_output_files() returns nothing.
"""

import sys

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from cell_annotations import annotation_lookup_key
from igvf_metadata.context import Context, IgvfConfig
from igvf_metadata.tables.prediction_tabular_files import build_alias as _ptf_build_alias

IGVF_CFG = IgvfConfig.from_dict(config.get("igvf", {}))


def portal_cell_metadata(dataset, cluster):
    """{"cell_annotation":..., "cl_id":..., "term_id":..., "term_name":...,
    "cell_qualifier":..., ...} for (dataset, cluster) -- see this module's
    docstring for which of these fields feeds which header field.

    Reads common.smk's CELL_ANNOTATIONS, the driver-written snapshot of
    state.db's cell_annotations (see scripts/cell_annotation_snapshot.py). No
    SQLite connection here any more: this function is called from `params`
    lambdas, which every Slurm worker evaluates after re-parsing the workflow,
    and multi-host WAL access to a Lustre-hosted DB is unsupported.

    Resolves via the cluster's cell_annotation_key override when set (ATAC-only
    variant clusters -- the real CellAnnotation lives under the base name, not
    the suffixed cluster key), the same lookup common.smk's
    REFORMAT_ELIGIBLE_CLUSTERS gate uses, so a cluster that passed that gate can
    never fail this lookup."""
    if PIPELINE_MODE == "local_only":
        # Unreachable via rule all (get_reformat_output_files returns [] in this
        # mode), so only an explicitly-named reformat target gets here. Say why
        # rather than looking like a cold cache.
        raise ValueError(
            f"{dataset}/{cluster}: reformatting requested while pipeline_mode=local_only, which "
            "deliberately makes no CellAnnotation available. Re-run in default mode via "
            "workflow/scripts/run_pipeline.py."
        )
    row = CELL_ANNOTATIONS.get(annotation_lookup_key(dataset, cluster, config["clusters"][dataset][cluster]))
    if row is None:
        raise ValueError(
            f"{dataset}/{cluster}: no Cell Annotation in the snapshot -- this cluster's portal "
            "primaries didn't resolve. See the driver's warm-stage status report for the reason."
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
    # REFORMAT_ELIGIBLE_CLUSTERS (common.smk), not UPLOAD_ELIGIBLE_CLUSTERS --
    # portal_cell_metadata() above raises if a (dataset, cluster) has no
    # Cell Annotation. Requesting reformat output for a quality-passing cluster
    # that simply hasn't been portal-annotated yet would crash `rule all`
    # outright; core generation (fragments, RNA matrices, predictions,
    # candidates, features) is never gated this way.
    #
    # Empty in local_only mode, by construction: REFORMAT_ELIGIBLE_CLUSTERS is
    # empty there. That is what makes `--omit-from reformat_predictions
    # reformat_predictions_thresholded reformat_lists` unnecessary -- and the
    # nargs='+' target-swallowing footgun that flag carries avoidable entirely.
    files = []
    for dataset, cluster in REFORMAT_ELIGIBLE_CLUSTERS:
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
