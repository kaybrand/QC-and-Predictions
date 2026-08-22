"""
Reformat scE2G predictions/lists into the standard sharing format.

Adapted from IGVF/workflow/rules/reformat.smk: reformat_predictions,
reformat_predictions_thresholded, and reformat_gene_list (originally one
{meta}-wildcarded reformat_lists rule shared with elements, split 2026-08-19
once elements stopped being an update_scE2G_pred_formats.R output -- see
reformat_element_list below) are unchanged in logic -- model is already a
wildcard and update_scE2G_pred_formats.R already branches on
`grepl("ATAC", opt$method)`, so no forking was needed here (unlike
candidates_features.smk). The `update_*` Synapse-upload rules from the
original file are intentionally not carried over -- Synapse upload for this
pipeline is handled by synapse_manifest.smk instead. Thresholds come from
get_model_threshold() (reads models/{model}/score_threshold_*) instead of a
static config['models'] dict.

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
from cell_annotations import annotation_lookup_key
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
    field. Resolves via the cluster's cell_annotation_key override when set
    (ATAC-only variant clusters -- the real CellAnnotation is cached under
    the base name, not the suffixed cluster key), same as common.smk's
    REFORMAT_ELIGIBLE_CLUSTERS gate, so a cluster that passed that gate can
    never fail this lookup."""
    if _STATE_CONN is None:
        raise ValueError(
            "igvf.state_db_path not set in the pipeline config -- required for "
            "SampleTermID/SampleSummaryShort, sourced from the IGVF Portal "
            "cell-metadata cache (see workflow/scripts/igvf_metadata/cell_metadata.py)"
        )
    lookup_dataset, lookup_cluster = annotation_lookup_key(dataset, cluster, config["clusters"][dataset][cluster])
    row = igvf_state.get_cell_annotation(_STATE_CONN, lookup_dataset, lookup_cluster)
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
    # REFORMAT_ELIGIBLE_CLUSTERS (common.smk), not UPLOAD_ELIGIBLE_CLUSTERS --
    # portal_cell_metadata() below raises if a (dataset, cluster) has no
    # cached Cell Annotation. Requesting reformat output for a
    # quality-passing cluster that simply hasn't been portal-annotated yet
    # would crash `rule all` outright; core generation (fragments, RNA
    # matrices, predictions, candidates, features) is never gated this way.
    files = []
    for dataset, cluster in REFORMAT_ELIGIBLE_CLUSTERS:
        models = config["clusters"][dataset][cluster]["models"]
        for model in models:
            threshold = get_model_threshold(config["scE2G_dir"], model)
            files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_scE2G_{model}.e2g.tsv.gz"))
            files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_scE2G_{model}_threshold{threshold}.e2g.tsv.gz"))
        # element_list sources from Neighborhoods/EnhancerList.bed (see
        # reformat_element_list below), which every reformat-eligible cluster
        # has regardless of model -- not gated on "multiome_powerlaw_v3 in
        # models" like gene_list below, which needs Multiome's own RNA-derived
        # scE2G_gene_list.tsv.gz.
        files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_element_list.bed.gz"))
        files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_element_list.bed.gz.tbi"))
        if "multiome_powerlaw_v3" in models:
            files.append(os.path.join(RESULTS_DIR_BASE, dataset, cluster, f"{dataset}_{cluster}_scE2G_multiome_v3_gene_list.tsv.gz"))
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


rule reformat_gene_list:
    input:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "multiome_powerlaw_v3", "scE2G_gene_list.tsv.gz"),
    output:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_scE2G_multiome_v3_gene_list.tsv.gz"),
    params:
        model="multiome_powerlaw_v3",
        version=config["scE2G_version"],
        name=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_name"],
        term_id=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_id"],
        summary=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["cell_annotation"],
        # "genes" is Prediction Tabular Files' own variant naming (plural).
        portal_link=lambda wildcards: portal_link_alias(
            wildcards.dataset, wildcards.cluster, "multiome_powerlaw_v3", "genes",
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


rule reformat_element_list:
    """Sources from Neighborhoods/EnhancerList.bed, NOT the multiome_powerlaw_v3
    predictions tree -- confirmed 2026-08-19 on the IGVF Data Portal, elements
    are represented as an INPUT to the prediction file (the candidate universe
    considered before scoring), not something derived from it, and
    EnhancerList.bed IS that pre-scoring candidate universe: cluster-level,
    model-agnostic (Neighborhoods/ is shared across every model for a
    cluster), already coordinate-sorted, and already in
    "chr\\tstart\\tend\\t{class}|{chr}:{start}-{end}" form -- byte-for-byte the
    same geometry as multiome_powerlaw_v3/scE2G_element_list.tsv.gz's own
    chr/start/end/class/name columns (confirmed same row count, same order),
    just without the model-specific quantitative annotation (RPM/RPKM/
    activity_base/etc) that file also carries and that this output never
    wanted anyway. Hence no data.table/dplyr transform is needed here (unlike
    reformat_gene_list) -- just the same portal-metadata header block
    update_scE2G_pred_formats.R assembles elsewhere, prepended directly, plus
    the "#"-commented BED4 column-name line tabix/IGV/UCSC all skip. Depends
    on the whole Neighborhoods directory rather than EnhancerList.bed
    directly -- that file is an undeclared side-effect of scE2G's own
    neighborhoods rule, so Snakemake can only resolve the dependency through
    directory() (confirmed empirically: requesting the .bed path directly
    raises MissingInputException).

    Deliberately NOT named with a model or "multiome_v3" tag --
    {dataset}_{cluster}_element_list, matching what it actually is.
    enhancer_list_packaging.smk's own EnhancerList.bed.gz (unindexed,
    header-less bgzip of the same file) was removed 2026-08-21 once this
    rule made it fully redundant -- nothing outside that file's own rule
    ever consumed it."""
    input:
        neighborhood_dir=os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "Neighborhoods"),
    output:
        os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_element_list.bed"),
    params:
        version=config["scE2G_version"],
        name=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_name"],
        term_id=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["term_id"],
        summary=lambda wildcards: portal_cell_metadata(wildcards.dataset, wildcards.cluster)["cell_annotation"],
        # "elements_bed" is Prediction Tabular Files' own variant naming
        # (renamed from "elements" 2026-08-21 once file_format became "bed").
        # Hardcoded to the Multiome family alias regardless of which models
        # this cluster actually ran, same as before this rule's source
        # changed -- elements_bed's dual-family manifest row (see Prediction
        # Tabular Files' module docstring) shares one physical file across
        # two portal aliases, and enabled_families defaults to Multiome-only,
        # so this hasn't yet had a scATAC-only cluster to get wrong.
        portal_link=lambda wildcards: portal_link_alias(
            wildcards.dataset, wildcards.cluster, "multiome_powerlaw_v3", "elements_bed",
        ),
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        {{
            echo "# Source: ABC_element_list"
            echo "# Version: {params.version}"
            echo "# GenomeReference: IGVFDS0280IQAI"
            echo "# URL: https://github.com/EngreitzLab/scE2G/tree/main"
            echo "# Assays: 10x Multiome"
            echo "# SampleAgnostic: False"
            echo "# SampleTermName: {params.name}"
            echo "# SampleTermID: {params.term_id}"
            echo "# CellAnnotation: {params.summary}"
            echo "# Metadata: https://data.igvf.org/tabular-files/{params.portal_link}"
            printf '#ElementChr\\tElementStart\\tElementEnd\\tElementName\\n'
            cat {input.neighborhood_dir}/EnhancerList.bed
        }} > {output}
        """


rule bgzip_index_element_list:
    """EnhancerList.bed is already coordinate-sorted (see reformat_element_list
    above), so no sort step is needed here."""
    input:
        bed=os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_element_list.bed"),
    output:
        gz=os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_element_list.bed.gz"),
        tbi=os.path.join(RESULTS_DIR_BASE, "{dataset}", "{cluster}", "{dataset}_{cluster}_element_list.bed.gz.tbi"),
    conda:
        "../envs/filter_multiome_env.yaml"  # has bioconda::htslib (bgzip, tabix)
    resources:
        mem_mb=determine_mem_mb,
    shell:
        """
        bgzip -c {input.bed} > {output.gz}
        tabix -p bed {output.gz}
        """
