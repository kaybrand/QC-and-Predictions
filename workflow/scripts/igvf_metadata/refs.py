"""Shared cross-table alias references -- so a reference format is defined
exactly once, not re-derived in every table module that needs to link to
it.

prediction_set_alias() references a table not built yet (Prediction Set)
-- raises a KeyError only if actually called before it's registered,
which depends_on in each caller correctly prevents in practice (see
orchestrator.py's plan_table).

atac_fragment_alias()/rna_matrix_alias() now delegate to the real
"filtered_atac_fragment_file"/"filtered_rna_count_matrix" TableSpecs --
every caller (Prediction Tabular Files' full/elements/genes rows, Signal
Files' derived_from, ATAC Index File's derived_from) needed the matching
("filtered_atac_fragment_file"|"filtered_rna_count_matrix", "") added to
their own depends_on once each table became real, same reasoning as
prediction_set -> principal_pseudobulk_set.

trained_model_file_alias() / trained_model_set_alias() name two different
objects, not one ambiguous alias. "scE2G_{Family}_trained_model" is the
model FILE (has submitted_file_name; used by derived_from-shaped fields,
e.g. Prediction Tabular Files' full row). "scE2G_{Family}_model" is the
model SET (referenced via some file's file_set; used by Sets-only fields,
e.g. Prediction Set's input_file_sets).

principal_pseudobulk_set_alias() references a table not built yet
(Principal Pseudobulk Set) -- same "raises only if called before
registered, which depends_on prevents" situation as prediction_set_alias
above.

primary_pseudobulk_metadata() is resolved (2026-07-21): delegates to
cell_metadata.get_metadata_for, which reads a 24h-TTL cache of one
PseudobulkSet-multireport GET per pipeline trigger (populated by
cell_metadata.refresh_if_stale, called once at the top of
orchestrator.run() -- one fetch per run, not per cluster, as originally
planned here). See cell_metadata.py's module docstring for the full
matching mechanism (join by `samples` identity against local subsample
sets, not alias parsing).

qc_thresholds_document_alias() is confirmed and deterministic (used by
both Principal Pseudobulk Set's and Filtered Barcode Lists' `documents`) --
scope="cluster", not cluster_model, same as Principal Pseudobulk Set.

per_cell_quality_report_aliases() / plotting_analysis_step_version_alias()
are true stubs -- Filtered Barcode Lists' derived_from (the per-cluster
per-cell QC metric files, content_type "per-cell quality report") and
analysis_step_version ("alias of plotting scripts analysis step version,
to be created") don't have formulas yet.
"""

from . import registry
from .context import make_alias


def prediction_set_alias(ctx):
    return registry.get("prediction_set").build_alias(ctx, "")


def bedpe_prediction_path(ctx):
    """Reaches into prediction_tabular_files' private _bedpe_path so other
    tables (e.g. bedpe_index_file) that need the bedpe file's own path
    don't duplicate that (still-unresolved) path formula."""
    from .tables.prediction_tabular_files import _bedpe_path

    return _bedpe_path(ctx)


def atac_fragment_alias(ctx):
    """The ATAC fragment FILE (submitted_file_name) -- for derived_from-shaped
    fields. Callers must also depends_on ("filtered_atac_fragment_file", "")."""
    return registry.get("filtered_atac_fragment_file").build_alias(ctx, "")


def rna_matrix_alias(ctx):
    """The RNA count matrix FILE (submitted_file_name) -- for derived_from-shaped
    fields. Callers must also depends_on ("filtered_rna_count_matrix", "")."""
    return registry.get("filtered_rna_count_matrix").build_alias(ctx, "")


def trained_model_file_alias(ctx):
    """The model FILE -- for derived_from-shaped fields (e.g. Prediction
    Tabular Files' full row)."""
    from .tables.prediction_tabular_files import family

    return make_alias(ctx.igvf, "scE2G", family(ctx.model), "trained_model")


def trained_model_set_alias(ctx):
    """The model SET -- for Sets-only fields (e.g. Prediction Set's
    input_file_sets)."""
    from .tables.prediction_tabular_files import family

    return make_alias(ctx.igvf, "scE2G", family(ctx.model), "model")


def principal_pseudobulk_set_alias(ctx):
    return registry.get("principal_pseudobulk_set").build_alias(ctx, "")


def primary_pseudobulk_metadata(ctx):
    """Returns {"cl_id": ..., "cell_qualifier": ...} from the cell_metadata
    cache -- see this module's docstring and cell_metadata.py."""
    from . import cell_metadata

    return cell_metadata.get_metadata_for(ctx)


def qc_thresholds_document_alias(ctx):
    """Cluster-scoped (not cluster_model), like Principal Pseudobulk Set."""
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "QC_thresholds")


def per_cell_quality_report_aliases(ctx):
    raise NotImplementedError(
        "Filtered Barcode Lists' derived_from names the per-cluster per-cell QC metric "
        "files (content_type 'per-cell quality report' in the pseudobulks) with no alias "
        "formula given yet -- fill in here once it is."
    )


def plotting_analysis_step_version_alias(ctx):
    raise NotImplementedError(
        "Filtered Barcode Lists' analysis_step_version names 'alias of plotting scripts "
        "analysis step version, to be created' -- not created/aliased yet."
    )


def qc_guide_to_atac_fragment_analysis_step_version_alias(ctx):
    raise NotImplementedError(
        "ATAC Index File's analysis_step_version names 'alias of QC guide to ATAC "
        "fragment file workflow version, to be created' -- not created/aliased yet."
    )


def annotation_table_alias(ctx):
    """Confirmed TBD: the annotation table's alias exists (or will exist)
    on the portal and can be computationally sourced -- but the mechanism
    isn't defined yet. Needed for Filtered ATAC Fragment Files'
    derived_from."""
    raise NotImplementedError(
        "Filtered ATAC Fragment Files' derived_from names 'the annotation table' -- "
        "confirmed to exist on the portal and be computationally sourceable, but no "
        "lookup mechanism/formula given yet."
    )


def primary_pseudobulk_atac_fragment_aliases(ctx):
    """Confirmed TBD: the ATAC fragment FILE aliases *within* each
    contributing primary pseudobulk (distinct from the pseudobulk SET
    alias itself, which IS known -- see
    tables/principal_pseudobulk_set.py's _primary_pseudobulk_aliases,
    "anshul-kundaje:{dataset}-{cluster}-{subsample}"). Also confirmed
    computationally sourceable, mechanism not defined yet. Needed for
    Filtered ATAC Fragment Files' derived_from."""
    raise NotImplementedError(
        "Filtered ATAC Fragment Files' derived_from names 'contributing ATAC fragment "
        "files in primary pseudobulks' -- confirmed to exist on the portal and be "
        "computationally sourceable, but no lookup mechanism/formula given yet."
    )


def qc_guide_to_rna_matrix_analysis_step_version_alias(ctx):
    raise NotImplementedError(
        "Filtered Matrix Files' analysis_step_version names 'alias of QC guide to RNA "
        "count matrix workflow version, to be created' -- not created/aliased yet."
    )


def primary_pseudobulk_h5ad_aliases(ctx):
    """Same situation as primary_pseudobulk_atac_fragment_aliases, but for
    the h5ad files within each contributing primary pseudobulk -- needed
    for Filtered Matrix Files' derived_from."""
    raise NotImplementedError(
        "Filtered Matrix Files' derived_from names 'contributing h5ads in primary "
        "pseudobulks' -- confirmed to exist on the portal and be computationally "
        "sourceable, but no lookup mechanism/formula given yet."
    )


def rna_matrix_file_format_specifications_alias(ctx):
    raise NotImplementedError(
        "Filtered Matrix Files' file_format_specifications is explicitly '{alias TBD}' "
        "-- not decided/given yet."
    )
