"""Shared cross-table alias references -- so a reference format is defined
exactly once, not re-derived in every table module that needs to link to
it.

prediction_set_alias()/atac_fragment_alias()/rna_matrix_alias() reference
tables not built yet (Prediction Set, Filtered ATAC Fragment Files,
Filtered Matrix Files) -- they raise a KeyError only if actually called
before those tables are registered, which depends_on in each caller
correctly prevents in practice (see orchestrator.py's plan_table).

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

primary_pseudobulk_metadata() is a true stub: cell_type/cell_qualifier on
Principal Pseudobulk Set need a live GET against data.igvf.org's
PseudobulkSet multireport endpoint, "indexed to the rows we need" -- but
which field to index/match on (aliases? cell_annotation? something else?)
hasn't been answered yet. Fill in here once it is; this is also where
request caching should live (one multireport fetch per run, not per
cluster).

qc_thresholds_document_alias() is confirmed and deterministic (used by
both Principal Pseudobulk Set's and a later table's `documents`) --
scope="cluster", not cluster_model, same as Principal Pseudobulk Set.
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
    """Should return {"cl_id": ..., "cell_qualifier": ...} looked up from
    data.igvf.org's PseudobulkSet multireport, matched to this (dataset,
    cluster)'s primary pseudobulk sets. Raises until the indexing/matching
    key is confirmed -- see this module's docstring."""
    raise NotImplementedError(
        "cell_type/cell_qualifier need a live PseudobulkSet multireport lookup, indexed "
        "on a key that hasn't been confirmed yet (aliases? cell_annotation? something "
        "else?) -- fill in here once it is."
    )


def qc_thresholds_document_alias(ctx):
    """Cluster-scoped (not cluster_model), like Principal Pseudobulk Set."""
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "QC_thresholds")
