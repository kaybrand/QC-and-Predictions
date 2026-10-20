"""Shared cross-table alias references -- so a reference format is defined
exactly once, not re-derived in every table module that needs to link to
it.

prediction_set_alias() references a table not built yet (Prediction Set)
-- raises a KeyError only if actually called before it's registered,
which depends_on in each caller correctly prevents in practice (see
orchestrator.py's plan_table).

atac_fragment_alias()/rna_matrix_alias() now delegate to the real
"filtered_atac_fragment_file"/"filtered_rna_count_matrix" TableSpecs --
every caller (Prediction Tabular Files' full/elements rows, Signal
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

gene_list_alias() / gene_list_curated_set_alias() (2026-09-25): fixed,
dataset-independent constants, same "already exists, never
created/uploaded by this table" status as trained_model_file_alias() /
trained_model_set_alias() above -- added when Prediction Tabular Files'
per-(dataset,cluster) "genes" variant (content_type "gene quantifications")
was retired. Every Multiome "full" row's derived_from now references
gene_list_alias() instead of its own now-removed "...predictions_genes"
alias, and every Prediction Set's input_file_sets gains
gene_list_curated_set_alias() alongside the model Set and Principal
Pseudobulk Set.

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

per_cell_quality_report_aliases() / primary_pseudobulk_atac_fragment_aliases() /
primary_pseudobulk_h5ad_aliases() (2026-08-03): resolved, per-subsample, in
the same "anshul-kundaje:{dataset}-{cluster}-{subsample}"-based namespace as
principal_pseudobulk_set.py's _primary_pseudobulk_aliases, just with a
"-per_cell_qc_tsv"/"-fragments_tsv_gz"/"-rna_counts_mtx_h5ad" suffix naming
which file within that primary pseudobulk.

annotation_table_alias() removed (2026-08-03): the annotation table is only
referenced by primary pseudobulks (Anshul Kundaje's lab, upstream of this
pipeline) -- Filtered ATAC Fragment Files' and Filtered Matrix Files' own
derived_from never pointed at it, so their two remaining parts (the primary
pseudobulk file aliases above + the filtered_barcode_list alias) are
sufficient.

plotting_analysis_step_version_alias() (2026-07-29): resolved to the real
portal analysis step version for Filtered Barcode Lists. Updated 2026-08-25
to /analysis-step-versions/0077c8e1-f3f7-4e4e-b79e-6e6560820c9b/, replacing
/analysis-step-versions/209d5c8e-8ccb-48c5-8b51-6919b426cbcb/.

qc_guide_to_atac_fragment_analysis_step_version_alias() /
qc_guide_to_rna_matrix_analysis_step_version_alias() (2026-07-29): both
resolved to the same real portal analysis step version,
_PRINCIPAL_PSEUDOBULK_ANALYSIS_STEP_VERSION below -- kept as two distinct
functions (one per calling table's own semantics: ATAC fragment/index vs.
RNA matrix workflow) even though the value happens to coincide today.
"""

import re

from . import registry, state, subsamples
from .context import make_alias

# Shared by qc_guide_to_atac_fragment_analysis_step_version_alias and
# qc_guide_to_rna_matrix_analysis_step_version_alias -- see module docstring.
# The principal-pseudobulk analysis step version, shared by Filtered ATAC
# Fragment File, Filtered RNA Count Matrix and ATAC Index File. Updated
# 2026-08-25 to 28d9aece-0621-458c-949f-0314c8301228, replacing 9ae05eb5-ab8e-4ee0-b537-ab0ae7a1cf44.
_PRINCIPAL_PSEUDOBULK_ANALYSIS_STEP_VERSION = "/analysis-step-versions/28d9aece-0621-458c-949f-0314c8301228/"

# A different lab's (Anshul Kundaje's) namespace, not ours.
_KUNDAJE_ALIAS_PREFIX = "anshul-kundaje"

# IGVFDS-accession submissions (2026-09 onward) changed the Kundaje alias shape:
#   set:  anshul-kundaje:pseudobulk-{IGVFDS...}-{annotation}-{subsample}
#   file: ...-per_cell_qc_tsv_gz / -fragments_tsv_gz / -rna_counts_mtx_h5ad
# (old igvfN: no "pseudobulk-" prefix, "-per_cell_qc_tsv"). For those datasets the
# aliases are READ from state.db rather than formulated: files from the
# portal_files ledger by (dataset, annotation, subsample, content_type,
# file_format), sets from the raw Cell Annotation cache. Both are filled by the
# existing wholesale PseudobulkSet GET, so this adds no Portal requests. A merged
# cluster lists every constituent's files (each subsample under the annotation it
# really belongs to). igvfN datasets keep the original formula byte-for-byte --
# they are live on the Portal, out of scope, and igvf18's manifests were corrected
# by hand.
_LEDGER_RESOLVED_DATASET = re.compile(r"^IGVFDS\w{8}$")
_FILE_KINDS = {
    "per_cell_qc": ("per-cell quality report", "tsv"),
    "fragments": ("fragments", "tsv"),
    "h5ad": ("cell by gene matrix", "h5ad"),
}


def _resolves_from_ledger(ctx):
    return bool(_LEDGER_RESOLVED_DATASET.match(ctx.dataset))


def _source_annotations(ctx):
    """The cluster's own name, or a merged cluster's constituent annotations --
    the same rule derive_scopes matches Cell Annotations on."""
    from .cell_metadata import _merge_sources

    return _merge_sources(ctx.cluster_cfg, ctx.cluster)


def _require_conn(ctx, what):
    if ctx.conn is None:
        raise ValueError(f"{ctx.dataset}/{ctx.cluster}: resolving {what} needs a state.db connection (ctx.conn)")


def _ledger_file_aliases(ctx, kind):
    """One alias per (contributing subsample, constituent annotation) that has
    this kind of file; raises if a contributing subsample has none, so the row
    is reported invalid instead of carrying a dangling reference."""
    content_type, file_format = _FILE_KINDS[kind]
    _require_conn(ctx, f"{content_type} aliases")
    names = _source_annotations(ctx)
    aliases = []
    for s in subsamples.unique_subsamples(ctx):
        hits = []
        for name in names:
            found = state.portal_file_aliases(ctx.conn, ctx.dataset, name, s, content_type, file_format)
            if len(found) > 1:
                raise ValueError(
                    f"{ctx.dataset}/{name}/{s}: {len(found)} {content_type} ({file_format}) files in the "
                    f"portal_files ledger, expected 1: {found}"
                )
            hits += found
        if not hits:
            raise ValueError(
                f"{ctx.dataset}/{ctx.cluster}: subsample {s} has no {content_type} ({file_format}) file under "
                f"{names} in the portal_files ledger -- re-run download_portal_pseudobulks.py discovery"
            )
        aliases += hits
    return aliases


def primary_pseudobulk_set_aliases(ctx):
    """The primary PseudobulkSets a principal pseudobulk is built from -- for
    Principal Pseudobulk Set's input_file_sets. One per contributing subsample
    (per constituent annotation, for a merged cluster)."""
    if not _resolves_from_ledger(ctx):
        return [f"{_KUNDAJE_ALIAS_PREFIX}:{ctx.dataset}-{ctx.cluster}-{s}" for s in subsamples.unique_subsamples(ctx)]
    _require_conn(ctx, "primary pseudobulk set aliases")
    names = _source_annotations(ctx)
    aliases = []
    for s in subsamples.unique_subsamples(ctx):
        hits = [
            a for a in (f"{_KUNDAJE_ALIAS_PREFIX}:pseudobulk-{ctx.dataset}-{name}-{s}" for name in names)
            if state.primary_pseudobulk_alias_exists(ctx.conn, a)
        ]
        if not hits:
            raise ValueError(
                f"{ctx.dataset}/{ctx.cluster}: subsample {s} has no primary pseudobulk set under {names} in the "
                "Cell Annotation cache -- refresh it (downloader discovery or the driver's warm stage)"
            )
        aliases += hits
    return aliases


def prediction_set_alias(ctx):
    return registry.get("prediction_set").build_alias(ctx, "")


def bedpe_prediction_path(ctx):
    """Reaches into prediction_tabular_files' private _bedpe_path so other
    tables (e.g. bedpe_index_file) that need the bedpe file's own path
    don't duplicate that (still-unresolved) path formula."""
    from .tables.prediction_tabular_files import _bedpe_path

    return _bedpe_path(ctx)


def elements_bed_path(ctx):
    """Reaches into prediction_tabular_files' private _elements_bed_path so
    other tables (e.g. elements_bed_index_file) that need the elements BED
    file's own path don't duplicate that formula."""
    from .tables.prediction_tabular_files import _elements_bed_path

    return _elements_bed_path(ctx)


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


def gene_list_alias(ctx):
    """The fixed, externally-curated scE2G candidate gene list FILE -- for
    derived_from-shaped fields (e.g. Prediction Tabular Files' full row).
    Like trained_model_file_alias, this names an object this pipeline never
    creates or uploads, so callers need no matching depends_on entry."""
    return "jesse-engreitz:sce2g_gencode_v43_gene_list"


def gene_list_curated_set_alias(ctx):
    """The fixed FileSet gene_list_alias's FILE belongs to -- for Sets-only
    fields (Prediction Set's input_file_sets). Same "already exists, never
    created/uploaded here" status as gene_list_alias."""
    return "jesse-engreitz:sce2g_gencode_v43_gene_list_curated_set"


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
    """One per contributing subsample's primary pseudobulk -- for Filtered
    Barcode Lists' derived_from."""
    if _resolves_from_ledger(ctx):
        return _ledger_file_aliases(ctx, "per_cell_qc")
    return [
        f"{_KUNDAJE_ALIAS_PREFIX}:{ctx.dataset}-{ctx.cluster}-{s}-per_cell_qc_tsv"
        for s in subsamples.unique_subsamples(ctx)
    ]


def plotting_analysis_step_version_alias(ctx):
    return "/analysis-step-versions/0077c8e1-f3f7-4e4e-b79e-6e6560820c9b/"


def qc_guide_to_atac_fragment_analysis_step_version_alias(ctx):
    return _PRINCIPAL_PSEUDOBULK_ANALYSIS_STEP_VERSION


def primary_pseudobulk_atac_fragment_aliases(ctx):
    """The ATAC fragment FILE aliases *within* each contributing primary
    pseudobulk (distinct from the pseudobulk SET alias itself -- see
    tables/principal_pseudobulk_set.py's _primary_pseudobulk_aliases,
    "anshul-kundaje:{dataset}-{cluster}-{subsample}"). One per subsample in
    the cluster's QC-filtered barcode list (cluster_cfg["qc_guide"]) -- NOT
    every subsample queryable on the portal. Needed for Filtered ATAC
    Fragment Files' derived_from."""
    if _resolves_from_ledger(ctx):
        return _ledger_file_aliases(ctx, "fragments")
    return [
        f"{_KUNDAJE_ALIAS_PREFIX}:{ctx.dataset}-{ctx.cluster}-{s}-fragments_tsv_gz"
        for s in subsamples.unique_subsamples(ctx)
    ]


def qc_guide_to_rna_matrix_analysis_step_version_alias(ctx):
    return _PRINCIPAL_PSEUDOBULK_ANALYSIS_STEP_VERSION


def primary_pseudobulk_h5ad_aliases(ctx):
    """Same situation as primary_pseudobulk_atac_fragment_aliases, but for
    the h5ad files within each contributing primary pseudobulk -- needed
    for Filtered Matrix Files' derived_from."""
    if _resolves_from_ledger(ctx):
        return _ledger_file_aliases(ctx, "h5ad")
    return [
        f"{_KUNDAJE_ALIAS_PREFIX}:{ctx.dataset}-{ctx.cluster}-{s}-rna_counts_mtx_h5ad"
        for s in subsamples.unique_subsamples(ctx)
    ]


def rna_matrix_file_format_specifications_alias(ctx):
    """Resolved 2026-07-29: the RNA count Matrix Market tar archive's own file
    format specification Document, uploaded as a one-off via
    igvf_manifests/rna_matrix_market_file_format_spec_uploader.tsv."""
    return "jesse-engreitz:rna_matrix_market_tar_archive_file_format"


def filtered_barcode_list_file_format_specifications_alias(ctx):
    """Resolved 2026-08-03: Filtered Barcode Lists' own file format
    specification Document alias."""
    return "jesse-engreitz:filtered_barcode_membership_file_format"
