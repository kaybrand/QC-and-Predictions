"""Filtered Matrix Files -- object_type "tabular_file" (TODO: confirm
portal profile id), scope "cluster" -- one per (dataset, cluster),
scE2G-Family-agnostic like everything else whose file_set is the Principal
Pseudobulk Set (same as Filtered ATAC Fragment Files).

This is the table refs.rna_matrix_alias delegates to -- see that
function's docstring for why every OTHER table calling it (Prediction
Tabular Files' full (Multiome-only) and genes rows) needed
("filtered_rna_count_matrix", "") added to their own depends_on once this
table became real.

derived_from mirrors Filtered ATAC Fragment Files' three-part shape:
  1. contributing h5ads in primary pseudobulks -- TBD, computationally
     sourceable, mechanism not yet defined (refs.primary_pseudobulk_h5ad_aliases).
  2. the annotation table -- TBD, same stub Filtered ATAC Fragment Files
     uses (refs.annotation_table_alias).
  3. jesse-engreitz:{dataset}_{cluster}_filtered_barcode_list -- known,
     references the real filtered_barcode_list table directly.

file_format_specifications is explicitly "{alias TBD}" (unlike Filtered
ATAC Fragment Files' already-known /documents/<uuid>/ path) --
refs.rna_matrix_file_format_specifications_alias, a stub.

analysis_step_version is its own stub, distinct from ATAC's
(refs.qc_guide_to_rna_matrix_analysis_step_version_alias) -- a different
workflow (QC guide -> RNA count matrix, not -> ATAC fragment file).

description reads "Matrix Market" (confirmed 2026-07-16) -- the standard
sparse-matrix file format name, matching file_format "tar" bundling a
10x-style matrix.mtx/barcodes/features archive.

submitted_file_name mirrors Filtered ATAC Fragment Files' own
multiome_data_cluster_dir convention, with its own filename.
"""

import os

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "filtered_rna_count_matrix"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_RNA_count_matrix")


def _path(ctx):
    return os.path.join(ctx.multiome_data_cluster_dir, f"rna_count_matrix_{ctx.dataset}_{ctx.cluster}.tar.gz")


def _enabled(ctx):
    return os.path.exists(_path(ctx))


def _row(ctx):
    derived_from_parts = [
        *refs.primary_pseudobulk_h5ad_aliases(ctx),
        refs.annotation_table_alias(ctx),
        registry.get("filtered_barcode_list").build_alias(ctx, ""),
    ]
    return {
        "file_format": "tar",
        "content_type": "cell by gene matrix",
        "reference_files": ["/reference-files/IGVFFI9561BASO/"],
        "analysis_step_version": refs.qc_guide_to_rna_matrix_analysis_step_version_alias(ctx),
        "derived_from": ",".join(derived_from_parts),
        "description": f"Filtered Matrix Market file containing RNA transcripts for cells annotated as {ctx.cluster}",
        "file_format_specifications": refs.rna_matrix_file_format_specifications_alias(ctx),
        "submitted_file_name": _path(ctx),
    }


def _scope_fields(ctx):
    return {
        "md5sum": None,  # left blank; igvf_utils computes+fills it from submitted_file_name
        "file_set": refs.principal_pseudobulk_set_alias(ctx),
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="tabular_file",  # TODO: confirm actual portal profile id
        scope="cluster",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "file_format", "file_set", "content_type", "controlled_access"],
        constant_fields={"filtered": True, "controlled_access": False, "derived_manually": False},
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="",
                build_row=_row,
                enabled=_enabled,
                depends_on=lambda ctx: [("principal_pseudobulk_set", ""), ("filtered_barcode_list", "")],
            ),
        ],
    )
)
