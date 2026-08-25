"""Filtered Matrix Files -- object_type "matrix_file" (confirmed 2026-08-03
-- not tabular_file as previously guessed), scope "cluster" -- one per (dataset, cluster),
scE2G-Family-agnostic like everything else whose file_set is the Principal
Pseudobulk Set (same as Filtered ATAC Fragment Files).

This is the table refs.rna_matrix_alias delegates to -- see that
function's docstring for why every OTHER table calling it (Prediction
Tabular Files' full (Multiome-only) and genes rows) needed
("filtered_rna_count_matrix", "") added to their own depends_on once this
table became real.

derived_from mirrors Filtered ATAC Fragment Files' two-part shape
(2026-08-03: the annotation table dropped -- only primary pseudobulks point
at it, not files derived from them):
  1. contributing h5ads in primary pseudobulks
     (refs.primary_pseudobulk_h5ad_aliases) -- one per subsample in this
     cluster's QC-filtered barcode list, resolved 2026-08-03.
  2. jesse-engreitz:{dataset}_{cluster}_filtered_barcode_list -- known,
     references the real filtered_barcode_list table directly.

file_format_specifications (refs.rna_matrix_file_format_specifications_alias,
resolved 2026-07-29) is jesse-engreitz:rna_matrix_market_tar_archive_file_format
-- a one-off Document uploaded via
igvf_manifests/rna_matrix_market_file_format_spec_uploader.tsv, unlike
Filtered ATAC Fragment Files' already-known /documents/<uuid>/ path.

analysis_step_version (refs.qc_guide_to_rna_matrix_analysis_step_version_alias,
resolved 2026-07-29) is /analysis-step-versions/9ae05eb5-ab8e-4ee0-b537-ab0ae7a1cf44/
-- the same analysis step version as ATAC's, per portal manager feedback.

description reads "Matrix Market" (confirmed 2026-07-16) -- the standard
sparse-matrix file format name, matching file_format "tar" bundling a
10x-style matrix.mtx/barcodes/features archive.

submitted_file_name mirrors Filtered ATAC Fragment Files' own
multiome_data_cluster_dir convention, with its own filename.

reference_files (2026-07-20 feedback): changed to /reference-files/IGVFFI9573KOZR/
-- a filtered RNA count matrix pipeline input, not the primary pseudobulk
set's own MatrixFiles' reference (which was wrong here, same class of
mistake as Filtered ATAC Fragment Files' reference_files).

controlled_access removed (2026-07-20 feedback): not a submittable field
for object type MatrixFile.
"""

import os

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "filtered_rna_count_matrix"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_RNA_count_matrix")


def _path(ctx):
    return os.path.join(ctx.multiome_data_cluster_dir, f"rna_count_matrix_{ctx.dataset}_{ctx.cluster}.tar.gz")


def _required_paths(ctx):
    return [_path(ctx)]


def _row(ctx):
    derived_from_parts = [
        *refs.primary_pseudobulk_h5ad_aliases(ctx),
        registry.get("filtered_barcode_list").build_alias(ctx, ""),
    ]
    return {
        "file_format": "tar",
        "content_type": "cell by gene matrix",
        "reference_files": ["/reference-files/IGVFFI9573KOZR/"],
        "analysis_step_version": refs.qc_guide_to_rna_matrix_analysis_step_version_alias(ctx),
        "derived_from": ",".join(derived_from_parts),
        "description": f"Filtered Matrix Market file containing RNA transcripts for cells annotated as {ctx.cluster}",
        "file_format_specifications": [refs.rna_matrix_file_format_specifications_alias(ctx)],
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
        object_type="matrix_file",  # confirmed 2026-08-03 -- not tabular_file as previously guessed
        scope="cluster",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "file_format", "file_set", "content_type"],
        constant_fields={"filtered": True, "derived_manually": False},
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="",
                build_row=_row,
                required_paths=_required_paths,
                depends_on=lambda ctx: [("principal_pseudobulk_set", ""), ("filtered_barcode_list", "")],
            ),
        ],
    )
)
