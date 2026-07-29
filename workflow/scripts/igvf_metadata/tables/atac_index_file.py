"""ATAC Index File -- object_type "tabular_file" (TODO: confirm portal
profile id), scope "cluster" -- one per (dataset, cluster), NOT split by
model: aliases/file_set/derived_from/submitted_file_name all have no
{Family}/model component, since the ATAC fragment file (and its index) is
shared across both Multiome and scATAC predictions for a cluster.

Resolved 2026-07-20: description is family-agnostic --
"ATAC fragment file index used as input for E2G predictions for {dataset}
{cluster} cells" -- no {Family}/model wording at all, since principal
pseudobulks now feed many E2G models, not just scE2G. This also matches
ctx.model being None here (same as Principal Pseudobulk Set/Filtered Barcode
Lists/QC Documents).

submitted_file_name mirrors workflow/rules/common.smk's own
multiome_data_dir(dataset)/cluster convention (see context.py's
Context.multiome_data_cluster_dir) -- same Synapse-side filtered_data
location the existing Synapse manifest system already uploads
atac_fragments_{dataset}_{cluster}.tsv.gz from.

analysis_step_version (refs.qc_guide_to_atac_fragment_analysis_step_version_alias,
resolved 2026-07-29) is /analysis-step-versions/9ae05eb5-ab8e-4ee0-b537-ab0ae7a1cf44/.

reference_files removed (2026-07-20 feedback): not a submittable field for
type IndexFile (content_type=index).
"""

import os

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "atac_index_file"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_ATAC_fragment_file_index")


def _path(ctx):
    return os.path.join(ctx.multiome_data_cluster_dir, f"atac_fragments_{ctx.dataset}_{ctx.cluster}.tsv.gz.tbi")


def _enabled(ctx):
    return os.path.exists(_path(ctx))


def _row(ctx):
    return {
        "file_format": "tbi",
        "content_type": "index",
        "analysis_step_version": refs.qc_guide_to_atac_fragment_analysis_step_version_alias(ctx),
        "derived_from": refs.atac_fragment_alias(ctx),
        "description": (
            f"ATAC fragment file index used as input for E2G predictions for {ctx.dataset} {ctx.cluster} cells"
        ),
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
        constant_fields={"controlled_access": False, "derived_manually": False},
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="",
                build_row=_row,
                enabled=_enabled,
                depends_on=lambda ctx: [("principal_pseudobulk_set", ""), ("filtered_atac_fragment_file", "")],
            ),
        ],
    )
)
