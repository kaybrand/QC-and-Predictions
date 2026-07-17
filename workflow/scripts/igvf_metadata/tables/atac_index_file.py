"""ATAC Index File -- object_type "tabular_file" (TODO: confirm portal
profile id), scope "cluster" -- one per (dataset, cluster), NOT split by
model: aliases/file_set/derived_from/submitted_file_name all have no
{Family}/model component, since the ATAC fragment file (and its index) is
shared across both Multiome and scATAC predictions for a cluster.

OPEN ISSUE (2026-07-16, unconfirmed): the given `description` template --
"ATAC fragment file index used as input for scE2G_{Family} predictions..."
-- references {Family}, but this table has no per-model context (ctx.model
is None here, same as Principal Pseudobulk Set/Filtered Barcode
Lists/Documents). Calling family(ctx.model) would raise (family() has no
mapping for None). Implemented as listing every family actually configured
for this cluster (from cluster_cfg["models"]), e.g. "scE2G Multiome, scE2G
scATAC predictions for..." when both are configured -- confirm this is the
intended wording, since the original text implies exactly one family.

submitted_file_name mirrors workflow/rules/common.smk's own
multiome_data_dir(dataset)/cluster convention (see context.py's
Context.multiome_data_cluster_dir) -- same Synapse-side filtered_data
location the existing Synapse manifest system already uploads
atac_fragments_{dataset}_{cluster}.tsv.gz from.

analysis_step_version is a true stub (refs.qc_guide_to_atac_fragment_analysis_step_version_alias)
-- "to be created," no alias given yet.
"""

import os

from .. import refs, registry
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "atac_index_file"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_ATAC_fragment_file_index")


def _path(ctx):
    return os.path.join(ctx.multiome_data_cluster_dir, f"atac_fragments_{ctx.dataset}_{ctx.cluster}.tsv.gz.tbi")


def _enabled(ctx):
    return os.path.exists(_path(ctx))


def _family_list_text(ctx):
    families = []
    for model in ctx.cluster_cfg["models"]:
        fam = family(model)
        if fam not in families:
            families.append(fam)
    return ", ".join(f"scE2G {f}" for f in families)


def _row(ctx):
    return {
        "file_format": "tbi",
        "content_type": "index",
        "reference_files": ["IGVFFI0653VCGH", "IGVFFI9573KOZR"],
        "analysis_step_version": refs.qc_guide_to_atac_fragment_analysis_step_version_alias(ctx),
        "derived_from": refs.atac_fragment_alias(ctx),
        "description": (
            f"ATAC fragment file index used as input for {_family_list_text(ctx)} predictions "
            f"for {ctx.dataset} {ctx.cluster} cells"
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
