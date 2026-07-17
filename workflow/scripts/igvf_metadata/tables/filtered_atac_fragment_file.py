"""Filtered ATAC Fragment Files -- object_type "tabular_file" (TODO:
confirm portal profile id), scope "cluster" -- one per (dataset, cluster),
scE2G-Family-agnostic like everything else whose file_set is the Principal
Pseudobulk Set (confirmed 2026-07-16).

This is the table refs.atac_fragment_alias delegates to -- see that
function's docstring for why every OTHER table calling it (Prediction
Tabular Files' full/elements rows, Signal Files' derived_from, ATAC Index
File's derived_from) needed ("filtered_atac_fragment_file", "") added to
their own depends_on once this table became real.

derived_from has three parts, per the 2026-07-16 conversation:
  1. contributing ATAC fragment files in primary pseudobulks -- TBD, exists
     (or will exist) on the portal, computationally sourceable, mechanism
     not yet defined (refs.primary_pseudobulk_atac_fragment_aliases, stub).
  2. the annotation table -- TBD, same situation
     (refs.annotation_table_alias, stub).
  3. jesse-engreitz:{dataset}_{cluster}_filtered_barcode_list -- known,
     references the real filtered_barcode_list table directly.
Because (1) and (2) raise, this table's own rows never reach 'uploaded'
until those lookup mechanisms are built -- which correctly keeps every
downstream table (everything depends_on-ing this one) deferred too.

analysis_step_version reuses the same stub as ATAC Index File
(refs.qc_guide_to_atac_fragment_analysis_step_version_alias) -- confirmed
"cluster agnostic," to be manually produced later, same object either way.

file_format_specifications is a raw portal path (/documents/<uuid>/), not
an alias -- an existing Document this pipeline doesn't create itself, so
just a literal constant, no depends_on needed for it.

submitted_file_name mirrors ATAC Index File's own multiome_data_cluster_dir
convention, minus the ".tbi" suffix.
"""

import os

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "filtered_atac_fragment_file"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_ATAC_fragment_file")


def _path(ctx):
    return os.path.join(ctx.multiome_data_cluster_dir, f"atac_fragments_{ctx.dataset}_{ctx.cluster}.tsv.gz")


def _enabled(ctx):
    return os.path.exists(_path(ctx))


def _row(ctx):
    derived_from_parts = [
        *refs.primary_pseudobulk_atac_fragment_aliases(ctx),
        refs.annotation_table_alias(ctx),
        registry.get("filtered_barcode_list").build_alias(ctx, ""),
    ]
    return {
        "file_format": "bed",
        "file_format_type": "bed3+",
        "content_type": "fragments",
        "file_format_specifications": "/documents/db2a6dd0-cc1d-439e-a610-f9f1d04cfd82/",
        "reference_files": ["IGVFFI0653VCGH", "IGVFFI9573KOZR"],
        "analysis_step_version": refs.qc_guide_to_atac_fragment_analysis_step_version_alias(ctx),
        "derived_from": ",".join(derived_from_parts),
        "description": f"Filtered ATAC fragment file containing reads from cells annotated as {ctx.cluster}",
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
