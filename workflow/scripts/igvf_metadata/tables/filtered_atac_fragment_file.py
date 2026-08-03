"""Filtered ATAC Fragment Files -- object_type "tabular_file" (TODO:
confirm portal profile id), scope "cluster" -- one per (dataset, cluster),
scE2G-Family-agnostic like everything else whose file_set is the Principal
Pseudobulk Set (confirmed 2026-07-16).

This is the table refs.atac_fragment_alias delegates to -- see that
function's docstring for why every OTHER table calling it (Prediction
Tabular Files' full/elements rows, Signal Files' derived_from, ATAC Index
File's derived_from) needed ("filtered_atac_fragment_file", "") added to
their own depends_on once this table became real.

derived_from has two parts (2026-08-03: the annotation table dropped -- only
primary pseudobulks point at it, not files derived from them):
  1. contributing ATAC fragment files in primary pseudobulks
     (refs.primary_pseudobulk_atac_fragment_aliases) -- one per subsample in
     this cluster's QC-filtered barcode list, resolved 2026-08-03.
  2. jesse-engreitz:{dataset}_{cluster}_filtered_barcode_list -- known,
     references the real filtered_barcode_list table directly.

analysis_step_version reuses the same resolved reference as ATAC Index File
(refs.qc_guide_to_atac_fragment_analysis_step_version_alias, resolved
2026-07-29 to /analysis-step-versions/9ae05eb5-ab8e-4ee0-b537-ab0ae7a1cf44/)
-- confirmed "cluster agnostic," same object either way.

file_format_specifications is a raw portal path (/documents/<uuid>/), not
an alias -- an existing Document this pipeline doesn't create itself, so
just a literal constant, no depends_on needed for it.

submitted_file_name mirrors ATAC Index File's own multiome_data_cluster_dir
convention, minus the ".tbi" suffix.

reference_files (2026-07-20 feedback): IGVFFI0653VCGH/IGVFFI9573KOZR were
wrong here -- those are the primary pseudobulk set ATAC fragment files' own
reference_files. This table's reference_files must instead name the file(s)
the code that GENERATES this data file directly used -- i.e. the
chrom-sizes file workflow/rules/filter_pseudobulks.smk passes as
filter_atac_fragments.py's --chrom-sizes
(reference/IGVF.DACC.GRCh38.chrom.sizes.tsv). Its IGVF accession/alias isn't
confirmed yet, so _CHROM_SIZES_REFERENCE_FILE below is an explicit
placeholder pending that.
"""

import os

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "filtered_atac_fragment_file"

# PLACEHOLDER: real IGVF accession/alias for reference/IGVF.DACC.GRCh38.chrom.sizes.tsv
# (the chrom-sizes file filter_atac_fragments.py's --chrom-sizes flag reads) not
# confirmed yet -- fill in once known, see module docstring.
_CHROM_SIZES_REFERENCE_FILE = "PLACEHOLDER_IGVF_alias_for_reference_IGVF.DACC.GRCh38.chrom.sizes.tsv"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_ATAC_fragment_file")


def _path(ctx):
    return os.path.join(ctx.multiome_data_cluster_dir, f"atac_fragments_{ctx.dataset}_{ctx.cluster}.tsv.gz")


def _enabled(ctx):
    return os.path.exists(_path(ctx))


def _row(ctx):
    derived_from_parts = [
        *refs.primary_pseudobulk_atac_fragment_aliases(ctx),
        registry.get("filtered_barcode_list").build_alias(ctx, ""),
    ]
    return {
        "file_format": "bed",
        "file_format_type": "bed3+",
        "content_type": "fragments",
        "file_format_specifications": "/documents/db2a6dd0-cc1d-439e-a610-f9f1d04cfd82/",
        "reference_files": [_CHROM_SIZES_REFERENCE_FILE],
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
