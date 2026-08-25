"""Filtered Barcode Lists -- object_type "tabular_file" (confirmed
2026-08-03), scope "cluster" (one per dataset/cluster, not split by
model -- same scope as Principal Pseudobulk Set, which this table's
file_set points at).

This is the object formerly informally called "the QC guide" in this
package's own local-file terminology (cluster_cfg["qc_guide"]) -- that
local path IS this table's submitted_file_name, "USUALLY" the standard
plots/{dataset}/{cluster}/filtered_barcodes_with_subsamples.tsv.gz path but
already free to point anywhere for a bespoke-filtered cluster, since
cluster_cfg["qc_guide"] is itself a per-cluster config value with no fixed
naming assumption. No new path builder needed.

Confirmed field-by-field 2026-07-16. One still-open stub, routed through
refs.py so there's exactly one place to fill in:
  - derived_from: the per-cluster per-cell QC metric files (content_type
    "per-cell quality report" in the pseudobulks) -- no alias formula given
    yet (refs.per_cell_quality_report_aliases).
analysis_step_version (refs.plotting_analysis_step_version_alias) is
/analysis-step-versions/0077c8e1-f3f7-4e4e-b79e-6e6560820c9b/ as of 2026-08-25
(previously .../209d5c8e-8ccb-48c5-8b51-6919b426cbcb/).

content_type is "barcode to sample mapping" (2026-07-20 feedback, explicitly
flagged by the user as subject to change -- was "filtered barcode list").

submitter_comment is genuinely optional/manual (per-cluster free text
"explaining any additional filtering not described by the QC document") --
read from cluster_cfg.get("submitter_comment"), omitted when absent since
it's not in required_columns. This already is the config field 2026-07-20
feedback asked for ("optional submitter comments describing custom filter
rules not covered by the QC_thresholds document") -- no further change
needed, a user just sets cluster_cfg["submitter_comment"] per cluster.

documents references the same QC_thresholds Document as Principal
Pseudobulk Set's own documents field (refs.qc_thresholds_document_alias) --
not yet a registered table, so depends_on lists it, same pattern as
elsewhere.
"""

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "filtered_barcode_list"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_barcode_list")


def _row(ctx):
    return {
        "file_format": "tsv",
        "content_type": "barcode to sample mapping",
        "documents": [refs.qc_thresholds_document_alias(ctx)],
        "file_format_specifications": [refs.filtered_barcode_list_file_format_specifications_alias(ctx)],
        "derived_from": ",".join(refs.per_cell_quality_report_aliases(ctx)),
        "description": (
            f"Filtered list of 16-bp cell barcodes defining membership in the {ctx.cluster}; "
            "see QC thresholds for a record of the filters applied"
        ),
        "submitter_comment": ctx.cluster_cfg.get("submitter_comment"),
        "analysis_step_version": refs.plotting_analysis_step_version_alias(ctx),
        "submitted_file_name": ctx.cluster_cfg["qc_guide"],  # "USUALLY" the standard path; already free to be bespoke
    }


def _scope_fields(ctx):
    return {
        "md5sum": None,  # left blank; igvf_utils computes+fills it from submitted_file_name
        "file_set": refs.principal_pseudobulk_set_alias(ctx),
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="tabular_file",  # confirmed 2026-08-03
        scope="cluster",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "file_format", "file_set", "content_type", "controlled_access"],
        constant_fields={"controlled_access": False, "filtered": True, "derived_manually": False},
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="",
                build_row=_row,
                depends_on=lambda ctx: [("principal_pseudobulk_set", ""), ("QC_documents", "")],
            ),
        ],
    )
)
