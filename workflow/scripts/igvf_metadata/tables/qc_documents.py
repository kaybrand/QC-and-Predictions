"""QC Documents -- the QC_thresholds Document referenced (as a forward
declaration, via depends_on) by both Principal Pseudobulk Set and Filtered
Barcode Lists' `documents` fields. object_type "document" (per igvf_utils'
own docstring, which uses exactly this profile as its example for the
`attachment` "path" trick). Scope "cluster" -- one per (dataset, cluster),
matching every existing reference to refs.qc_thresholds_document_alias.

Named "QC_documents" (not just "documents") since more Document type(s) are
expected to be added later (2026-07-20 feedback) -- keeping this table's own
name scoped to what it actually is now avoids a collision/rename later.

attachment MUST be {"path": "<local file path>"} -- not a raw path string.
This is igvf_utils.connection.Connection's own documented special case:
given this single-key dict, it builds the real `attachment` schema object
itself (base64-encoding the file, computing its own metadata) rather than
us constructing that object. Serializes correctly through the existing TSV
writer with no changes needed: portal_client._tsv_cell already JSON-dumps
dict values, matching iu_register.py's own requirement that an
object-typed TSV field be valid JSON.

Confirmed 2026-07-16: attachment path is
"{QC_CLUSTERS_PLOTS_DIR}/{dataset}/{cluster}/qc_thresholds.tsv" --
QC_CLUSTERS_PLOTS_DIR is a genuinely separate directory from this repo
(QC_pseudobulks), not a typo.

Still open, per the user's own flag: whether this needs a .txt extension
instead of .tsv, pending what the portal's attachment validation actually
accepts.

derived_manually removed from constant_fields (2026-07-20 feedback): not a
submittable field for Document objects.
"""

import os

from .. import registry
from ..context import make_alias

TABLE_NAME = "QC_documents"

# A sibling project directory, NOT this repo's own plots/ dir -- confirmed
# 2026-07-16. Kept as a module constant (rather than hardcoded inline)
# so it's easy to find/override; could be promoted into IgvfConfig if it
# ever needs to vary per environment.
QC_CLUSTERS_PLOTS_DIR = "/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_clusters/plots"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "QC_thresholds")


def _attachment_path(ctx):
    return os.path.join(QC_CLUSTERS_PLOTS_DIR, ctx.dataset, ctx.cluster, "qc_thresholds.tsv")


def _enabled(ctx):
    return os.path.exists(_attachment_path(ctx))


def _row(ctx):
    return {
        "description": f"Quality Control thresholds applied to each cell in {ctx.cluster}",
        "document_type": "pipeline parameters",
        "attachment": {"path": _attachment_path(ctx)},
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="document",
        scope="cluster",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "document_type", "attachment"],
        variants=[
            registry.VariantSpec(name="", build_row=_row, enabled=_enabled),
        ],
    )
)
