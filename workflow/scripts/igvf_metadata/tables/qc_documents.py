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
"{QC_PSEUDOBULKS_PLOTS_DIR}/{dataset}/{cluster}/qc_thresholds.tsv".

Corrected 2026-08-03: QC_PSEUDOBULKS_PLOTS_DIR is the `main`-branch
QC_pseudobulks worktree's plots/ dir (a sibling worktree of this repo, per
`git worktree list` -- NOT a separate "QC_clusters" project, which doesn't
exist anywhere on disk). The constant was previously misnamed
QC_CLUSTERS_PLOTS_DIR and pointed at a nonexistent path, silently disabling
this table (and everything depends_on-ing it) for every dataset, not just
one -- caught when a dry run for igvf4 showed zero QC_documents rows ever
enabled despite qc_thresholds.tsv actually existing on disk.

Still open, per the user's own flag: whether this needs a .txt extension
instead of .tsv, pending what the portal's attachment validation actually
accepts.

derived_manually removed from constant_fields (2026-07-20 feedback): not a
submittable field for Document objects.
"""

import os

from .. import registry
from ..context import QC_PSEUDOBULKS_WDIR, make_alias

TABLE_NAME = "QC_documents"

# The `main`-branch QC_pseudobulks worktree's plots/ dir -- a sibling
# worktree of this repo (this repo is the igvf-portal-submission worktree),
# NOT this checkout's own plots/ dir. Derived from context.QC_PSEUDOBULKS_WDIR
# (the one place that worktree's root is defined) rather than its own literal,
# so there's exactly one path to change if that worktree ever moves --
# same reasoning as context.py's own multiome_data_cluster_dir.
QC_PSEUDOBULKS_PLOTS_DIR = os.path.join(QC_PSEUDOBULKS_WDIR, "plots")


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "QC_thresholds")


def _attachment_path(ctx):
    return os.path.join(QC_PSEUDOBULKS_PLOTS_DIR, ctx.dataset, ctx.cluster, "qc_thresholds.tsv")


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
