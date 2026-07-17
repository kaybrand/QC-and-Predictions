"""Prediction Set -- the FileSet object binding together every scE2G model
output for one (dataset, cluster, model) on the IGVF Data Portal. Every
other table's file_set points at this table's alias; see
igvf_metadata.refs.prediction_set_alias, which delegates here.

input_file_sets/derived_from are both fully resolved: exactly
[trained_model_set_alias, principal_pseudobulk_set_alias]. _row() itself no
longer raises -- but depends_on lists ("principal_pseudobulk_set", "") so
real uploads still correctly wait for that object to actually exist on the
portal before this one links to it.

Still open: samples (computed correctly -- unique "subsamples" column
values from the cluster's local barcode-list file -- but unconfirmed
whether the portal wants these raw labels directly or references to actual
Sample objects).
"""

import csv
import gzip

from .. import refs, registry
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "prediction_set"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "predictions")


def _unique_subsamples(ctx):
    """The list of unique values in the cluster's local barcode-list file's
    "subsamples" column -- one such file defines the barcodes present in
    this cluster."""
    with gzip.open(ctx.cluster_cfg["qc_guide"], "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return sorted({row["subsamples"] for row in reader})


def _row(ctx):
    # derived_from: same two Sets as input_file_sets, but as a pre-joined
    # no-spaces comma string, matching every other derived_from field in
    # this package (not a Python list left for the TSV writer to join).
    derived_from_parts = [refs.trained_model_set_alias(ctx), refs.principal_pseudobulk_set_alias(ctx)]
    return {
        "input_file_sets": [
            refs.trained_model_set_alias(ctx),  # the model Set, not the FILE derived_from uses
            refs.principal_pseudobulk_set_alias(ctx),  # itself contains the RNA/ATAC files, "and more"
        ],
        "derived_from": ",".join(derived_from_parts),
        "documents": [make_alias(ctx.igvf, "E2G_prediction_set_files")],
        "description": f"scE2G {family(ctx.model)} predictions for {ctx.dataset} {ctx.cluster} cells",
        "samples": _unique_subsamples(ctx),  # TODO: raw labels vs Sample-object aliases -- see module docstring
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="prediction_set",  # TODO: confirm actual portal profile id
        scope="cluster_model",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab"],
        constant_fields={"file_set_type": "functional effect", "scope": "genome-wide"},
        variants=[
            # name="" (not e.g. "default") deliberately -- every other table's
            # depends_on=[("prediction_set", "")] must match this exactly;
            # state.py normalizes a missing/None variant to "" the same way.
            registry.VariantSpec(
                name="",
                build_row=_row,
                depends_on=lambda ctx: [("principal_pseudobulk_set", "")],
            ),
        ],
    )
)
