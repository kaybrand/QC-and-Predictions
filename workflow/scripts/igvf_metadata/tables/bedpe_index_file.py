"""BEDPE Index File -- object_type "tabular_file" (TODO: confirm -- mirrors
Prediction Tabular Files' shape closely enough that this is a reasonable
guess, but could be a distinct index/reference-file profile). Scope
cluster_model, one row per (dataset, cluster, model): the .tbi index for
the thresholded BEDPE prediction file.

Confirmed 2026-07-13: aliases is
"jesse-engreitz:{dataset}_{cluster}_scE2G_{Family}_predictions_bedpe_index"
(the initially-given template named an unrelated ATAC fragment file and
also collided with the BEDPE prediction file's own alias -- corrected).

submitted_file_name wasn't given -- inferred as the bedpe prediction
file's own (still-unresolved) path + ".tbi" via refs.bedpe_prediction_path,
so this table unblocks automatically alongside Prediction Tabular Files'
bedpe variant once that naming convention lands.
"""

import os

from .. import refs, registry
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "bedpe_index_file"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "predictions", "bedpe", "index")


def _path(ctx):
    return refs.bedpe_prediction_path(ctx) + ".tbi"


def _enabled(ctx):
    try:
        return os.path.exists(_path(ctx))
    except NotImplementedError:
        return False


def _row(ctx):
    prediction_table = registry.get("prediction_tabular_files")
    return {
        "file_format": "tbi",
        "content_type": "index",
        "reference_files": ["IGVFDS0280IQAI"],
        "derived_from": prediction_table.build_alias(ctx, "bedpe"),
        "description": f"index file for thresholded scE2G {family(ctx.model)} predictions BEDPE file",
        "submitted_file_name": _path(ctx),
    }


def _scope_fields(ctx):
    return {
        "md5sum": None,  # left blank; igvf_utils computes+fills it from submitted_file_name
        "file_set": refs.prediction_set_alias(ctx),
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="tabular_file",
        scope="cluster_model",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "file_format", "file_set", "content_type", "controlled_access"],
        constant_fields={
            "controlled_access": False,
            "derived_manually": False,
            "analysis_step_version": "jesse-engreitz:analysis_step_v1_run_scE2G",
        },
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="bedpe_index",
                build_row=_row,
                enabled=_enabled,
                depends_on=lambda ctx: [("prediction_set", ""), ("prediction_tabular_files", "bedpe")],
            ),
        ],
    )
)
