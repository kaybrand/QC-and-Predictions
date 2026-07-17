"""Signal Files -- object_type "signal_file" (TODO: confirm this is the
actual portal profile id used for --profile_id), scope cluster_model.
Exactly one row per (dataset, cluster, model): the ATAC read-depth bigWig.

No model gate here (or anywhere else in this package): which models get
uploaded is entirely a function of each cluster's `models` list in the
pipeline config (multiome_powerlaw_v3-only in 2026, a mixture of both once
scATAC is added in 2027) -- the same iteration this table shares with every
other cluster_model-scoped table already handles that correctly with zero
code changes when a cluster's config grows a second model.

file_set reuses refs.prediction_set_alias -- confirmed formula (2026-07-13).

derived_from (the filtered ATAC fragment TabularFile) is a true stub via
refs.atac_fragment_alias, same as Prediction Tabular Files' full/elements
rows -- fill in refs.py once that table's alias format is known.

submitted_file_name does NOT vary by model in the given path
(.../{dataset}/{cluster}/ATAC_norm.bw) even though this table is
model-scoped -- fine while only one model is configured per cluster, but
will collide if a cluster ever has two models' signal files landing at the
same path. Flagging for when scATAC signal files are added.
"""

import os

from .. import refs, registry
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "signal_files"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "ATAC_bw")


def _path(ctx):
    return os.path.join(ctx.cluster_dir, "ATAC_norm.bw")


def _enabled(ctx):
    return os.path.exists(_path(ctx))


def _row(ctx):
    return {
        "file_format": "bigWig",
        "content_type": "read-depth signal",
        "reference_files": ["IGVFDS0280IQAI"],
        "strand_specificity": "unstranded",
        "derived_from": refs.atac_fragment_alias(ctx),
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
        object_type="signal_file",
        scope="cluster_model",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "file_format", "file_set", "content_type"],
        constant_fields={
            "assembly": "GRCh38",
            "normalized": True,
            "analysis_step_version": "jesse-engreitz:analysis_step_v1_run_scE2G",
        },
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="atac_bw",
                build_row=_row,
                enabled=_enabled,
                depends_on=lambda ctx: [("prediction_set", ""), ("filtered_atac_fragment_file", "")],
            ),
        ],
    )
)
