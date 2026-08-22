"""Elements BED Index File -- object_type "index_file", scope cluster_model,
one row per (dataset, cluster, model): the .tbi index for Prediction Tabular
Files' elements_bed variant ({dataset}_{cluster}_element_list.bed.gz.tbi).
Same shape as BEDPE Index File (bedpe_index_file.py), including its aliases
formula "jesse-engreitz:{dataset}_{cluster}_scE2G_{Family}_predictions_elements_bed_index".

submitted_file_name isn't given directly -- inferred as the elements BED
file's own (still-unresolved) path + ".tbi" via refs.elements_bed_path, same
"unblocks automatically alongside Prediction Tabular Files' variant" reasoning
as BEDPE Index File.

reference_files omitted, same as BEDPE/ATAC Index File (2026-07-20 feedback):
not a submittable field for type IndexFile (content_type=index).

Family-gating needs no code here, same as every other scope="cluster_model"
table (enforced once, centrally, in orchestrator._iter_scopes via
IgvfConfig.enabled_families) -- unlike BEDPE Index File (Multiome-only in
practice today), elements_bed itself is required for both families per
Prediction Tabular Files' own docstring, so this table's rows follow
whichever families are enabled.
"""

import os

from .. import refs, registry
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "elements_bed_index_file"


def build_alias(ctx, variant_name):
    return make_alias(
        ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "predictions", "elements_bed", "index",
    )


def _path(ctx):
    return refs.elements_bed_path(ctx) + ".tbi"


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
        "derived_from": prediction_table.build_alias(ctx, "elements_bed"),
        "description": (
            f"index file for annotated candidate elements in scE2G ({family(ctx.model)}) "
            f"predictions for {ctx.dataset} {ctx.cluster} cells"
        ),
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
        object_type="index_file",
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
                name="elements_bed_index",
                build_row=_row,
                enabled=_enabled,
                depends_on=lambda ctx: [("prediction_set", ""), ("prediction_tabular_files", "elements_bed")],
            ),
        ],
    )
)
