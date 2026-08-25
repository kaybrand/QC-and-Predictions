"""BEDPE Index File -- object_type "index_file" (confirmed 2026-08-03, same
as ATAC Index File -- not tabular_file as previously guessed). Scope
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

reference_files removed (2026-07-20 feedback): not a submittable field for
type IndexFile (content_type=index).

Family-gating (2026-07-20 feedback: only Family=multiome rows, scATAC
infrastructure stays in place but isn't triggered) needs no code here --
it's enforced once, centrally, in orchestrator._iter_scopes via
IgvfConfig.enabled_families, shared by every scope="cluster_model" table.
"""


from .. import refs, registry
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "bedpe_index_file"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "predictions", "bedpe", "index")


def _path(ctx):
    return refs.bedpe_prediction_path(ctx) + ".tbi"


def _required_paths(ctx):
    return [_path(ctx)]


def _row(ctx):
    prediction_table = registry.get("prediction_tabular_files")
    return {
        "file_format": "tbi",
        "content_type": "index",
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
                name="bedpe_index",
                build_row=_row,
                required_paths=_required_paths,
                depends_on=lambda ctx: [("prediction_set", ""), ("prediction_tabular_files", "bedpe")],
            ),
        ],
    )
)
