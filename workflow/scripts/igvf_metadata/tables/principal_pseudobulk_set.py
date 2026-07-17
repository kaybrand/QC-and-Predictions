"""Principal Pseudobulk Set -- object_type "pseudobulk_set" (TODO: confirm
portal profile id; multireport query used type=PseudobulkSet). Scope
"cluster" -- one per (dataset, cluster), NOT split by model: up to 10
models' Prediction Sets can point at the same Principal Pseudobulk Set via
refs.principal_pseudobulk_set_alias (see prediction_set.py's
input_file_sets).

input_file_sets here points to primary pseudobulk sets contributed by a
DIFFERENT lab (Anshul Kundaje's), in THEIR alias namespace -- confirmed
2026-07-13: "anshul-kundaje:{dataset}-{cluster}-{subsample}" (hyphens, not
this lab's usual underscore/jesse-engreitz convention). One entry per
relevant subsample. Hardcoded as a literal, not derived from
ctx.igvf.alias_prefix, since it names someone else's namespace, not ours.

cell_type/cell_qualifier still need the live-portal-lookup mechanism
discussed but not yet resolved (which field to index the multireport
response on) -- see refs.primary_pseudobulk_metadata, a stub.

samples: the given spec says "extract from primary pseudobulk sets or
filtered QC TabularFiles" -- read as two alternative sources, not two
combined ones. For now reuses the same local subsamples.unique_subsamples
list that both this field and input_file_sets need anyway (and that
Prediction Set's own `samples` field already uses) rather than round-
tripping through the live portal query for the same information -- flag if
that's wrong and it should instead come from the primary pseudobulk sets'
own `samples` field via the multireport lookup.

documents references a "QC_thresholds" Document per (dataset, cluster) --
not yet a registered table (no "Documents" module built), so, like
"prediction_set" before it existed, depends_on lists it so real uploads
correctly wait once it is.
"""

import os

from .. import refs, registry, subsamples
from ..context import make_alias

TABLE_NAME = "principal_pseudobulk_set"

_KUNDAJE_ALIAS_PREFIX = "anshul-kundaje"  # a different lab's namespace, not ours -- see module docstring


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "filtered_pseudobulk_set")


def _primary_pseudobulk_aliases(ctx):
    return [f"{_KUNDAJE_ALIAS_PREFIX}:{ctx.dataset}-{ctx.cluster}-{s}" for s in subsamples.unique_subsamples(ctx)]


def _row(ctx):
    metadata = refs.primary_pseudobulk_metadata(ctx)  # raises until the indexing key is resolved
    return {
        "cell_type": f"/sample-terms/{metadata['cl_id']}/",
        "cell_qualifier": metadata["cell_qualifier"],
        "documents": [refs.qc_thresholds_document_alias(ctx)],  # the QC_thresholds Document (not yet built)
        "samples": subsamples.unique_subsamples(ctx),
        "input_file_sets": _primary_pseudobulk_aliases(ctx),
        "description": (
            f"Filtered datafiles describing a single annotated cell cluster ({ctx.cluster}); "
            "these datafiles are inputs to E2G predictive models"
        ),
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="pseudobulk_set",  # TODO: confirm actual portal profile id
        scope="cluster",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "file_set_type"],
        constant_fields={"file_set_type": "pseudobulk analysis", "merged": True, "controlled_access": False},
        variants=[
            # name="" to match every dependent table's depends_on=[("principal_pseudobulk_set", "")]
            registry.VariantSpec(
                name="",
                build_row=_row,
                depends_on=lambda ctx: [("documents", "")],
            ),
        ],
    )
)
