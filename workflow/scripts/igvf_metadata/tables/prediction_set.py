"""Prediction Set -- the FileSet object binding together every scE2G model
output for one (dataset, cluster, model) on the IGVF Data Portal. Every
other table's file_set points at this table's alias; see
igvf_metadata.refs.prediction_set_alias, which delegates here.

input_file_sets is fully resolved: exactly [trained_model_set_alias,
principal_pseudobulk_set_alias]. _row() itself no longer raises -- but
depends_on lists ("principal_pseudobulk_set", "") so real uploads still
correctly wait for that object to actually exist on the portal before this
one links to it.

No derived_from field: the Prediction Set profile doesn't have one at all
(confirmed 2026-08-05) -- it was previously duplicated here as a
pre-joined-string copy of input_file_sets, which was both invalid (not a
real property on this schema) and redundant even if it had been.

Resolved 2026-07-20: samples' raw "subsample" column values (e.g.
"IGVFSM6456LUAO") ARE literal portal Sample accessions already -- on the
IGVF Data Portal, "Sample" is the In-Vitro-System object type, distinct from
"biosample" (this pipeline's "cluster") and from this pipeline's own
"subsample"/in-vitro-system vocabulary. Verifiable directly from the
accession itself: every IGVF accession is IGVF + a 2-character object-type
code + an 8-character alphanumeric ID, and "SM" is the Sample type code --
"IGVFSM6456LUAO" decomposes as IGVF/SM/6456LUAO. No further alias-resolution
needed; passing the raw values straight through is correct. Ordered by
descending barcode-count frequency (2026-07-20 feedback: most-contributing
subsample first), via subsamples.subsamples_by_frequency.

Family-gating ("only Multiome unless scATAC is configured," 2026-07-20
feedback) needs no code here -- enforced once, centrally, in
orchestrator._iter_scopes via IgvfConfig.enabled_families, shared by every
scope="cluster_model" table.
"""

from .. import refs, registry, subsamples
from ..context import make_alias
from .prediction_tabular_files import family

TABLE_NAME = "prediction_set"


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "predictions")


def _row(ctx):
    return {
        "input_file_sets": [
            refs.trained_model_set_alias(ctx),  # the model Set, not the FILE derived_from-shaped fields use
            refs.principal_pseudobulk_set_alias(ctx),  # itself contains the RNA/ATAC files, "and more"
        ],
        "documents": [
            make_alias(ctx.igvf, "E2G_prediction_set_file_format_specs_pdf"),
            make_alias(ctx.igvf, "E2G_prediction_set_file_format_specs_md"),
        ],
        "description": f"scE2G {family(ctx.model)} predictions for {ctx.dataset} {ctx.cluster} cells",
        "samples": subsamples.subsamples_by_frequency(ctx),  # raw values ARE Sample accessions already -- see module docstring
    }


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="prediction_set",  # confirmed 2026-08-03
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
