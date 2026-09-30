"""Prediction Tabular Files -- object_type "tabular_file", scope
cluster_model. Four content variants per (dataset, cluster, model): full,
thresholded, bedpe, elements_bed -- 8 rows per cluster once both scE2G
Multiome and scATAC are configured for it (bedpe/elements_bed are
required for both models, per the 2026-07-10 design conversation).

2026-09-25: the "genes" variant (content_type "gene quantifications") is
retired -- gene quantifications are no longer a pipeline project. scE2G
predictions already score the full candidate-gene universe, and E2G Scores
for candidate links involving a gene with < 1 TPM are zero-ed out
downstream, so a separate per-cluster gene list file added nothing IGVF
didn't already have via the fixed, externally-curated
jesse-engreitz:sce2g_gencode_v43_gene_list. Every Multiome "full" row's
derived_from now references that fixed alias (refs.gene_list_alias) in
place of the old per-(dataset,cluster) "...predictions_genes" alias -- see
_full_row. The reformat.smk rule that produced
{cluster_dir}/{dataset}_{cluster}_scE2G_multiome_v3_gene_list.tsv.gz
(reformat_gene_list) is removed along with it; that file is no longer a
pipeline output target.

Upload order within this table, driven by depends_on:
    elements_bed  (no dependency on sibling prediction rows)
              |
              v
             full  (derived_from unconditionally lists elements_bed's
                     alias now that it's required content, not
                     optional -- confirmed 2026-07-16 -- so full's
                     depends_on must actually wait on it being uploaded
                     first, not just check enabled(); the module docstring
                     described this ordering from the start but the
                     dependency itself was missing until now)
              |
              v
        thresholded
              |
              v
           bedpe

All four also depend on the "prediction_set" table (needed for file_set) --
see igvf_metadata.refs.prediction_set_alias for why using its (provisional,
unconfirmed) alias formula now is still safe even though that table isn't
registered yet: depends_on still gates real uploads on it, so this only
replaces an "invalid" validation failure with a correctly-labeled
"deferred" outcome.

full/elements_bed also need the ATAC fragment TabularFile / RNA count
matrix MatrixFile tables (for derived_from) -- resolved, routed through
igvf_metadata.refs (refs.atac_fragment_alias / refs.rna_matrix_alias),
which delegate to the real filtered_atac_fragment_file / filtered_rna_count_matrix
TableSpecs (see refs.py's own docstring). refs.gene_list_alias, by contrast,
is a fixed constant (like refs.trained_model_file_alias) -- no depends_on
entry needed, since it names an already-existing, externally-managed Portal
object this pipeline never creates.

bedpe/elements_bed path builders (2026-07-21 update): both now resolve to
real, already-produced pipeline outputs, and declare the same required_paths
that full/thresholded do -- no other change needed.
- elements_bed (renamed from "elements" 2026-08-21, content itself reworked
  2026-08-19): NOT derived from any prediction -- confirmed on the IGVF Data
  Portal, elements are represented as an INPUT to the prediction file (the
  pre-scoring candidate universe), not something derived from it.
  reformat.smk's reformat_element_list + bgzip_index_element_list rules now
  source {cluster_dir}/Neighborhoods/EnhancerList.bed (scE2G's own ABC
  candidate-element output, cluster-level and model-agnostic -- confirmed
  byte-for-byte the same coordinates/class, same order, as the
  multiome_powerlaw_v3-scoped scE2G_element_list.tsv.gz this used to read
  instead) into
  {cluster_dir}/{dataset}_{cluster}_element_list.bed.gz(.tbi) -- a
  bgzip+tabix-indexed, IGV/UCSC-loadable BED4 (ElementChr/Start/End/Name,
  Name = "{class}|{chr}:{start}-{end}", column header row "#"-commented).
  Deliberately NOT named with "multiome_v3" or any model tag, and NOT gated
  on any model being present in get_reformat_output_files() -- EnhancerList.bed
  exists for every reformat-eligible cluster regardless of model, matching
  this table's "elements required for both
  families" design (previously ungenerateable for scATAC-only clusters,
  since scE2G never writes scE2G_element_list.tsv.gz outside
  multiome_powerlaw_v3/). Variant renamed "elements" -> "elements_bed" (and
  its alias's trailing "...predictions_elements" ->
  "...predictions_elements_bed") once file_format became "bed" instead of
  "tsv" -- file_format_specifications aliases renamed to
  elements_reference_bed_file_format_specification_{pdf,md} to match, and
  submitted_assembly (never actually set on this row, so nothing to remove
  in code) is not a field of this row.
- bedpe: scE2G's own write_sc_e2g_predictions_bedpe rule only emits a plain,
  uncompressed, unindexed .bedpe -- IGVF Portal submission requires
  .bedpe.gz + .bedpe.gz.tbi. workflow/rules/qc_stats.smk's new
  bgzip_index_bedpe rule (gated on make_IGV_tracks, same as the raw .bedpe
  itself) produces
  {cluster_dir}/{dataset}_{cluster}_scE2G_{model}_threshold{threshold}.bedpe.gz;
  BEDPE Index File's row (refs.bedpe_prediction_path(ctx) + ".tbi") follows
  automatically.

Family-gating ("only Multiome unless scATAC is configured," 2026-07-20
feedback) needs no code here -- enforced once, centrally, in
orchestrator._iter_scopes via IgvfConfig.enabled_families, shared by every
scope="cluster_model" table.

reference_files (added 2026-08-06): confirmed via a live upload test against
the real portal to resolve an upload error on this object type. Fixed,
content_type-keyed accessions -- not derived per-row -- so every
"element to gene interactions"/"elements reference" row (full, thresholded,
bedpe, elements_bed) gets the same 3.
"""

import glob
import os
import re

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "prediction_tabular_files"

_REFERENCE_FILES_ELEMENT_TO_GENE = ["IGVFFI7969JLFC", "IGVFFI0653VCGH", "IGVFFI9573KOZR"]

# Extend as new models are added -- deliberately a lookup, not a string
# transform, since "multiome_powerlaw_v3" -> "Multiome" but "scATAC_powerlaw_v3"
# -> "scATAC" don't follow one consistent capitalization rule.
FAMILY_DISPLAY = {
    "multiome_powerlaw_v3": "Multiome",
    "multiome_megamap_v3": "Multiome",
    "scATAC_powerlaw_v3": "scATAC",
    "scATAC_megamap_v3": "scATAC",
}


def family(model):
    try:
        return FAMILY_DISPLAY[model]
    except KeyError:
        raise ValueError(f"no display-family mapping for model {model!r} -- add it to FAMILY_DISPLAY")


def build_alias(ctx, variant_name):
    return make_alias(ctx.igvf, ctx.dataset, ctx.cluster, "scE2G", family(ctx.model), "predictions", variant_name)


def _score_threshold(ctx):
    """Memoized per model (not per cluster) -- confirmed constant across every
    run of a given model. Reads models/{model}/score_threshold_<decimal>, an
    empty marker file scE2G writes once per model."""
    cache_key = ("score_threshold", ctx.model)
    if cache_key not in ctx.cache:
        pattern = os.path.join(ctx.scE2G_dir, "models", ctx.model, "score_threshold_*")
        matches = glob.glob(pattern)
        if len(matches) != 1:
            raise ValueError(f"expected exactly one score_threshold_* file matching {pattern}, found {len(matches)}")
        m = re.search(r"score_threshold_([0-9.]+)$", os.path.basename(matches[0]))
        if not m:
            raise ValueError(f"couldn't parse a threshold decimal out of {matches[0]}")
        threshold = m.group(1)
        # The marker file itself omits the leading zero (e.g. "score_threshold_.177"),
        # but the e2g output filename includes it ("..._threshold0.177.e2g.tsv.gz") --
        # confirmed against models/multiome_powerlaw_v3/score_threshold_.177 on disk.
        if threshold.startswith("."):
            threshold = "0" + threshold
        ctx.cache[cache_key] = threshold
    return ctx.cache[cache_key]


def _full_path(ctx):
    return os.path.join(ctx.cluster_dir, f"{ctx.dataset}_{ctx.cluster}_scE2G_{ctx.model}.e2g.tsv.gz")


def _thresholded_path(ctx):
    threshold = _score_threshold(ctx)
    return os.path.join(ctx.cluster_dir, f"{ctx.dataset}_{ctx.cluster}_scE2G_{ctx.model}_threshold{threshold}.e2g.tsv.gz")


def _bedpe_path(ctx):
    threshold = _score_threshold(ctx)
    return os.path.join(
        ctx.cluster_dir, f"{ctx.dataset}_{ctx.cluster}_scE2G_{ctx.model}_threshold{threshold}.bedpe.gz",
    )


def _elements_bed_path(ctx):
    return os.path.join(ctx.cluster_dir, f"{ctx.dataset}_{ctx.cluster}_element_list.bed.gz")


def _requires(path_fn):
    """Declares the one pipeline output this variant describes. Replaces the old
    _existing_file_enabled() gate (2026-08-24): a missing file is now reported by
    the orchestrator as skipped-missing-file WITH the path, instead of collapsing
    into the same silent skipped-disabled bucket as a deliberate semantic
    exclusion (a variant's enabled() gate returning False).

    The old wrapper also carried an `except NotImplementedError: return False`
    clause, which was dead code -- refs.py's stubs were all resolved, so nothing
    in this chain raises it. Worse, it never caught what actually escapes here:
    _score_threshold raises ValueError on 0 or >=2 score_threshold_* markers,
    which is a REAL misconfiguration and must stay loud."""

    def _fn(ctx):
        return [path_fn(ctx)]

    return _fn


def _full_row(ctx):
    # Confirmed 2026-07-16: RNA count matrix + gene list are Multiome-only data
    # points -- a scATAC full's derived_from omits both, even though ATAC
    # fragment file + elements apply to both families. depends_on below
    # mirrors this (only waits on the RNA count matrix being uploaded for
    # Multiome; refs.gene_list_alias is a fixed external constant, so it
    # needs no depends_on entry -- see its own docstring in refs.py).
    parts = [
        refs.trained_model_file_alias(ctx),  # derived_from wants the FILE, not the Set
        refs.atac_fragment_alias(ctx),
        build_alias(ctx, "elements_bed"),
    ]
    if family(ctx.model) == "Multiome":
        parts.append(refs.rna_matrix_alias(ctx))
        parts.append(refs.gene_list_alias(ctx))
    return {
        "content_type": "element to gene interactions",
        "file_format": "tsv",
        # "full" is every candidate link scE2G scored, not a threshold-filtered
        # subset -- unlike thresholded/bedpe below. (elements_bed also omits
        # this comparison: it doesn't set "filtered" at all -- see its own
        # build_row.)
        "filtered": False,
        "description": f"Full scE2G ({family(ctx.model)}) predictions for {ctx.dataset} {ctx.cluster} cells",
        "derived_from": ",".join(parts),
        # scE2G_feature_table_file_format_specification_{pdf,md} added
        # 2026-09-30: "full" also carries the raw scE2G feature-table columns
        # (unlike thresholded/bedpe, which only ever show the standard
        # E2G_prediction_file_format columns), so it needs that second spec
        # alongside element_to_gene_interaction_predictions_file_format.
        # Matches file_format_spec_igvf_manifests/prediction_tabular_files_full_gene_list_patch.tsv,
        # already applied to the 138 already-live "full" records -- this
        # keeps future uploads consistent with those.
        "file_format_specifications": [
            make_alias(ctx.igvf, "element_to_gene_interaction_predictions_file_format_pdf"),
            make_alias(ctx.igvf, "element_to_gene_interaction_predictions_file_format_md"),
            make_alias(ctx.igvf, "scE2G_feature_table_file_format_specification_pdf"),
            make_alias(ctx.igvf, "scE2G_feature_table_file_format_specification_md"),
        ],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _full_path(ctx),
    }


def _thresholded_row(ctx):
    return {
        "content_type": "element to gene interactions",
        "file_format": "tsv",
        # Score >= threshold has been applied to "full" to produce this row
        # (reformat_predictions_thresholded's own --threshold flag) -- unlike
        # full/elements_bed, which carry every candidate.
        "filtered": True,
        "description": f"Thresholded scE2G ({family(ctx.model)}) predictions for {ctx.dataset} {ctx.cluster} cells",
        "derived_from": build_alias(ctx, "full"),
        "file_format_specifications": [
            make_alias(ctx.igvf, "element_to_gene_interaction_predictions_file_format_pdf"),
            make_alias(ctx.igvf, "element_to_gene_interaction_predictions_file_format_md"),
        ],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _thresholded_path(ctx),
    }


def _bedpe_row(ctx):
    return {
        "content_type": "element to gene interactions",
        "file_format": "bedpe",
        # Same filtered subset as "thresholded" (derived_from it below), just
        # repackaged as bedpe for genome-browser visualization.
        "filtered": True,
        "description": (
            f"Bedpe file for genome browser visualization of thresholded scE2G ({family(ctx.model)}) "
            f"predictions for {ctx.dataset} {ctx.cluster} cells"
        ),
        "derived_from": build_alias(ctx, "thresholded"),
        "file_format_specifications": [make_alias(ctx.igvf, "E2G_bedpe_file_format")],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _bedpe_path(ctx),
    }


def _elements_bed_row(ctx):
    return {
        "content_type": "elements reference",
        "file_format": "bed",
        # REQUIRED for file_format bed/bigBed -- tabular_file.json says "Bed files
        # must specify their file format type", and omitting it is not a soft
        # warning: the whole row fails validation. The failure message is
        # misleading, because the constraint is an anyOf and jsonschema reports the
        # OTHER branch, so you get "'bed' is not one of ['bedpe', 'csv', 'gtf',
        # 'gvcf', 'tar', 'tsv', 'vcf']" -- which reads as "bed is not allowed at
        # all" when bed IS allowed and the real problem is the missing companion
        # field. Caught by --mode validate on 2026-08-26; this row had never
        # successfully uploaded.
        #
        # bed3+ (not bed4, which the enum does not offer) matches
        # filtered_atac_fragment_file.py's established choice for the same shape:
        # three standard BED columns plus extras. This file is
        # chr/start/end/ElementName.
        "file_format_type": "bed3+",
        # No "filtered" key at all here (2026-09-30) -- not True, not False,
        # entirely absent. Elements are the pre-scoring candidate universe
        # (see module docstring), so "filtered" doesn't cleanly apply either
        # way; removed live from every already-uploaded elements_bed record
        # via igvf_utils.connection.Connection.remove_props (see
        # workflow/scripts/remove_filtered_from_elements_bed.py) rather than
        # set to False, so new uploads must omit the key too or a future
        # reconciliation PATCH would silently reintroduce it.
        "description": (
            f"Annotated elements in scE2G ({family(ctx.model)}) predictions for {ctx.dataset} {ctx.cluster} "
            "cells for IGV visualization"
        ),
        "derived_from": refs.atac_fragment_alias(ctx),
        "file_format_specifications": [
            make_alias(ctx.igvf, "elements_reference_bed_file_format_specification_pdf"),
            make_alias(ctx.igvf, "elements_reference_bed_file_format_specification_md"),
        ],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _elements_bed_path(ctx),
    }


def _scope_fields(ctx):
    return {
        "md5sum": None,  # left blank; igvf_utils computes+fills it from submitted_file_name
        "file_set": refs.prediction_set_alias(ctx),  # provisional formula -- see refs.py
    }


def _dep(*extra):
    """Every variant depends on this scope's Prediction Set row (for
    file_set), plus whatever intra-table ordering the variant itself needs."""
    return lambda ctx: [("prediction_set", "")] + list(extra)


def _full_depends_on(ctx):
    deps = [("prediction_set", ""), ("filtered_atac_fragment_file", ""), (TABLE_NAME, "elements_bed")]
    if family(ctx.model) == "Multiome":
        deps.append(("filtered_rna_count_matrix", ""))
    return deps


TABLE = registry.register(
    registry.TableSpec(
        name=TABLE_NAME,
        object_type="tabular_file",
        scope="cluster_model",
        build_alias=build_alias,
        required_columns=["aliases", "award", "lab", "content_type", "controlled_access", "file_format", "file_set"],
        constant_fields={
            "controlled_access": False,
            # filtered is NOT here: thresholded/bedpe set it True, full sets
            # it False, and elements_bed doesn't set it at all (2026-09-30 --
            # not True, not False, entirely absent) -- no single default, or
            # even a single True/False split, applies to all four.
            "derived_manually": False,
            "analysis_step_version": "jesse-engreitz:analysis_step_v1_run_scE2G",
        },
        scope_fields=_scope_fields,
        variants=[
            registry.VariantSpec(
                name="full",
                build_row=_full_row,
                required_paths=_requires(_full_path),
                depends_on=_full_depends_on,
            ),
            registry.VariantSpec(
                name="thresholded",
                build_row=_thresholded_row,
                required_paths=_requires(_thresholded_path),
                depends_on=_dep((TABLE_NAME, "full")),
            ),
            registry.VariantSpec(
                name="bedpe",
                build_row=_bedpe_row,
                required_paths=_requires(_bedpe_path),
                depends_on=_dep((TABLE_NAME, "thresholded")),
            ),
            registry.VariantSpec(
                name="elements_bed",
                build_row=_elements_bed_row,
                required_paths=_requires(_elements_bed_path),
                depends_on=_dep(("filtered_atac_fragment_file", "")),
            ),
        ],
    )
)
