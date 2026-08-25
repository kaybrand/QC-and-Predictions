"""Prediction Tabular Files -- object_type "tabular_file", scope
cluster_model. Five content variants per (dataset, cluster, model): full,
thresholded, bedpe, elements, genes -- 10 rows per cluster once both scE2G
Multiome and scATAC are configured for it (bedpe/elements/genes are
required for both models, per the 2026-07-10 design conversation).

Upload order within this table, driven by depends_on:
    elements, genes  (no dependency on sibling prediction rows)
              |
              v
             full  (derived_from unconditionally lists elements+genes'
                     aliases now that they're required content, not
                     optional -- confirmed 2026-07-16 -- so full's
                     depends_on must actually wait on them being uploaded
                     first, not just check enabled(); the module docstring
                     described this ordering from the start but the
                     dependency itself was missing until now)
              |
              v
        thresholded
              |
              v
           bedpe

All five also depend on the "prediction_set" table (needed for file_set) --
see igvf_metadata.refs.prediction_set_alias for why using its (provisional,
unconfirmed) alias formula now is still safe even though that table isn't
registered yet: depends_on still gates real uploads on it, so this only
replaces an "invalid" validation failure with a correctly-labeled
"deferred" outcome.

full/elements/genes also need the ATAC fragment TabularFile / RNA count
matrix MatrixFile tables (for derived_from) -- resolved, routed through
igvf_metadata.refs (refs.atac_fragment_alias / refs.rna_matrix_alias),
which delegate to the real filtered_atac_fragment_file / filtered_rna_count_matrix
TableSpecs (see refs.py's own docstring).

bedpe/elements/genes path builders (2026-07-21 update): all three now resolve
to real, already-produced pipeline outputs, and declare the same required_paths
that full/thresholded do -- no other change needed.
- elements/genes: reformat.smk's own reformat_lists rule (already in this
  pipeline) reformats the raw {cluster_dir}/{model}/scE2G_{element,gene}_list.tsv.gz
  scE2G output into the consortium-standard
  {cluster_dir}/{dataset}_{cluster}_scE2G_multiome_v3_{element,gene}_list.tsv.gz
  -- there was no missing reformatting script after all, just a path builder
  that hadn't been pointed at it yet. Same as scE2G's raw list output, this
  is produced once per cluster regardless of which model triggered the row
  (only ever from multiome_powerlaw_v3), matching this table's "elements
  required for both families, genes Multiome-only" design.
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
bedpe, elements) gets the same 3, and every "gene quantifications" row
(genes) gets the 1 that applies to it.
"""

import glob
import os
import re

from .. import refs, registry
from ..context import make_alias

TABLE_NAME = "prediction_tabular_files"

_REFERENCE_FILES_ELEMENT_TO_GENE = ["IGVFFI7969JLFC", "IGVFFI0653VCGH", "IGVFFI9573KOZR"]
_REFERENCE_FILES_GENE_QUANTIFICATIONS = ["IGVFFI9573KOZR"]

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


def _elements_path(ctx):
    return os.path.join(ctx.cluster_dir, f"{ctx.dataset}_{ctx.cluster}_scE2G_multiome_v3_element_list.tsv.gz")


def _genes_path(ctx):
    return os.path.join(ctx.cluster_dir, f"{ctx.dataset}_{ctx.cluster}_scE2G_multiome_v3_gene_list.tsv.gz")


def _requires(path_fn):
    """Declares the one pipeline output this variant describes. Replaces the old
    _existing_file_enabled() gate (2026-08-24): a missing file is now reported by
    the orchestrator as skipped-missing-file WITH the path, instead of collapsing
    into the same silent skipped-disabled bucket as a deliberate semantic
    exclusion like "genes under scATAC".

    The old wrapper also carried an `except NotImplementedError: return False`
    clause, which was dead code -- refs.py's stubs were all resolved, so nothing
    in this chain raises it. Worse, it never caught what actually escapes here:
    _score_threshold raises ValueError on 0 or >=2 score_threshold_* markers,
    which is a REAL misconfiguration and must stay loud."""

    def _fn(ctx):
        return [path_fn(ctx)]

    return _fn


def _full_row(ctx):
    # Confirmed 2026-07-16: RNA count matrix + genes are Multiome-only data
    # points -- a scATAC full's derived_from omits both, even though ATAC
    # fragment file + elements apply to both families. depends_on below
    # mirrors this (only waits on "genes" being uploaded for Multiome).
    parts = [
        refs.trained_model_file_alias(ctx),  # derived_from wants the FILE, not the Set
        refs.atac_fragment_alias(ctx),
        build_alias(ctx, "elements"),
    ]
    if family(ctx.model) == "Multiome":
        parts.append(refs.rna_matrix_alias(ctx))
        parts.append(build_alias(ctx, "genes"))
    return {
        "content_type": "element to gene interactions",
        "file_format": "tsv",
        "description": f"Full scE2G ({family(ctx.model)}) predictions for {ctx.dataset} {ctx.cluster} cells",
        "derived_from": ",".join(parts),
        "file_format_specifications": [
            make_alias(ctx.igvf, "element_to_gene_interaction_predictions_file_format_pdf"),
            make_alias(ctx.igvf, "element_to_gene_interaction_predictions_file_format_md"),
        ],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _full_path(ctx),
    }


def _thresholded_row(ctx):
    return {
        "content_type": "element to gene interactions",
        "file_format": "tsv",
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
        "description": (
            f"Bedpe file for genome browser visualization of thresholded scE2G ({family(ctx.model)}) "
            f"predictions for {ctx.dataset} {ctx.cluster} cells"
        ),
        "derived_from": build_alias(ctx, "thresholded"),
        "file_format_specifications": [make_alias(ctx.igvf, "E2G_bedpe_file_format")],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _bedpe_path(ctx),
    }


def _elements_row(ctx):
    return {
        "content_type": "elements reference",
        "file_format": "tsv",
        "description": f"Annotated elements in scE2G ({family(ctx.model)}) predictions for {ctx.dataset} {ctx.cluster} cells",
        "derived_from": refs.atac_fragment_alias(ctx),
        "file_format_specifications": [
            make_alias(ctx.igvf, "elements_reference_file_format_specification_pdf"),
            make_alias(ctx.igvf, "elements_reference_file_format_specification_md"),
        ],
        "reference_files": _REFERENCE_FILES_ELEMENT_TO_GENE,
        "submitted_file_name": _elements_path(ctx),
    }


def _genes_row(ctx):
    return {
        "content_type": "gene quantifications",
        "file_format": "tsv",
        "description": f"Annotated genes in scE2G ({family(ctx.model)}) predictions for {ctx.dataset} {ctx.cluster} cells",
        "derived_from": refs.rna_matrix_alias(ctx),
        "file_format_specifications": [
            make_alias(ctx.igvf, "gene_quantifications_file_format_specification_pdf"),
            make_alias(ctx.igvf, "gene_quantifications_file_format_specification_md"),
        ],
        "reference_files": _REFERENCE_FILES_GENE_QUANTIFICATIONS,
        "submitted_file_name": _genes_path(ctx),
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


def _genes_enabled(ctx):
    # RNA count matrix + genes are Multiome-only data points (2026-07-16) --
    # a real semantic gate, so a stray file at the expected path could never
    # produce a "genes" row under scATAC. File existence is declared separately
    # via required_paths, so "no gene list on disk" and "scATAC, so no gene list
    # should exist" are reported as different outcomes rather than one bucket.
    return family(ctx.model) == "Multiome"


def _full_depends_on(ctx):
    deps = [("prediction_set", ""), ("filtered_atac_fragment_file", ""), (TABLE_NAME, "elements")]
    if family(ctx.model) == "Multiome":
        deps.append((TABLE_NAME, "genes"))
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
            "filtered": True,
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
                name="elements",
                build_row=_elements_row,
                required_paths=_requires(_elements_path),
                depends_on=_dep(("filtered_atac_fragment_file", "")),
            ),
            registry.VariantSpec(
                name="genes",
                build_row=_genes_row,
                enabled=_genes_enabled,
                required_paths=_requires(_genes_path),
                depends_on=_dep(("filtered_rna_count_matrix", "")),
            ),
        ],
    )
)
