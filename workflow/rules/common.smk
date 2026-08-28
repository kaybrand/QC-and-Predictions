"""
Shared config loading, exclusion resolution, and helper functions used across
all rule files in this pipeline. Included once from the top-level
workflow/Snakefile before any other rule file.

Supports clusters from multiple datasets in a single run (see "Multi-dataset
support" in the design plan): `config["clusters"]` is nested by dataset, and
every internal structure here is keyed by the (dataset, cluster) tuple, not
bare cluster name -- cluster names are only guaranteed unique within their
own dataset, not globally.
"""

import glob
import os
import re
import sys

import yaml

# Repo root (the directory containing this repo's own Snakefile/workflow/...),
# derived from the location of the top-level workflow/Snakefile rather than the
# current working directory. Code directory only -- for the pipeline's own
# reference/ files etc. NEVER for datafiles (plots/, datatables/,
# multiome_data/): the code checkout and the data location are independent,
# since a worktree checks the same code out at a different path than wherever
# the data actually lives. Datafile roots come from config["data_dir"] instead
# (see multiome_data_dir/resolve_exclusions below) -- a self-consistent path
# the user sets once, in the pipeline config.
WDIR = os.path.dirname(workflow.basedir)

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from resolve_exclusions import (  # noqa: E402
    resolve_exclusions,
    write_cluster_stats_table,
    manifest_eligible_clusters,
    predictions_on_everything,
)
from write_scE2G_config import (  # noqa: E402
    write_cell_clusters_table,
    write_cluster_metadata_table,
    load_cluster_metadata,
)
from pipeline_paths import resolve_repo_relative  # noqa: E402
from cell_annotations import annotation_lookup_key  # noqa: E402

# Read-only snapshot of state.db's cell_annotations, written by the driver
# before Snakemake starts. Replaces the former `from igvf_metadata import state`
# parse-time SQLite connection -- every Slurm worker re-parses this file, and
# multi-host WAL access to a Lustre-hosted SQLite DB is not supported. Nothing
# under workflow/rules/ opens state.db any more.
import cell_annotation_snapshot as cas  # noqa: E402

# Consolidated OUTPUT root for everything this pipeline WRITES (filtered
# ATAC/RNA outputs, scE2G's own results tree, Synapse manifests, cluster-stats
# tables, the coverage report). config["data_dir"] (above) is now READ-ONLY --
# used only to locate pre-existing QC-guide plots/datatables inputs, never
# written to. Relative paths resolve against WDIR (this repo's own root), not
# data_dir and not the current working directory.
OUTPUT_DIR = resolve_repo_relative(config.get("output_dir", "./results"), WDIR)

MAX_MEM_MB = 250 * 1000  # 250GB, matches scE2G_options.max_memory_allocation_mb default


def determine_mem_mb(wildcards, input, attempt, min_gb=8):
    """ABC's memory resource calculator, reused verbatim (see IGVF/workflow/rules/utils.smk)."""
    input_size_mb = input.size_mb
    if ".gz" in str(input):
        input_size_mb *= 8  # assume gz compressed the file <= 8x
    attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
    mem_to_use_mb = attempt_multiplier * max(2 * input_size_mb, min_gb * 1000)
    return min(mem_to_use_mb, MAX_MEM_MB)


def resolve_primary_model(models):
    """Exactly one candidates file and one feature table per cluster: prefer
    multiome_powerlaw_v3 when both models were run, else scATAC_powerlaw_v3."""
    return "multiome_powerlaw_v3" if "multiome_powerlaw_v3" in models else "scATAC_powerlaw_v3"


def get_model_threshold(scE2G_dir, model_name):
    """Read a model's score threshold from its models/{model_name}/score_threshold_*
    directory name, instead of hardcoding thresholds that change when a model updates."""
    matches = glob.glob(os.path.join(scE2G_dir, "models", model_name, "score_threshold_*"))
    if not matches:
        raise ValueError(f"No score_threshold_* directory found under models/{model_name}")
    m = re.search(r"score_threshold_(\.?\d+\.?\d*)", os.path.basename(matches[0]))
    if not m:
        raise ValueError(f"Could not parse threshold from {matches[0]}")
    return float(m.group(1))


def _make_paths_absolute(obj, base_path):
    """Recursively resolve relative paths in a config dict against base_path.
    Reused pattern from IGVF/workflow/rules/utils.smk: scE2G's own default
    config.yaml has paths (encode_re2g_dir, model_dir, gene_annotations, ...)
    that are relative to scE2G_dir, meant to be resolved when scE2G is run
    from its own directory -- but importing it as a Snakemake `module` from
    this pipeline does NOT rebase those paths automatically, so they must be
    absolutized here or every scE2G-internal rule that reads them breaks."""
    if isinstance(obj, dict):
        return {k: _make_paths_absolute(v, base_path) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_make_paths_absolute(v, base_path) for v in obj]
    if isinstance(obj, str):
        candidate = os.path.join(base_path, obj)
        if os.path.exists(candidate):
            return candidate
    return obj


def build_scE2G_config(config, scE2G_dir, cell_clusters_table, results_dir):
    """
    Load scE2G's own config/config.yaml as the default, then overlay every
    key set under the master config's `scE2G_options` -- every OPTIONS /
    REFERENCE FILES / WORKFLOW PARAMS field scE2G supports is reachable this
    way, not just a fixed subset (generalizes IGVF/workflow's get_scE2G_config()).
    Called once per dataset -- cell_clusters_table and results_dir are that
    dataset's own partitioned table/output directory.
    """
    with open(os.path.join(scE2G_dir, "config", "config.yaml")) as f:
        scE2G_config = yaml.safe_load(f)

    scE2G_config = _make_paths_absolute(scE2G_config, scE2G_dir)

    scE2G_config["cell_clusters"] = cell_clusters_table
    scE2G_config["results_dir"] = results_dir
    scE2G_config["IGV_dir"] = results_dir

    scE2G_config.update(_make_paths_absolute(config.get("scE2G_options", {}), scE2G_dir))
    return scE2G_config


def multiome_data_dir(dataset):
    return os.path.join(OUTPUT_DIR, "multiome_data", dataset)


def pseudobulks_dir(dataset):
    return os.path.join(config["pseudobulks_root"], dataset, "pseudobulks")


# ---------------------------------------------------------------------------
# Exclusion resolution -- runs once at parse time, before any rule is defined.
# Every set below is keyed by (dataset, cluster) tuples.
# ---------------------------------------------------------------------------
INCLUDED_CLUSTERS, UPLOAD_ELIGIBLE_CLUSTERS, EXCLUDED_CLUSTERS, CLUSTER_STATS = resolve_exclusions(config, config["data_dir"])

if EXCLUDED_CLUSTERS:
    print(f"[QC-and-Predictions] Excluding clusters this run: {sorted(EXCLUDED_CLUSTERS)}")
if not predictions_on_everything(config.get("exclusion", {})):
    print(f"[QC-and-Predictions] Processing clusters: {sorted(INCLUDED_CLUSTERS)}")
else:
    print(f"[QC-and-Predictions] Processing (incl. excluded, no-upload) clusters: {sorted(INCLUDED_CLUSTERS)}")
    print(f"[QC-and-Predictions] Upload-eligible clusters: {sorted(UPLOAD_ELIGIBLE_CLUSTERS)}")

# Upload-eligible MINUS any cluster flagged `igvf_manifest_excluded: true` (the 4
# ATAC-only variant clusters). Enforced here and in the driver that builds
# manage_igvf_metadata.py's --cluster-keys, from the one shared helper -- until
# now the flag was written into configs and reported by generate_report.py but
# never actually kept anything out of a manifest. Deliberately NOT folded into
# REFORMAT_ELIGIBLE_CLUSTERS below: those clusters are still reformatted, just
# never submitted.
MANIFEST_ELIGIBLE_CLUSTERS = manifest_eligible_clusters(config, UPLOAD_ELIGIBLE_CLUSTERS)
if MANIFEST_ELIGIBLE_CLUSTERS != UPLOAD_ELIGIBLE_CLUSTERS:
    print(
        "[QC-and-Predictions] Excluded from IGVF manifests by igvf_manifest_excluded: "
        f"{sorted(UPLOAD_ELIGIBLE_CLUSTERS - MANIFEST_ELIGIBLE_CLUSTERS)}"
    )

# Cluster-stats persistence -- CLUSTER_STATS above is otherwise computed and
# discarded every parse. Write it out at PARSE time (including on `-n` dry
# runs) so the coverage report (generate_report.py) can consume real
# quality-gating data for every configured cluster without recomputing it or
# requiring real execution. Every configured dataset gets a table, not just
# ones with passing clusters -- a dataset whose every cluster failed quality
# still needs its stats reported.
CLUSTER_STATS_DIR = os.path.join(OUTPUT_DIR, "cluster_stats")
for _stats_dataset in config["clusters"]:
    write_cluster_stats_table(_stats_dataset, CLUSTER_STATS, CLUSTER_STATS_DIR)

# ---------------------------------------------------------------------------
# Pipeline mode. Top-level (not nested under a section) specifically so it can be
# overridden per-invocation with `--config pipeline_mode=local_only`.
#
#   local_only -- everything Phase 1 + Phase 2 produce: filtered fragments/RNA
#                 matrices, scE2G predictions/candidates/features, QC report, IGV
#                 tracks, EnhancerList + RNA matrix packaging. NEVER opens
#                 state.db, never reads the CellAnnotation snapshot, never
#                 requests a reformat target. Not "reads the cache and finds
#                 nothing" -- genuinely no contact, so a cold/absent/locked cache
#                 cannot affect it at all. Formalises what used to be an
#                 `--omit-from reformat_predictions ...` CLI convention, which
#                 was fragile: --omit-from takes nargs='+' and silently swallowed
#                 a following target path (a real 5.5-hour incident).
#   default    -- the above PLUS portal-format reformatting, for clusters the
#                 driver's warm stage resolved a CellAnnotation for.
# ---------------------------------------------------------------------------
PIPELINE_MODE = config.get("pipeline_mode", "default")
if PIPELINE_MODE not in ("default", "local_only"):
    raise ValueError(f"pipeline_mode must be 'default' or 'local_only', got {PIPELINE_MODE!r}")

# CellAnnotation-based reformat eligibility: a quality-passing cluster missing
# CellAnnotation still gets full predictions/candidates/features -- it's just not
# reformatted into portal format.
#
# Read from the per-dataset snapshot the driver writes (see
# scripts/cell_annotation_snapshot.py for why this is a file and not the
# state.db query it replaced: every Slurm worker re-parses this file, and
# multi-host WAL access to a Lustre-hosted SQLite DB is unsupported). The
# snapshot is timestamp- and digest-checked on read, so unlike the
# manually-refreshed preview TSV this pipeline used to gate on -- which once
# claimed a cluster was annotated while state.db was cold, crashing a reformat
# rule for real (commit b6e62e0) -- a stale or foreign one raises instead of
# quietly answering the wrong question.
#
# ATAC-only variant clusters resolve via their cell_annotation_key override
# (base name), not their suffixed key.
REFORMAT_ELIGIBLE_CLUSTERS = set()
CELL_ANNOTATIONS = {}
if PIPELINE_MODE == "local_only":
    print(
        "[QC-and-Predictions] local_only mode: portal reformatting and IGVF manifest generation are "
        f"DISABLED by config. {len(UPLOAD_ELIGIBLE_CLUSTERS)} quality-passing cluster(s) still get full "
        "scE2G output. No state.db or CellAnnotation-cache access at all."
    )
else:
    _max_age = config.get("igvf", {}).get("cache_max_age_hours", cas.DEFAULT_MAX_AGE_HOURS)
    _digest = cas.cluster_set_digest(
        {(d, c) for d, clusters in config["clusters"].items() for c in clusters}
    )
    for _dataset in config["clusters"]:
        CELL_ANNOTATIONS.update(
            cas.read_snapshot(
                cas.snapshot_path(OUTPUT_DIR, _dataset),
                expected_digest=_digest,
                max_age_hours=_max_age,
                fix_hint=(
                    "run the pipeline through workflow/scripts/run_pipeline.py (it warms the cache and "
                    "writes this snapshot before invoking Snakemake), or pass "
                    "--config pipeline_mode=local_only to skip everything that needs the portal."
                ),
            )
        )
    for dataset, cluster in UPLOAD_ELIGIBLE_CLUSTERS:
        if annotation_lookup_key(dataset, cluster, config["clusters"][dataset][cluster]) in CELL_ANNOTATIONS:
            REFORMAT_ELIGIBLE_CLUSTERS.add((dataset, cluster))

_missing_annotation = UPLOAD_ELIGIBLE_CLUSTERS - REFORMAT_ELIGIBLE_CLUSTERS
if _missing_annotation and PIPELINE_MODE != "local_only":
    print(
        "[QC-and-Predictions] Quality-passing but no resolvable CellAnnotation (predictions only, no "
        f"reformat -- see the driver's warm-stage status report for the per-cluster reason): "
        f"{sorted(_missing_annotation)}"
    )

# Distinct datasets actually needed this run (a dataset whose clusters were all
# excluded gets no scE2G module instance at all -- no wasted setup).
DATASETS = sorted({dataset for dataset, _ in INCLUDED_CLUSTERS})

# Cluster names contain underscores, and are only unique within their own
# dataset, so filenames like {dataset}_{cluster}_... are genuinely ambiguous
# to split without this: Snakemake would otherwise happily (and wrongly)
# parse "igvf10_telohaec_crispri_..." as dataset=igvf10_telohaec,
# cluster=crispri. Pin both wildcards globally to the known, fixed sets of
# values instead of guessing from the string. The (dataset, cluster) PAIR
# still disambiguates correctly even if the same cluster name recurs in two
# different datasets, since both wildcards are resolved jointly per rule.
ALL_CLUSTER_NAMES = sorted({cluster for _, cluster in INCLUDED_CLUSTERS} | {cluster for _, cluster in EXCLUDED_CLUSTERS})
wildcard_constraints:
    dataset="|".join(re.escape(d) for d in config["clusters"]) if config["clusters"] else "NOMATCH",
    cluster="|".join(re.escape(c) for c in ALL_CLUSTER_NAMES) if ALL_CLUSTER_NAMES else "NOMATCH",
    # Also needed to disambiguate reformat_predictions from reformat_predictions_thresholded:
    # with `model` unconstrained, "multiome_powerlaw_v3_threshold0.177" parses equally well as
    # model="multiome_powerlaw_v3_threshold0.177" (reformat_predictions) or as
    # model="multiome_powerlaw_v3", threshold="0.177" (reformat_predictions_thresholded).
    model="multiome_powerlaw_v3|scATAC_powerlaw_v3",

# ---------------------------------------------------------------------------
# scE2G config table generation -- ALSO happens at parse time (including on
# `snakemake -n` dry runs), not as a lazy rule, and ONCE PER DATASET. This is
# required, not a stylistic choice: scE2G's own Snakefile reads its
# `cell_clusters` TSV at ITS parse time (when each dataset's `module scE2G_*`
# statement is evaluated) to build its biosample expansion list, so each
# dataset's table must already exist on disk before that happens. The
# filtered ATAC/RNA file paths it references don't need to exist yet --
# Snakemake's DAG matches them by path against this pipeline's own
# atac_fragment_file/rna_count_matrix rule outputs and sequences execution
# accordingly.
# ---------------------------------------------------------------------------
CONFIG_DIR = os.path.join(workflow.basedir, "config")
RESULTS_DIR_BASE = os.path.join(OUTPUT_DIR, "uniformly_processed")

CELL_CLUSTERS_TABLES = {}
RESULTS_DIRS = {}
SCE2G_CONFIGS = {}
CLUSTER_METADATA = {}

for _dataset in DATASETS:
    _included_names = {cluster for ds, cluster in INCLUDED_CLUSTERS if ds == _dataset}

    CELL_CLUSTERS_TABLES[_dataset] = write_cell_clusters_table(
        _dataset, config["clusters"][_dataset], _included_names,
        CONFIG_DIR, multiome_data_dir(_dataset), config["scE2G_dir"],
    )
    _metadata_path = write_cluster_metadata_table(
        _dataset, config["clusters"][_dataset], _included_names, CONFIG_DIR,
    )
    for _cluster, _row in load_cluster_metadata(_metadata_path).items():
        CLUSTER_METADATA[(_dataset, _cluster)] = _row

    RESULTS_DIRS[_dataset] = os.path.join(RESULTS_DIR_BASE, _dataset)
    SCE2G_CONFIGS[_dataset] = build_scE2G_config(
        config, config["scE2G_dir"], CELL_CLUSTERS_TABLES[_dataset], RESULTS_DIRS[_dataset],
    )

# ---------------------------------------------------------------------------
# KNOWN LIMITATION -- read before touching: two of ABC's own rule files
# (ENCODE_rE2G/ABC/workflow/rules/macs2.smk, predictions.smk) render shell
# commands with a plain, un-namespaced `{RESULTS_DIR}`/`{SCRIPTS_DIR}`
# placeholder rather than a proper `{params...}`/`{output...}` reference.
# Snakemake's shell-string formatting for those resolves them from this
# pipeline's TOP-LEVEL namespace (empirically confirmed: removing our own
# same-named global here changed a working `--outdir` into a hard
# `NameError`, not a fix), not from a copy of ABC's own module-local scope --
# so nested `use rule ... as` re-exporting (scE2G's own import of
# ENCODE_rE2G, plus this pipeline's own per-dataset re-export on top of that)
# does not properly isolate this variable per dataset. Since every Slurm
# worker node re-parses this entire Snakefile (looping over every dataset in
# DATASETS), a run spanning 2+ datasets would have every macs2/predictions
# job pick up whichever dataset happened to be LAST in that loop -- silently
# writing peaks/predictions to another dataset's directory for every other
# dataset. This is upstream scE2G's bug, not fixable from here; tracked in
# the scE2G-cwd-footgun-task.md doc alongside the SCRIPTS_DIR/os.getcwd()
# issue (same root cause: plain globals referenced via shell placeholders
# don't survive being imported as a module).
#
# Workaround, safe ONLY for single-dataset REAL EXECUTION (dry-run DAG
# building across multiple datasets is unaffected -- shell commands are
# never rendered/executed during a dry run): expose RESULTS_DIR as the last
# dataset processed above. Do not rely on this for a real multi-dataset run
# until the upstream fix lands -- run those one dataset at a time for now.
# ---------------------------------------------------------------------------
if len(DATASETS) > 1:
    print(
        f"[QC-and-Predictions] WARNING: {len(DATASETS)} datasets in this run "
        f"({DATASETS}). Real execution of macs2/predictions rules is only "
        "verified correct for a single dataset per run -- see the "
        "KNOWN LIMITATION comment in common.smk before running this for real "
        "(not just -n) with multiple datasets."
    )
if DATASETS:
    RESULTS_DIR = RESULTS_DIRS[DATASETS[-1]]

# ---------------------------------------------------------------------------
# One scE2G `module` instance per dataset, generated as literal per-dataset
# .smk files rather than a Python for-loop directly containing `module`/`use
# rule` statements: those are Snakemake DSL keywords with their own grammar
# (not plain Python expressions), and neither the module name nor the `use
# rule ... as` prefix accepts an f-string/computed value -- confirmed by
# Snakemake 9.16.3 rejecting both `module f"scE2G_{dataset}":` ("Expected
# name or colon after module keyword") and `use rule * from X as
# f"sce2g_{dataset}_*"` ("Expecting rulename modifying pattern"). `include:`
# has no such restriction (it just textually splices in a file), so each
# dataset gets a small generated file with the dataset name baked in as a
# literal identifier, and the main Snakefile includes one per dataset.
# ---------------------------------------------------------------------------
# sce2g_modules: false skips importing scE2G entirely -- for a PORTAL-ONLY pass
# (`--config sce2g_modules=false` with the `portal_reformat` target), where the
# only rules that run are the 3 reformat rules, whose inputs are scE2G outputs
# that already exist on disk as files. Snakemake matches those by path, so it
# needs no rule capable of producing them.
#
# Worth being precise about what this does and doesn't buy, because it is easy to
# mistake for a fix it isn't. It is a PARSE-COST optimisation: importing scE2G
# once per dataset is slow, and it drags in the ABC_BIOSAMPLES mtime-preservation
# dance in the top-level Snakefile (scE2G rewrites tmp/config_abc_biosamples.tsv
# at parse time on every run, which fools an mtime DAG into rerunning ~140 jobs).
# Skipping the import removes that hazard by construction rather than working
# around it. It is NOT the fix for spurious rebuilds in a full run -- that is
# --rerun-triggers mtime, which run_pipeline.py always passes.
#
# Safe only AFTER Phase 2 is complete for these clusters: with no scE2G rules
# loaded, a genuinely missing prediction file has no producing rule and Snakemake
# fails with MissingInputException. That is the desired behaviour -- loud, not a
# silent skip -- but it does mean this is strictly a post-Phase-2 operation.
def _config_bool(key, default):
    """Boolean config value that survives `--config key=false`.

    Snakemake parses --config values with a literal_eval and falls back to the
    raw string, so `--config sce2g_modules=false` yields the STRING "false" --
    which is truthy. A bare `config.get("sce2g_modules", True)` therefore silently
    ignored the flag: scE2G was imported anyway, and a pass meant to touch only
    the reformat rules planned the whole generate_atac_matrix -> compute_kendall
    -> arc_e2g cascade instead (measured 2026-08-28: 127 jobs for igvf3, 21 of
    them compute_kendall, where 1 was expected). Only `False` or `0` happened to
    work, which is not a contract anyone should have to know.

    Raises on an unrecognised value rather than guessing: a typo silently
    re-enabling a multi-hour cascade is exactly the failure this exists to stop.
    """
    value = config.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    text = str(value).strip().lower()
    if text in ("false", "0", "no", "off"):
        return False
    if text in ("true", "1", "yes", "on"):
        return True
    raise ValueError(
        f"config key {key!r} must be a boolean-ish value "
        f"(true/false/1/0/yes/no/on/off), got {value!r}"
    )


SCE2G_MODULES_ENABLED = _config_bool("sce2g_modules", True)

GENERATED_RULES_DIR = os.path.join(workflow.basedir, "rules", "generated")
os.makedirs(GENERATED_RULES_DIR, exist_ok=True)
SCE2G_MODULE_FILES = []
# Explicit path->dataset lookup (not implicit zip(DATASETS, SCE2G_MODULE_FILES)
# ordering) so the Snakefile's own inclusion loop can compute each dataset's
# ABC_BIOSAMPLES path without assuming these two lists stay in lockstep --
# see that loop's own comment for why it needs this at all.
SCE2G_MODULE_FILE_DATASETS = {}
for _dataset in DATASETS if SCE2G_MODULES_ENABLED else []:
    _module_name = f"scE2G_{_dataset}"
    _rule_prefix = f"sce2g_{_dataset}"
    _generated_path = os.path.join(GENERATED_RULES_DIR, f"{_module_name}.smk")
    with open(_generated_path, "w") as f:
        f.write(
            f'module {_module_name}:\n'
            f'    snakefile:\n'
            f'        os.path.join(config["scE2G_dir"], "workflow", "Snakefile")\n'
            f'    config:\n'
            f'        SCE2G_CONFIGS["{_dataset}"]\n\n'
            # plot_stats is excluded: it rebuilds all_qc_stats.tsv from scratch every
            # run, scoped only to clusters in THIS dataset's cell_clusters.tsv (i.e.
            # only clusters ever run through this pipeline) -- it can't preserve rows
            # for clusters processed outside this pipeline, or clusters untouched in
            # this particular invocation. qc_stats.smk's aggregate_qc_stats rule
            # replaces it, writing to the exact same output path
            # (RESULTS_DIR/qc_plots/all_qc_stats.tsv), so the imported (unexcluded)
            # `hover_plots` rule -- whose only input is that same path -- transparently
            # depends on our replacement instead, with no changes needed on its end.
            #
            # `exclude` MUST precede `as` here: Snakemake 9.16.3's parser only
            # recognizes the `exclude` keyword directly after `from MODULE`;
            # once inside the `as`-clause state, "exclude" is swallowed as if it
            # were literal text in the rename pattern instead of a keyword,
            # silently producing a mangled name modifier like
            # `sce2g_igvf10_*exclude plot_stats` and NOT actually excluding
            # anything (confirmed empirically -- `as PREFIX_* exclude RULE`
            # renamed every single imported rule with an "excludeplot_stats"
            # suffix and left plot_stats' own output paths still claimed by
            # the (renamed) original rule, silently shadowing our replacement
            # rules without even raising Snakemake's usual AmbiguousRuleException).
            f'use rule * from {_module_name} exclude plot_stats as {_rule_prefix}_*\n'
        )
    SCE2G_MODULE_FILES.append(_generated_path)
    SCE2G_MODULE_FILE_DATASETS[_generated_path] = _dataset
