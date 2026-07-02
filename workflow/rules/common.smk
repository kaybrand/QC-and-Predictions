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

# Repo root (the directory containing this repo's Snakefile/plots/datatables/...),
# derived from the location of the top-level workflow/Snakefile rather than the
# current working directory, so `plots/` and `datatables/` resolve correctly
# regardless of where `snakemake` is invoked from.
WDIR = os.path.dirname(workflow.basedir)

sys.path.insert(0, os.path.join(workflow.basedir, "scripts"))
from resolve_exclusions import resolve_exclusions  # noqa: E402
from write_scE2G_config import (  # noqa: E402
    write_cell_clusters_table,
    write_cluster_metadata_table,
    load_cluster_metadata,
)

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
    return os.path.join(WDIR, "multiome_data", dataset)


def pseudobulks_dir(dataset):
    return os.path.join(config["pseudobulks_root"], dataset, "pseudobulks")


# ---------------------------------------------------------------------------
# Exclusion resolution -- runs once at parse time, before any rule is defined.
# Every set below is keyed by (dataset, cluster) tuples.
# ---------------------------------------------------------------------------
INCLUDED_CLUSTERS, UPLOAD_ELIGIBLE_CLUSTERS, EXCLUDED_CLUSTERS, CLUSTER_STATS = resolve_exclusions(config, WDIR)

if EXCLUDED_CLUSTERS:
    print(f"[QC-and-Predictions] Excluding clusters this run: {sorted(EXCLUDED_CLUSTERS)}")
if not config.get("exclusion", {}).get("process_excluded_no_upload", False):
    print(f"[QC-and-Predictions] Processing clusters: {sorted(INCLUDED_CLUSTERS)}")
else:
    print(f"[QC-and-Predictions] Processing (incl. excluded, no-upload) clusters: {sorted(INCLUDED_CLUSTERS)}")
    print(f"[QC-and-Predictions] Upload-eligible clusters: {sorted(UPLOAD_ELIGIBLE_CLUSTERS)}")

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
RESULTS_DIR_BASE = os.path.join(config["scE2G_dir"], "results", "uniformly_processed")

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
GENERATED_RULES_DIR = os.path.join(workflow.basedir, "rules", "generated")
os.makedirs(GENERATED_RULES_DIR, exist_ok=True)
SCE2G_MODULE_FILES = []
for _dataset in DATASETS:
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
            f'use rule * from {_module_name} as {_rule_prefix}_*\n'
        )
    SCE2G_MODULE_FILES.append(_generated_path)
