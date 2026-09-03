"""Shared repo-relative path resolution for output_dir.

Both workflow/rules/common.smk (Snakemake context, has workflow.basedir) and
workflow/scripts/manage_igvf_metadata.py (standalone CLI, no Snakemake
context) need to resolve a possibly-relative output_dir against the repo
root identically -- this is the one place that resolution rule is defined,
so the two call sites can't drift apart.
"""

import os


def resolve_repo_relative(path, repo_root):
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(repo_root, path))


def repo_root_from_script(script_path):
    # workflow/scripts/<this file> -> repo root is two directories up.
    return os.path.abspath(os.path.join(os.path.dirname(script_path), "..", ".."))


def results_dir_base(config, repo_root):
    """The dataset-AGNOSTIC results parent, "{output_dir}/uniformly_processed".

    Added here rather than recomputed per caller because there are now two:
    common.smk's RESULTS_DIR_BASE, and resolve_exclusions.py, which needs it to
    find a dataset's qc_plots/all_qc_stats.tsv. If those two disagreed, the
    quality gate would read a different tree than the pipeline writes -- and it
    would fail open (no stats file found -> nothing gated), which is the silent
    direction.
    """
    return os.path.join(
        resolve_repo_relative(config.get("output_dir", "./results"), repo_root),
        "uniformly_processed",
    )
