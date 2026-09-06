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


def scE2G_chr_sizes(config):
    """Absolute path to the chromosome sizes file scE2G ITSELF uses.

    DERIVED from scE2G's own configs, never taken from our pipeline config.
    Sorting fragments into a different chromosome order than scE2G's own
    frag_to_tagAlign uses is precisely the bug class this whole sorting detour
    exists to avoid, so a second independently-settable copy of this path would
    be a way to reintroduce it silently.

    The lookup mirrors scE2G's: `chr_sizes` is declared in
    ENCODE_rE2G/config/config.yaml as a path relative to `encode_re2g_dir`
    ("reference/GRCh38_EBV.no_alt.chrom.sizes.tsv"), and `encode_re2g_dir` is
    declared in scE2G's config/config.yaml relative to scE2G_dir.
    """
    import yaml

    scE2G_dir = config["scE2G_dir"]
    with open(os.path.join(scE2G_dir, "config", "config.yaml")) as f:
        sce2g_cfg = yaml.safe_load(f) or {}
    encode_dir = sce2g_cfg.get("encode_re2g_dir", "ENCODE_rE2G")
    if not os.path.isabs(encode_dir):
        encode_dir = os.path.join(scE2G_dir, encode_dir)

    with open(os.path.join(encode_dir, "config", "config.yaml")) as f:
        encode_cfg = yaml.safe_load(f) or {}
    rel = encode_cfg.get("chr_sizes")
    if not rel:
        raise ValueError(
            f"no `chr_sizes` in {encode_dir}/config/config.yaml -- cannot determine the "
            "chromosome order scE2G expects, and guessing it would risk sorting fragments "
            "into an order scE2G's own rules disagree with."
        )
    path = rel if os.path.isabs(rel) else os.path.join(encode_dir, rel)
    if not os.path.exists(path):
        raise ValueError(f"scE2G's chr_sizes points at {path}, which does not exist")
    return path


def bed_split_sort_script(config):
    """scE2G's own bedSplitSort.sh -- split by chromosome, sort each by position,
    concatenate in the sizes file's order. Reused rather than reimplemented so
    fragments get the SAME ordering guarantee scE2G's frag_to_tagAlign already
    relies on (added there in scE2G 299ca9e / PR #94).

    Read-only use of a file in the scE2G checkout. That checkout is a separate
    shared worktree and must never be modified from here.
    """
    path = os.path.join(config["scE2G_dir"], "workflow", "scripts", "bedSplitSort.sh")
    if not os.path.exists(path):
        raise ValueError(
            f"scE2G's bedSplitSort.sh not found at {path}. It was added in scE2G PR #94; "
            "an older checkout will not have it."
        )
    return path


def sorted_fragments_dir(config):
    """Where coordinate-sorted copies of `atac_frag_file` fragments live.

    Defaults to "{download_root}_sorted" so it sits beside the download tree on
    the same (scratch) filesystem, and is overridable via
    config["sorted_fragments_dir"]. Returns None when neither is available --
    only datasets with atac_frag_file clusters need this at all.
    """
    explicit = config.get("sorted_fragments_dir")
    if explicit:
        return explicit
    root = config.get("download_root")
    return f"{root.rstrip('/')}_sorted" if root else None


def sorted_frag_file(config, dataset, cluster):
    """The coordinate-sorted fragments file the pipeline hands to scE2G for an
    `atac_frag_file` cluster, in place of the raw download.

    WHY THIS EXISTS. scE2G's frag_to_norm_bigWig pipes fragments straight into
    `bedtools genomecov -bg`, which requires chromosome-GROUPED input, and only
    sorts the resulting bedGraph afterwards -- which cannot repair overlapping
    intervals genomecov has already emitted from interleaved input. The CATlas
    fragment files are inconsistently sorted upstream: positions are sorted
    within each run of a chromosome, but ~18% of files revisit chromosomes
    later, and those fail with "Error - overlapping regions in bedGraph".
    scE2G's `fragments_preprocessed: False` does NOT help -- that branch
    (process_fragment_file) filters to chromosomes in chr_sizes and never
    sorts. Verified 2026-09-03; scE2G's own frag_to_tagAlign DOES sort the same
    input via bedSplitSort.sh, so this is an inconsistency in scE2G, tracked for
    an upstream PR.

    Defined here rather than in either caller because two places must agree
    exactly: the rule that WRITES this file (prefiltered_fragments.smk's
    sort_atac_fragments) and write_scE2G_config.py, which puts the path into
    scE2G's cell_clusters TSV as its atac_frag_file. A mismatch would leave
    Snakemake unable to find a producing rule for scE2G's declared input, i.e. a
    MissingInputException at DAG build.
    """
    base = sorted_fragments_dir(config)
    if not base:
        raise ValueError(
            f"{dataset}/{cluster} uses atac_frag_file, so it needs sorted fragments, but neither "
            "config['sorted_fragments_dir'] nor config['download_root'] is set to derive a "
            "location from."
        )
    return os.path.join(base, dataset, f"{cluster}_fragments.sorted.tsv.gz")


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
