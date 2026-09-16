#!/usr/bin/env python3
"""Generate a CATlas *_pipeline_config.yaml from the download ledger.

    python workflow/scripts/generate_catlas_config.py \
        --state-db resources/igvf_metadata_state.db \
        --dataset catlas-human \
        --out workflow/config/catlas-human_pipeline_config.yaml

Hand-authoring 257 cluster entries is not viable, and the ledger already holds
exactly what a cluster entry needs: the cluster name (parsed from the Portal
alias) and the verified local path of its fragments file. This reads that rather
than globbing the download tree, so a file that is present but NOT md5-verified
can never silently become a pipeline input.

WHAT IT SKIPS, and says so on stderr:
  - download_state != 'done'    -- not verified, so not an input.
  - review_reasons containing "upload_status:" -- the Portal REJECTED the file.
    For CATlas that means content_md5sum duplication: the bytes are fine and the
    file downloads cleanly, but another cluster already has identical fragments.
    Including it would produce a second prediction set over the same data.
    `human-supercluster-Bergmann_glia` is exactly this case (duplicate of
    `Bgl_51`), which is why 258 downloaded files yield 257 clusters.

Deliberately NOT written: `models` is fixed to scATAC_powerlaw_v3 for every
cluster, because CATlas is scATAC-only by construction (no RNA anywhere in the
Yang Li pseudobulks). If that ever changes it becomes a per-cluster field.
"""

import argparse
import os
import sqlite3
import sys

TEMPLATE_HEADER = '''\
# QC-and-Predictions pipeline config for "{dataset}" -- GENERATED FILE.
#
# Regenerate with:
#     python workflow/scripts/generate_catlas_config.py \\
#         --state-db {state_db_rel} \\
#         --dataset {dataset} \\
#         --out {out_rel}
#
# Do not hand-edit the clusters block: it is derived from the download ledger's
# portal_files table and will be overwritten. Edit the generator instead.
#
# {n_clusters} scATAC-only CATlas (Yang Li lab) pseudobulks, fragments pulled
# pre-filtered/QC'd from the IGVF Data Portal as curated-set PseudobulkSets --
# NOT this pipeline's usual raw-archive + QC-guide-barcode-list flow. Every
# cluster therefore sets atac_frag_file directly instead of
# pseudobulk_annotation/qc_guide. See
# workflow/config/example_pipeline_config.yaml for the general schema.
{skipped_note}
scE2G_dir: {sce2g_dir}

# hg38 ONLY (via scE2G_options.gene_TSS500/genes below). This is why mouse is
# out of scope: 369 of the 627 CATlas pseudobulks are Mus musculus and cannot be
# run against these annotations. A catlas-mouse dataset needs its own build.
#
# NO `gtf` KEY, deliberately. It is read only by filter_pseudobulks.smk's
# rna_count_matrix rule, for Ensembl-ID -> gene-symbol conversion on an RNA count
# matrix. That rule is gated on MULTIOME_CLUSTERS, which is empty for any dataset
# whose every cluster is models: [scATAC_powerlaw_v3] -- so it can never run here
# and the GTF would never be opened. (Note the ACCESSION IGVFFI9573KOZR is still
# cited in this dataset's manifests as a reference_files value, but that comes
# from constants in igvf_metadata/tables/, not from this config.)

# CATlas has no QC-guide step, so none of the QC-guide inputs data_dir normally
# supplies (plots/, datatables/, multiome_data/) exist or are read for this
# dataset. The quality gate reads scE2G's own qc_plots/all_qc_stats.tsv under
# output_dir instead. Still required as a top-level key, and still WRITTEN to:
# rule prefiltered_cell_subsample_metrics puts its per-cluster intermediate at
# {{data_dir}}/plots/{dataset}/{{cluster}}/filtered_cell_subsample_metrics.tsv,
# which aggregate_qc_stats consumes to fill cell_count/fragments_total into
# all_qc_stats.tsv. That write is why data_dir points inside THIS worktree --
# it must never be the QC_pseudobulks tree.
data_dir: {data_dir}

# WHERE THE CATLAS FRAGMENTS LIVE, and the single source of truth for it.
# download_portal_pseudobulks.py --pipeline-config reads this to decide where to
# PUT downloads, and every atac_frag_file below resolves under it -- so the
# fetch destination and the pipeline's inputs cannot drift apart.
#
# NOT `pseudobulks_root`: that key feeds pseudobulks_dir() -> the QC-guide
# atac_fragment_file/rna_count_matrix rules, which are never instantiated for
# atac_frag_file clusters and never will be for CATlas. It is deliberately
# absent from this config rather than set to a plausible-looking value, because
# pseudobulks_dir() appends "/{{dataset}}/pseudobulks" and CATlas's layout is
# "/{{dataset}}/fragments" -- a value here would resolve to a directory that
# does not exist while looking live.
download_root: {download_root}

# Where coordinate-SORTED copies of the fragments go. On $SCRATCH beside the
# download tree -- job data I/O, and these are full-size copies (~776 GB for all
# 257 human clusters), marked temp() so Snakemake reclaims them.
#
# Needed because scE2G does NOT sort fragments and assumes it has been done:
# frag_to_norm_bigWig pipes them straight into `bedtools genomecov -bg`, which
# requires chromosome-GROUPED input, and only sorts the resulting bedGraph --
# too late. ~18% of the CATlas files revisit chromosomes later in the file and
# fail with "Error - overlapping regions in bedGraph".
# `fragments_preprocessed: False` does not help; that branch only filters to
# chromosomes in chr_sizes. scE2G's own frag_to_tagAlign DOES sort the same
# input (via bedSplitSort.sh), so this is an inconsistency upstream -- tracked
# for a PR, after which prefiltered_fragments.smk's sort_atac_fragments and this
# key can both go away.
sorted_fragments_dir: {sorted_fragments_dir}

# Replaces the pre-rebase `results_dir_base`, which nothing reads any more --
# common.smk now derives RESULTS_DIR_BASE = {{output_dir}}/uniformly_processed.
# Keeps results under THIS worktree rather than the shared scE2G_dir tree; the
# three worktrees are separate forks with separate submission targets.
output_dir: {output_dir}

# preferred_assay_titles on all 627 yang-li CATlas PseudobulkSets, unanimous
# across both species. Feeds the "# Assays:" header of every reformatted file;
# without it they would claim "10x Multiome".
assay_title: snATAC-seq

scE2G_version: "{sce2g_version}"

# The conda env cache Snakemake must reuse. run_pipeline.py reads this key
# (falling back to its own DEFAULT_CONDA_PREFIX) and hard-refuses to run if the
# directory is absent -- Snakemake must NEVER be left to build envs itself: a
# fresh prefix reads as "every env needs building from scratch", and concurrent
# runs then race and corrupt each other's downloads.
#
# Stated here rather than left to that default so this dataset names its own
# prefix explicitly, instead of depending on a hardcoded path in a script shared
# with the other worktrees.
snakemake_conda_prefix: {conda_prefix}

# {pipeline_mode_note}
pipeline_mode: {pipeline_mode}

# CATlas-local state.db -- NEVER the shared one under
# QC_pseudobulks-igvf-portal-submission.
igvf:
  state_db_path: {state_db}

scE2G_options:
  max_cell_count: {max_cell_count}
  make_IGV_tracks: true
  gene_TSS500: resources/genome_annotations/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.TSS500bp.bed
  genes: resources/genome_annotations/CollapsedGeneBound.hg38.GENCODEv43GeneSymbol.bed

exclusion:
  user_specified: []
  # Renamed from the deprecated `process_excluded_no_upload` (still accepted
  # with a notice, but the old name now prints one on every parse).
  #
  # true: PROCESS every cluster -- predictions, ATAC_norm.bw IGV tracks -- and
  # let the thresholds below decide only what is UPLOAD-eligible. Wanted here
  # because the Prediction QC Metrics are the basis for deciding which CATlas
  # clusters to keep, and a cluster that was never processed produces no
  # metrics to judge.
  #
  # It is also what makes the first run work at all: these clusters' stats come
  # from scE2G's own {{output_dir}}/uniformly_processed/{dataset}/qc_plots/
  # all_qc_stats.tsv, which does not exist until scE2G has run -- so on a first
  # parse no cluster has stats and every one reports pending_qc_stats. With this
  # flag false, nothing would ever be processed and nothing would ever be
  # measured.
  predictions_on_everything_but_do_not_upload: true
  auto_thresholds:
    # Same bar as every other dataset in this project. No min_umi_count: these
    # clusters are ATAC-only, there is no RNA anywhere in the CATlas pseudobulks,
    # and resolve_exclusions reports umi_count as None for them so the threshold
    # would never be compared. Naming it would imply a measurement that does not
    # exist here.
    min_cell_count: {min_cell_count}
    min_fragments_total: {min_fragments_total}

clusters:
  {dataset}:
'''

CLUSTER_ENTRY = """\
    {cluster}:
      models: [scATAC_powerlaw_v3]
      atac_frag_file: {path}
"""


def log(msg):
    print(f"[generate_catlas_config] {msg}", file=sys.stderr)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--state-db", required=True)
    p.add_argument("--dataset", default="catlas-human",
                   help="ledger dataset label AND the config's dataset key (default catlas-human)")
    p.add_argument("--out", required=True)
    p.add_argument(
        "--include-portal-flagged", action="store_true",
        help="also include clusters whose fragments file the Portal REJECTED "
             "(review_reasons contains 'upload_status:'). For CATlas that rejection is "
             "content_md5sum duplication -- the bytes downloaded and md5-verified fine, so the "
             "file is usable input. Off by default because a future flag might mean something "
             "genuinely broken rather than merely duplicated. Turned ON for catlas-human as of "
             "2026-09-04: duplicates are accepted (the submitting lab says either version may be "
             "used), so excluding a usable supercluster while keeping a comparable one was an "
             "inconsistency in the delivered set. See "
             ".baton/2026-09-04-duplicate-fragment-files.md.",
    )
    p.add_argument(
        "--clusters", default=None,
        help="comma-separated cluster names to restrict to, for a small rehearsal config. The "
             "dataset KEY is unchanged, so a test config writes into the same results tree as the "
             "full one and its outputs carry forward instead of being thrown away. Every name must "
             "exist in the ledger; a typo is an error, not a silently smaller batch.",
    )
    p.add_argument(
        "--download-root", default=None,
        help="where the CATlas fragments live; written into the config as `download_root`, which "
             "is what download_portal_pseudobulks.py --pipeline-config then reads to decide where "
             "to fetch INTO. Defaults to the root already recorded in --out (regeneration) or, "
             "failing that, inferred from the ledger's own local_path values.",
    )
    p.add_argument("--repo", default=None, help="repo root; defaults to two levels above this script")
    p.add_argument(
        "--pipeline-mode", default="local_only", choices=("local_only", "default", "reformat_only"),
        help="local_only (default) makes no Portal contact and requests no reformat targets -- the "
             "right mode for a first predictions-only run. reformat_only warms the CellAnnotation "
             "cache and reformats like default, but skips IGVF manifest generation.",
    )
    p.add_argument("--min-cell-count", type=int, default=100)
    p.add_argument("--min-fragments-total", type=int, default=2_000_000)
    p.add_argument("--max-cell-count", type=int, default=20000)
    p.add_argument(
        "--conda-prefix",
        default="/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/.snakemake/conda",
        help="written into the config as snakemake_conda_prefix. Must be an ALREADY-POPULATED env "
             "cache; Snakemake is never allowed to build envs itself.",
    )
    p.add_argument("--scE2G-version", default="1.3.1")
    p.add_argument("--scE2G-dir",
                   default="/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G")
    return p.parse_args(argv)


def _existing_download_root(path):
    """`download_root` from an already-generated config, or None.

    Read with a line scan rather than a YAML parse so this stays usable when the
    file is absent, half-written, or has a clusters block too large to be worth
    loading just to recover one string.
    """
    if not os.path.exists(path):
        return None
    with open(path) as f:
        for line in f:
            if line.startswith("download_root:"):
                return line.split(":", 1)[1].strip() or None
            if line.startswith("clusters:"):
                break
    return None


def load_clusters(conn, dataset, include_portal_flagged=False):
    """(clusters, skipped) from the ledger. clusters is [(name, local_path)]."""
    rows = conn.execute(
        "SELECT annotation, local_path, download_state, review_reasons, accession "
        "FROM portal_files WHERE dataset = ? AND content_type = 'fragments' "
        "ORDER BY annotation",
        (dataset,),
    ).fetchall()
    clusters, skipped = [], []
    for annotation, local_path, state, reasons, accession in rows:
        reasons = reasons or ""
        if state != "done":
            skipped.append((annotation, accession, f"download_state={state}"))
            continue
        if "upload_status:" in reasons and not include_portal_flagged:
            # Portal-rejected. For CATlas this is content_md5sum duplication --
            # the file is byte-valid but duplicates another cluster's fragments.
            skipped.append((annotation, accession, f"portal-rejected ({reasons})"))
            continue
        if not local_path:
            skipped.append((annotation, accession, "no local_path"))
            continue
        clusters.append((annotation, local_path))
    return clusters, skipped


def main(argv=None):
    args = parse_args(argv)
    repo = args.repo or os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
    conn = sqlite3.connect(args.state_db)
    clusters, skipped = load_clusters(conn, args.dataset, args.include_portal_flagged)
    conn.close()

    if not clusters:
        raise SystemExit(
            f"no usable clusters for dataset {args.dataset!r} in {args.state_db}. "
            "Has the download run? (download_portal_pseudobulks.py --family catlas)"
        )

    if args.clusters:
        wanted = [c.strip() for c in args.clusters.split(",") if c.strip()]
        available = {c for c, _p in clusters}
        missing = [c for c in wanted if c not in available]
        if missing:
            # Loud, because a silently-smaller rehearsal batch would still LOOK
            # like it passed while never touching the cluster you meant to test.
            raise SystemExit(
                f"--clusters names {len(missing)} cluster(s) not usable in the ledger for "
                f"{args.dataset!r}: {missing}. They may be misspelled, not downloaded, or "
                "portal-rejected (see the skip list above)."
            )
        order = {c: i for i, c in enumerate(wanted)}
        clusters = sorted((c for c in clusters if c[0] in order), key=lambda cp: order[cp[0]])
        log(f"--clusters restricted the batch to {len(clusters)}: {[c for c, _ in clusters]}")

    log(f"{len(clusters)} cluster(s) from the ledger; {len(skipped)} skipped")
    for name, accession, why in skipped:
        log(f"  SKIPPED {name} ({accession}): {why}")

    missing = [(c, p) for c, p in clusters if not os.path.exists(p)]
    if missing:
        # The ledger says done, so a vanished file means $SCRATCH was purged or
        # the tree moved. Loud, because the alternative is a MissingInputException
        # hundreds of jobs into a run.
        raise SystemExit(
            f"{len(missing)} ledger-verified fragments file(s) are not on disk, e.g. "
            f"{missing[0][1]} -- was $SCRATCH purged? Re-run the downloader (it resumes)."
        )

    # download_root precedence: explicit flag > the root already in --out (so
    # regenerating never silently relocates a dataset) > inferred from the
    # ledger. Whatever wins is then CHECKED against the ledger, because the
    # config is what the downloader will fetch into and the atac_frag_file paths
    # below have to resolve under it.
    #
    # Ledger layout is {root}/{dataset}/fragments/{file} (see
    # portal_files.catlas_resolve_scope), so the root is three levels up.
    ledger_roots = {os.path.dirname(os.path.dirname(os.path.dirname(p))) for _c, p in clusters}
    if len(ledger_roots) != 1:
        raise SystemExit(
            f"ledger local_paths span {len(ledger_roots)} download roots "
            f"({sorted(ledger_roots)}); download_root can only name one. "
            "Re-download into a single tree."
        )
    ledger_root = ledger_roots.pop()
    download_root = args.download_root or _existing_download_root(args.out) or ledger_root
    if download_root != ledger_root:
        raise SystemExit(
            f"download_root {download_root!r} disagrees with where the ledger says the fragments "
            f"actually are ({ledger_root!r}). The config's download_root is BOTH the fetch "
            "destination and the root of every atac_frag_file path, so a mismatch would generate "
            "a config pointing at files that are not there. Re-download into the new root, or "
            "pass --download-root to match the ledger."
        )

    if args.pipeline_mode == "local_only":
        mode_note = (
            "local_only: no IGVF Portal contact at any stage, and no reformat targets "
            "(REFORMAT_ELIGIBLE_CLUSTERS is empty by construction). Predictions + QC only."
        )
    elif args.pipeline_mode == "reformat_only":
        mode_note = (
            "reformat_only: warms the CellAnnotation cache and reformats predictions into "
            "portal format, same as default -- but the driver (run_pipeline.py) skips IGVF "
            "manifest generation and audit (stages 3-4)."
        )
    else:
        mode_note = (
            "default: warms the CellAnnotation cache and builds manifests. Requires "
            "every cluster to be Portal-annotated."
        )
    skipped_note = ""
    if skipped:
        skipped_note = "#\n# SKIPPED, and why (see the generator's stderr for the full list):\n" + "".join(
            f"#   {name}: {why}\n" for name, _acc, why in skipped
        )

    body = TEMPLATE_HEADER.format(
        dataset=args.dataset,
        n_clusters=len(clusters),
        skipped_note=skipped_note,
        state_db=args.state_db if os.path.isabs(args.state_db) else os.path.join(repo, args.state_db),
        state_db_rel=args.state_db,
        out_rel=args.out,
        sce2g_dir=args.scE2G_dir,
        data_dir=os.path.join(repo, "data"),
        download_root=download_root,
        sorted_fragments_dir=f"{download_root.rstrip('/')}_sorted",
        output_dir=os.path.join(repo, "results"),
        sce2g_version=args.scE2G_version,
        conda_prefix=args.conda_prefix,
        pipeline_mode=args.pipeline_mode,
        pipeline_mode_note=mode_note,
        max_cell_count=args.max_cell_count,
        min_cell_count=args.min_cell_count,
        min_fragments_total=args.min_fragments_total,
    )
    body += "".join(CLUSTER_ENTRY.format(cluster=c, path=p) for c, p in clusters)

    with open(args.out, "w") as f:
        f.write(body)
    log(f"wrote {args.out}: dataset={args.dataset}, {len(clusters)} clusters, "
        f"pipeline_mode={args.pipeline_mode}, "
        f"gate >={args.min_cell_count} cells / >={args.min_fragments_total} fragments")
    return 0


if __name__ == "__main__":
    sys.exit(main())
