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
gtf: {gtf}

# hg38 ONLY. This is why mouse is out of scope: 369 of the 627 CATlas
# pseudobulks are Mus musculus and cannot be run against these annotations.
# A catlas-mouse dataset needs its own build before it can be generated.

# data_dir is READ-ONLY on this fork (QC-guide plots/datatables inputs). CATlas
# has no QC-guide step, so in practice only
# {{data_dir}}/plots/{dataset}/{{cluster}}/filtered_cell_subsample_metrics.tsv
# is read from here; multiome_data/datatables are never populated.
data_dir: {data_dir}

# Never invoked: every cluster sets atac_frag_file, so pseudobulks_dir()/
# atac_fragment_file are unreachable for this dataset. Required top-level key.
pseudobulks_root: {pseudobulks_root}

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
  # from filtered_cell_subsample_metrics.tsv, which this pipeline's own
  # prefiltered_cell_subsample_metrics rule produces -- so on a first parse no
  # cluster has stats yet.
  predictions_on_everything_but_do_not_upload: true
  auto_thresholds:
    # Same bar as every other dataset in this project. min_umi_count is
    # irrelevant for ATAC-only clusters (umi_count is None and the comparison
    # is skipped), kept only so the block matches the shared convention.
    min_cell_count: {min_cell_count}
    min_fragments_total: {min_fragments_total}
    min_umi_count: 1000000

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
    p.add_argument("--repo", default=None, help="repo root; defaults to two levels above this script")
    p.add_argument("--pipeline-mode", default="local_only", choices=("local_only", "default"),
                   help="local_only (default) makes no Portal contact and requests no reformat "
                        "targets -- the right mode for a first predictions-only run")
    p.add_argument("--min-cell-count", type=int, default=100)
    p.add_argument("--min-fragments-total", type=int, default=2_000_000)
    p.add_argument("--max-cell-count", type=int, default=20000)
    p.add_argument("--scE2G-version", default="1.2")
    p.add_argument("--scE2G-dir",
                   default="/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G")
    p.add_argument("--gtf", default="/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/"
                                   "uniformly_processed/GENCODEv43/IGVFFI9573KOZR.gtf.gz")
    return p.parse_args(argv)


def load_clusters(conn, dataset):
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
        if "upload_status:" in reasons:
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
    clusters, skipped = load_clusters(conn, args.dataset)
    conn.close()

    if not clusters:
        raise SystemExit(
            f"no usable clusters for dataset {args.dataset!r} in {args.state_db}. "
            "Has the download run? (download_portal_pseudobulks.py --family catlas)"
        )

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

    mode_note = (
        "local_only: no IGVF Portal contact at any stage, and no reformat targets "
        "(REFORMAT_ELIGIBLE_CLUSTERS is empty by construction). Predictions + QC only."
        if args.pipeline_mode == "local_only" else
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
        gtf=args.gtf,
        data_dir=os.path.join(repo, "data"),
        pseudobulks_root=f"/scratch/users/{os.environ.get('USER','USER')}/UNUSED_no_pseudobulks_root_for_catlas",
        output_dir=os.path.join(repo, "results"),
        sce2g_version=args.scE2G_version,
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
