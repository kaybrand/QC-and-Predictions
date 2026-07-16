"""
Rebuilds all_qc_stats.tsv for one dataset by scanning that dataset's own
results directory fresh every run -- NOT by merging with the previous
all_qc_stats.tsv. Replaces scE2G's own `plot_stats` rule (excluded when
importing scE2G -- see common.smk), which instead rebuilds strictly from
whatever clusters are in THIS pipeline's own cell_clusters.tsv (i.e. only
clusters ever run through this pipeline), so it can't reflect clusters
processed some other way, or preserve clusters untouched in a given
invocation of this pipeline.

Row uniqueness is (cluster, model_name); scE2G's own
get_stats_per_model_per_cluster rule already tags every per-model stats.tsv
with those two columns (see
scE2G/workflow/scripts/prediction_qc/get_stats_per_cluster.R), so grouping
is done on file CONTENT, not on path parsing.

Every cluster directory on disk (any sibling of qc_plots/config/tmp under
the dataset's results directory) contributes its most-recently-generated
per-model stats.tsv, whether that file was just written by this run or is
untouched from a previous one:
  - clusters processed this run get fresh rows (their stats.tsv was just rewritten)
  - clusters NOT processed this run but still present on disk keep their existing
    rows automatically (their stats.tsv wasn't touched, so the same content
    is picked up by the scan)
  - clusters whose directory no longer exists disappear entirely (nothing to scan)

If a model's score_threshold ever changes, stale threshold-suffixed stats.tsv
files can accumulate on disk alongside the current one for the same
(cluster, model) -- only the most recently modified file per (cluster,
model) is kept, so a superseded threshold's stats never leak back in.

scE2G's own get_stats_per_cluster.R writes cell_count=0 whenever a cluster's
feature-generation method isn't Kendall/ARC (see that script's `cell_path`
directory-check) -- which is every ATAC-only cluster, since cell_count.txt is
only ever produced alongside Kendall/ARC's RNA-based features. That 0 is not
a missing measurement -- the actual cell count is already known, from the
same QC guide / filtered_cell_subsample_metrics.tsv resolve_exclusions.py
uses for exclusion thresholds -- so it's patched in here from disk, purely by
re-deriving it from that same metrics file, with no dependency on this run's
own config (keeps this script's existing "works for any cluster dir on disk,
not just this run's configured ones" property).
"""

import csv
import glob
import os
import sys

EXCLUDED_DIRS = {"qc_plots", "config", "tmp"}


def _cell_count_from_metrics(plots_dir, cluster):
    """Sum n_cells across subsamples in filtered_cell_subsample_metrics.tsv --
    same file + same column resolve_exclusions.py sums for exclusion thresholds
    (validated there to exactly match scE2G's own cell_count for clusters where
    scE2G does compute it). Returns None if the file doesn't exist (e.g. a
    non-default QC guide was used, or the cluster predates this pipeline)."""
    metrics_path = os.path.join(plots_dir, cluster, "filtered_cell_subsample_metrics.tsv")
    if not os.path.exists(metrics_path):
        return None
    with open(metrics_path) as f:
        return sum(int(row["n_cells"]) for row in csv.DictReader(f, delimiter="\t"))


def _cluster_dirs(dataset_dir):
    """Every sibling of qc_plots/config/tmp under a dataset's results directory,
    with NO quality/threshold filtering of any kind. This is a sanity-check view:
    a cluster excluded elsewhere in this pipeline for low cell_count/fragments/umi
    (see resolve_exclusions.py) must still show up here for comparison as long as
    it has ever actually been processed (has a directory + stats.tsv on disk) --
    the whole point is being able to see, e.g., a 15-cell cluster sitting visibly
    below every threshold line in these plots, not have it silently vanish. Do not
    add cell_count/fragments_total/etc. filtering here."""
    return [
        d for d in os.listdir(dataset_dir)
        if d not in EXCLUDED_DIRS and os.path.isdir(os.path.join(dataset_dir, d))
    ]


def _latest_stats_row_per_cluster_model(dataset_dir, plots_dir):
    latest = {}  # (cluster, model_name) -> (mtime, row)
    for cluster_dir in _cluster_dirs(dataset_dir):
        pattern = os.path.join(dataset_dir, cluster_dir, "*", "scE2G_predictions_threshold*_stats.tsv")
        for path in glob.glob(pattern):
            mtime = os.path.getmtime(path)
            with open(path) as f:
                row = next(csv.DictReader(f, delimiter="\t"), None)
            if row is None:
                continue
            if int(row.get("cell_count", 0)) == 0:
                patched = _cell_count_from_metrics(plots_dir, row["cluster"])
                if patched:
                    row["cell_count"] = str(patched)
            key = (row["cluster"], row["model_name"])
            if key not in latest or mtime > latest[key][0]:
                latest[key] = (mtime, row)
    return [row for _, row in latest.values()]


def build_all_qc_stats(dataset_dir, plots_dir, out_path):
    rows = _latest_stats_row_per_cluster_model(dataset_dir, plots_dir)
    header = list(rows[0].keys()) if rows else []

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp_path = f"{out_path}.tmp.{os.getpid()}"
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["cluster"], r["model_name"])):
            writer.writerow(row)
    os.replace(tmp_path, out_path)  # atomic: concurrent readers never see a partial file
    return out_path


if __name__ == "__main__":
    dataset_dir, plots_dir, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
    build_all_qc_stats(dataset_dir, plots_dir, out_path)
