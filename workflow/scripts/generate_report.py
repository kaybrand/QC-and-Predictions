#!/usr/bin/env python
"""Coverage report: one row per cluster covering every dataset that has a
cluster_stats table (written at Snakemake parse time by
resolve_exclusions.write_cluster_stats_table, common.smk) -- this is the
authoritative in-scope universe, since igvf7/igvf9/igvf13/the future-format
row never get a config and therefore never get a stats file.

has_cell_annotation is read from the SAME live igvf_metadata.state cache
common.smk's REFORMAT_ELIGIBLE_CLUSTERS gate and reformat.smk's
portal_cell_metadata() use -- not cell_annotations_by_dataset_cluster.tsv,
which is a manual preview snapshot, not an authoritative source (see
cell_annotations.py's docstring). This keeps the report's has_cell_annotation
column always consistent with what will actually get reformatted.

Safe to re-run at any point (before any real execution -- cluster_stats is
written even on `snakemake -n` -- or partway through the 13 per-dataset real
runs) since it only reads existing files, never mutates anything.
"""

import argparse
import csv
import glob
import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cell_annotations import annotation_lookup_key  # noqa: E402
from igvf_metadata import state as igvf_state  # noqa: E402

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DEFAULT_CONFIG_DIR = os.path.join(SCRIPT_DIR, "..", "config")
DEFAULT_STATE_DB = os.path.join(REPO_ROOT, "resources", "igvf_metadata_state.db")

REPORT_HEADER = [
    "dataset", "cluster", "cell_count", "fragments_total", "umi_count",
    "quality_pass", "exclusion_reason", "has_cell_annotation",
    "igvf_manifest_excluded", "predictions_generated", "reformatted", "status",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", default=os.path.join(REPO_ROOT, "results"))
    p.add_argument(
        "--state-db", default=DEFAULT_STATE_DB,
        help="igvf.state_db_path -- the live IGVF-Portal-fetch cache (cell_metadata.refresh_if_stale), "
        "NOT cell_annotations_by_dataset_cluster.tsv, which is a manual preview snapshot only",
    )
    p.add_argument("--configs-dir", default=DEFAULT_CONFIG_DIR)
    p.add_argument("--out", default=None, help="default: {output_dir}/report.tsv")
    return p.parse_args()


def load_cluster_stats(output_dir):
    """{(dataset, cluster): {cell_count, fragments_total, umi_count, reason}}
    across every {dataset}_cluster_stats.tsv found -- the authoritative
    162(-ish)-cluster in-scope universe for this round."""
    rows = {}
    stats_dir = os.path.join(output_dir, "cluster_stats")
    for path in sorted(glob.glob(os.path.join(stats_dir, "*_cluster_stats.tsv"))):
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                rows[(row["dataset"], row["cluster"])] = row
    return rows


def load_cluster_configs(configs_dir, datasets):
    """{(dataset, cluster): cluster_cfg} for every dataset that has a
    cluster_stats table -- needed for cell_annotation_key/igvf_manifest_excluded
    overrides, which live in the pipeline config, not the stats table."""
    cfgs = {}
    for dataset in sorted(datasets):
        path = os.path.join(configs_dir, f"{dataset}_pipeline_config.yaml")
        if not os.path.exists(path):
            print(f"[generate_report] WARNING: no pipeline config found for {dataset} at {path} -- skipping its rows", file=sys.stderr)
            continue
        with open(path) as f:
            config = yaml.safe_load(f)
        for cluster, cluster_cfg in config.get("clusters", {}).get(dataset, {}).items():
            cfgs[(dataset, cluster)] = cluster_cfg
    return cfgs


def predictions_generated(output_dir, dataset, cluster):
    path = os.path.join(
        output_dir, "uniformly_processed", "candidate_e2g_pairs", f"{dataset}_{cluster}_candidate_e2g_pairs.tsv.gz"
    )
    return os.path.exists(path)


def reformatted(output_dir, dataset, cluster):
    pattern = os.path.join(
        output_dir, "uniformly_processed", dataset, cluster, f"{dataset}_{cluster}_scE2G_*.e2g.tsv.gz"
    )
    return bool(glob.glob(pattern))


def derive_status(exclusion_reason, quality_pass, has_annotation, has_predictions, is_reformatted):
    if exclusion_reason in ("missing_qc_guide", "missing_per_cell_qc_table"):
        return "missing-qc-guide"
    if not quality_pass:
        return "excluded-quality"
    if not has_predictions:
        return "quality-pass-not-yet-processed"
    if not has_annotation:
        return "predictions-only-missing-annotation"
    if not is_reformatted:
        return "predictions-generated-reformat-pending"
    if has_predictions and is_reformatted and has_annotation and quality_pass:
        return "fully-processed"
    return "unknown"  # logic-gap signal, should never actually appear


def main():
    args = parse_args()
    out_path = args.out or os.path.join(args.output_dir, "report.tsv")

    cluster_stats = load_cluster_stats(args.output_dir)
    datasets = {dataset for dataset, _ in cluster_stats}
    cluster_configs = load_cluster_configs(args.configs_dir, datasets)
    state_conn = igvf_state.connect(args.state_db)

    rows = []
    status_counts = {}
    for (dataset, cluster), stats_row in sorted(cluster_stats.items()):
        cluster_cfg = cluster_configs.get((dataset, cluster))
        if cluster_cfg is None:
            # Config missing/dataset not loadable -- report what we can, skip
            # annotation/prediction lookups that need cluster_cfg.
            print(f"[generate_report] WARNING: no cluster_cfg for {dataset}/{cluster} -- annotation/manifest columns will be blank", file=sys.stderr)
            cluster_cfg = {}

        exclusion_reason = stats_row["reason"]
        quality_pass = exclusion_reason == "pass"
        lookup_dataset, lookup_cluster = annotation_lookup_key(dataset, cluster, cluster_cfg)
        has_annotation = igvf_state.get_cell_annotation(state_conn, lookup_dataset, lookup_cluster) is not None
        manifest_excluded = bool(cluster_cfg.get("igvf_manifest_excluded", False))
        has_predictions = predictions_generated(args.output_dir, dataset, cluster)
        is_reformatted = reformatted(args.output_dir, dataset, cluster)
        status = derive_status(exclusion_reason, quality_pass, has_annotation, has_predictions, is_reformatted)
        status_counts[status] = status_counts.get(status, 0) + 1

        rows.append({
            "dataset": dataset,
            "cluster": cluster,
            "cell_count": stats_row["cell_count"],
            "fragments_total": stats_row["fragments_total"],
            "umi_count": stats_row["umi_count"],
            "quality_pass": "y" if quality_pass else "n",
            "exclusion_reason": exclusion_reason,
            "has_cell_annotation": "y" if has_annotation else "n",
            "igvf_manifest_excluded": "y" if manifest_excluded else "n",
            "predictions_generated": "y" if has_predictions else "n",
            "reformatted": "y" if is_reformatted else "n",
            "status": status,
        })

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=REPORT_HEADER, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    state_conn.close()

    print(f"[generate_report] wrote {len(rows)} row(s) to {out_path}", file=sys.stderr)
    for status, count in sorted(status_counts.items(), key=lambda kv: -kv[1]):
        print(f"[generate_report]   {status}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
