#!/usr/bin/env python
"""CLI entrypoint for the IGVF metadata uploader -- shared by the
pipeline-integrated Snakemake rule (one cluster at a time, low latency) and
the standalone scanner/backfill/reconciliation script (a full run's worth of
clusters at once). Both just assemble a different --cluster-keys list and
call into igvf_metadata.orchestrator.run; see that module's docstring for
the plan-then-execute mechanism and igvf_metadata/portal_client.py for why
real submission shells out to igvf_utils' own iu_register.py rather than
calling the API directly.

Runs alongside (not instead of) manage_synapse_manifest.py: separate
destination, separate state DB, same UPLOAD_ELIGIBLE_CLUSTERS feed from
resolve_exclusions.py.

--mode defaults to "preview": no contact with iu_register.py at all, just
the per-table post/patch TSVs written to --manifest-dir for review. Nothing
here ever uploads for real unless --mode upload is passed explicitly.
"""

import argparse
import sys

import yaml

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from igvf_metadata import orchestrator, tables  # noqa: F401 (tables import registers table modules)
from igvf_metadata.context import IgvfConfig
from pipeline_paths import resolve_repo_relative, repo_root_from_script


def parse_cluster_keys(raw):
    keys = set()
    for token in raw.split(","):
        if not token:
            continue
        dataset, _, cluster = token.partition("/")
        keys.add((dataset, cluster))
    return keys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pipeline-config", required=True, help="path to a *_pipeline_config.yaml")
    p.add_argument(
        "--cluster-keys", required=True, help='comma-separated "dataset/cluster" tokens, upload-eligible this run'
    )
    p.add_argument(
        "--excluded-cluster-keys",
        default="",
        help='comma-separated "dataset/cluster" tokens excluded this run (recorded, never uploaded)',
    )
    p.add_argument("--state-db", required=True)
    p.add_argument(
        "--manifest-dir",
        required=True,
        help="where per-dataset subfolders are written: <dataset>/round{N}_{table}[_{variant}]_"
        "{post,patch}.tsv for review/upload (ephemeral, shrinks as pieces go live), plus "
        "<dataset>/<object_type>.tsv, a durable accumulator of every alias ever confirmed live",
    )
    p.add_argument(
        "--mode",
        choices=["preview", "validate", "upload"],
        default="preview",
        help="preview (default): write TSVs only, no iu_register.py call. "
        "validate: also run iu_register.py --dry-run (real schema validation, zero writes). "
        "upload: run iu_register.py for real -- requires explicitly choosing this.",
    )
    p.add_argument("--tables", default="", help="comma-separated table names to restrict to (default: all registered)")
    p.add_argument(
        "--igvf-mode", default="prod",
        help="passed through to igvf_utils, e.g. sandbox/prod -- defaults to prod (the real production portal)",
    )
    p.add_argument("--iu-register-path", default=None, help="override the default iu_register.py path")
    args = p.parse_args()

    with open(args.pipeline_config) as f:
        config = yaml.safe_load(f)

    igvf_cfg = IgvfConfig.from_dict(config.get("igvf", {}))
    cluster_keys = parse_cluster_keys(args.cluster_keys)
    excluded = parse_cluster_keys(args.excluded_cluster_keys)
    table_names = [t for t in args.tables.split(",") if t] or None

    cluster_configs = {}
    for dataset, cluster in cluster_keys:
        cluster_configs[(dataset, cluster)] = config["clusters"][dataset][cluster]

    # Every (dataset, cluster) in the whole pipeline config, NOT just this invocation's
    # --cluster-keys -- the cell_metadata cache (built from one wholesale multireport GET
    # covering every primary pseudobulk on the portal) must cover every cluster we know
    # about, or a later invocation for a different --cluster-keys subset, within the same
    # 24h TTL, finds nothing cached for its own clusters and can do nothing about it.
    all_cluster_configs = {
        (dataset, cluster): cfg
        for dataset, clusters in config["clusters"].items()
        for cluster, cfg in clusters.items()
    }

    output_dir = resolve_repo_relative(config.get("output_dir", "./results"), repo_root_from_script(__file__))

    run_kwargs = dict(
        cluster_keys=cluster_keys,
        cluster_configs=cluster_configs,
        all_cluster_configs=all_cluster_configs,
        igvf_cfg=igvf_cfg,
        scE2G_dir=config["scE2G_dir"],
        data_dir=config["data_dir"],
        output_dir=output_dir,
        state_db_path=args.state_db,
        manifest_dir=args.manifest_dir,
        mode=args.mode,
        table_names=table_names,
        excluded=excluded,
        igvf_mode=args.igvf_mode,
    )
    if args.iu_register_path:
        run_kwargs["iu_register_path"] = args.iu_register_path

    report = orchestrator.run(**run_kwargs)
    print(f"[manage_igvf_metadata] mode={args.mode}", file=sys.stderr)
    for table_name, table_report in report.items():
        print(f"[manage_igvf_metadata] {table_name}: {table_report['counts']}", file=sys.stderr)


if __name__ == "__main__":
    main()
