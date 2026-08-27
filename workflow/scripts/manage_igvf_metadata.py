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


def parse_cluster_keys(raw, clusters=None, strict=False, flag="--cluster-keys"):
    """Parse comma-separated scope tokens into a set of (dataset, cluster).

    Two accepted token forms:
      "dataset/cluster"  -- one specific cluster
      "dataset"          -- EVERY cluster that dataset has in the config, which
                            needs `clusters` (the config's "clusters" mapping)

    The bare-dataset form exists because a plain `--cluster-keys igvf4` used to
    partition() into ("igvf4", "") and then die on `config["clusters"]["igvf4"][""]`
    with a bare `KeyError: ''` -- which says nothing about the real mistake. It is
    also the obvious way to ask for a whole dataset, so it now means that.

    strict=True (for --cluster-keys, whose keys must index the config) validates
    every token and exits with an actionable message listing what IS available.
    Left lenient for --excluded-cluster-keys, which is only recorded in the
    ledger and never used to look anything up.
    """
    clusters = clusters or {}
    keys = set()
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        dataset, sep, cluster = token.partition("/")

        if not sep or not cluster:
            # Bare dataset -> expand to all of its clusters.
            known = clusters.get(dataset)
            if known:
                keys.update((dataset, c) for c in known)
                continue
            if strict:
                sys.exit(
                    f"{flag}: {token!r} names no dataset in this config.\n"
                    f"  datasets available: {', '.join(sorted(clusters)) or '(none)'}\n"
                    f"  use 'dataset/cluster' for one cluster, or 'dataset' for all of them."
                )
            keys.add((dataset, cluster))
            continue

        if strict:
            if dataset not in clusters:
                sys.exit(
                    f"{flag}: unknown dataset {dataset!r} in token {token!r}.\n"
                    f"  datasets available: {', '.join(sorted(clusters)) or '(none)'}"
                )
            if cluster not in clusters[dataset]:
                available = sorted(clusters[dataset])
                shown = ", ".join(available[:12]) + (f", ... (+{len(available)-12} more)" if len(available) > 12 else "")
                sys.exit(
                    f"{flag}: unknown cluster {cluster!r} in dataset {dataset!r}.\n"
                    f"  {len(available)} cluster(s) available: {shown}"
                )
        keys.add((dataset, cluster))
    return keys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pipeline-config", required=True, help="path to a *_pipeline_config.yaml")
    p.add_argument(
        "--cluster-keys", required=True,
        help='comma-separated scope tokens, upload-eligible this run. "dataset/cluster" for one '
             'cluster; a bare "dataset" means every cluster that dataset has in the config. '
             'e.g. "igvf4/wtc11_macrophage_m0" or "igvf4"',
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
    all_clusters = config.get("clusters") or {}
    cluster_keys = parse_cluster_keys(args.cluster_keys, all_clusters, strict=True)
    excluded = parse_cluster_keys(
        args.excluded_cluster_keys, all_clusters, flag="--excluded-cluster-keys"
    )
    table_names = [t for t in args.tables.split(",") if t] or None

    if not cluster_keys:
        sys.exit("--cluster-keys resolved to no clusters; nothing to do.")
    scopes = ", ".join(f"{d}/{c}" for d, c in sorted(cluster_keys))
    print(
        f"[manage_igvf_metadata] mode={args.mode}  {len(cluster_keys)} cluster(s): "
        f"{scopes if len(cluster_keys) <= 12 else scopes[:400] + ' ...'}",
        file=sys.stderr,
    )

    # parse_cluster_keys(strict=True) already guaranteed both levels exist, so this
    # can no longer raise a bare KeyError on a mistyped or half-written token.
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

    # Deliberately no per-table count dump here. orchestrator.run already logs each
    # table's counts as it goes and then prints one SUMMARY block (see
    # orchestrator._print_run_summary) covering what was already there, what this
    # pass uploaded, and what is still pending. Repeating the same dicts underneath
    # that block -- which is what this used to do -- just doubled every line and
    # buried the summary.
    summary = report.get("_summary") or {}
    if args.mode == "upload":
        failed = summary.get("failed") or []
        if failed:
            print(
                f"[manage_igvf_metadata] {len(failed)} row(s) could not be verified after submission; "
                "see the FAILED section above",
                file=sys.stderr,
            )
            return 1
    return 0


if __name__ == "__main__":
    # main() returns non-zero when an upload pass submitted rows it could not read
    # back afterwards -- exit code carries that, so a wrapper script can react
    # without parsing the log.
    sys.exit(main())
