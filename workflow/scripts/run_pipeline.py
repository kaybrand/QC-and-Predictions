#!/usr/bin/env python
"""Single entry point for one dataset, end to end.

    python workflow/scripts/run_pipeline.py \
        --pipeline-config workflow/config/igvf0_pipeline_config.yaml [-n]

Five stages, in order, stopping according to the failure policy below:

  [0] preflight  Local only. Load the config, resolve exclusions, compute the
                 cluster sets, write a plan TSV, and say exactly which stages
                 this mode will and won't run.
  [1] warm       DEFAULT MODE ONLY. Under an flock: fetch_if_stale (one wholesale
                 portal GET, 24h TTL) then derive_scopes (pure local) to build
                 this dataset's cell_annotations rows. Writes the CellAnnotation
                 snapshot Snakemake reads, plus a per-cluster status TSV.
  [2] snakemake  ONE invocation. --conda-prefix and --rerun-triggers mtime are
                 baked in because forgetting either has caused real incidents.
                 No --omit-from, ever: local_only makes the reformat targets
                 empty by construction instead.
  [3] manifest   DEFAULT MODE ONLY. manage_igvf_metadata.py in preview (or
                 validate) mode, scoped to MANIFEST_READY.
  [4] audit      DEFAULT MODE ONLY. Reads manifest_coverage.tsv, regenerates
                 report.tsv, and decides the exit code.

WHY A DRIVER AND NOT A SNAKEMAKE CHECKPOINT. The cache warm consumes nothing the
DAG produces -- it is a prerequisite, not an intermediate -- so modelling it
mid-DAG would add machinery for a dependency that doesn't exist, and `snakemake
-n` cannot execute a checkpoint, which would mean a dry run could never again
enumerate the reformat jobs. Manifest generation can't be a Snakemake rule at all:
its output filenames are round{N}_{table}_{variant}_{post,patch}.tsv, and
post-vs-patch comes from a live portal lookup, so the filename set is unknowable
at parse time. An external final stage is required regardless.

NEVER UPLOADS. --manifest-mode accepts only preview and validate. "upload" is
rejected with a pointer to the standalone command. Uploading to the IGVF Data
Portal is a human action, typed deliberately, every time.

EXIT CODES. 0 means complete and verified, and nothing else does.
  0  every manifest-eligible cluster resolved and every expected row is present
  1  the run produced output but something is incomplete (portal unreachable,
     a cluster's annotation unresolved, a manifest gap) -- always itemised
  2  preflight/config error: nothing ran

MULTI-DATASET. One invocation per dataset; they are independent and need no
ordering, because derive_scopes works from the raw cache regardless of whose
fetch populated it:

    for ds in igvf0 igvf1 igvf2; do sbatch resources/run_pipeline.sbatch "$ds"; done

The first driver to find a stale TTL does the single wholesale GET; the rest reuse
it. The flock keeps their state.db writes serialised, and no Snakemake worker ever
opens state.db (see cell_annotation_snapshot.py).
"""

import argparse
import csv
import fcntl
import os
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone

import yaml

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
sys.path.insert(0, SCRIPT_DIR)

import cell_annotation_snapshot as cas  # noqa: E402
from resolve_exclusions import resolve_exclusions, manifest_eligible_clusters  # noqa: E402
from cell_annotations import annotation_lookup_key  # noqa: E402
from igvf_metadata import orchestrator  # noqa: E402

# Never omit this: Snakemake folds the target env directory's absolute path into
# each conda env hash, so without it a byte-identical env computes a different
# hash, reads as "needs building from scratch", and concurrent invocations race
# and corrupt each other's downloads. That happened for real and cost hours.
DEFAULT_CONDA_PREFIX = (
    "/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/.snakemake/conda"
)
WARM_STATUS_NAME = "{dataset}_cell_annotation_status.tsv"
PLAN_NAME = "{dataset}_pipeline_plan.tsv"


def log(msg):
    print(f"[run_pipeline] {msg}", file=sys.stderr, flush=True)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pipeline-config", required=True)
    p.add_argument(
        "--mode", choices=["default", "local_only"], default=None,
        help="overrides the config's pipeline_mode",
    )
    p.add_argument(
        "--manifest-mode", choices=["preview", "validate"], default="preview",
        help="preview (default) writes TSVs only; validate additionally runs "
        "iu_register.py --dry-run for real schema checking. Uploading is NOT "
        "reachable from here -- run manage_igvf_metadata.py --mode upload yourself.",
    )
    p.add_argument("-n", "--dry-run", action="store_true",
                   help="preview every stage: no portal GET, no derive, snakemake -n, no manifest/audit")
    p.add_argument("--conda-prefix", default=None)
    p.add_argument("--jobs", default="4")
    p.add_argument("--snakemake-arg", action="append", default=[],
                   help="extra argument passed through to snakemake (repeatable)")
    p.add_argument("--skip-snakemake", action="store_true",
                   help="stages 0/1/3/4 only -- for regenerating a manifest against an existing results tree")
    p.add_argument("--sce2g-modules", choices=["true", "false"], default="true",
                   help="false skips importing scE2G per dataset: a fast portal-only pass. "
                   "Only valid once Phase 2 is complete (see common.smk).")
    return p.parse_args()


def stage_header(n, name, detail=""):
    log("")
    log(f"===== stage {n}: {name} {'-- ' + detail if detail else ''}")


@contextmanager
def state_db_lock(state_db, why):
    """Exclusive access to state.db for the duration of the block.

    Held by EVERY stage that touches state.db, not just the warm stage: stage 3
    writes to it too, because orchestrator.run() calls refresh_if_stale (and thus
    derive_scopes) at the top. Without this, one driver's stage-3 writes could
    overlap another driver's stage-1 writes, which is exactly what the lock
    exists to prevent.

    Deliberately NOT held across stage 2 -- that's hours of scE2G compute, and
    holding it there would serialise every dataset.

    Taken even on a dry run: state.connect() runs CREATE TABLE IF NOT EXISTS plus
    the column migrations and commits, so opening the DB at all is a write in
    principle, and the invariant is "one accessor at any instant", not "one
    writer".

    flock, not SQLite's busy timeout: state.db is a WAL database on Lustre, and
    WAL's shared-memory index is only supported when every connection is on one
    host. flock blocks in the kernel (no busy-wait) and is confirmed working on
    this mount.
    """
    lock_path = f"{state_db}.warmlock"
    fd = os.open(lock_path, os.O_CREAT | os.O_RDWR, 0o644)
    try:
        log(f"waiting for state.db lock ({why})")
        fcntl.flock(fd, fcntl.LOCK_EX)
        log(f"state.db lock acquired ({why})")
        yield
    finally:
        fcntl.flock(fd, fcntl.LOCK_UN)
        os.close(fd)
        log(f"state.db lock released ({why})")


# ---------------------------------------------------------------------------
# [0] preflight
# ---------------------------------------------------------------------------
def preflight(args):
    with open(args.pipeline_config) as f:
        config = yaml.safe_load(f)

    mode = args.mode or config.get("pipeline_mode", "default")
    if mode not in ("default", "local_only"):
        raise SystemExit(f"pipeline_mode must be 'default' or 'local_only', got {mode!r}")

    datasets = sorted(config["clusters"])
    if len(datasets) != 1:
        log(f"NOTE: config spans {len(datasets)} datasets ({datasets}). Real execution of "
            "macs2/predictions is only verified for one dataset per run -- see common.smk's "
            "KNOWN LIMITATION.")

    output_dir = config.get("output_dir", "./results")
    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(os.path.join(REPO_ROOT, output_dir))

    included, upload_eligible, excluded, stats = resolve_exclusions(config, config["data_dir"])
    manifest_eligible = manifest_eligible_clusters(config, upload_eligible)
    all_clusters = {(d, c) for d, clusters in config["clusters"].items() for c in clusters}

    log(f"mode                : {mode}")
    log(f"output_dir          : {output_dir}")
    log(f"configured clusters : {len(all_clusters)}")
    log(f"included (processed): {len(included)}")
    log(f"upload-eligible     : {len(upload_eligible)}")
    log(f"manifest-eligible   : {len(manifest_eligible)}"
        + (f"  (igvf_manifest_excluded: {sorted(upload_eligible - manifest_eligible)})"
           if manifest_eligible != upload_eligible else ""))
    log(f"excluded            : {len(excluded)}")
    if mode == "local_only":
        log("stages              : [0] preflight -> [2] snakemake. "
            "Warm/manifest/audit SKIPPED (no portal contact of any kind).")
    else:
        log("stages              : [0] preflight -> [1] warm -> [2] snakemake -> "
            f"[3] manifest ({args.manifest_mode}) -> [4] audit")

    plan_dir = os.path.join(output_dir, cas.SNAPSHOT_DIR_NAME)
    os.makedirs(plan_dir, exist_ok=True)
    for dataset in datasets:
        path = os.path.join(plan_dir, PLAN_NAME.format(dataset=dataset))
        with open(path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["dataset", "cluster", "included", "upload_eligible", "manifest_eligible",
                        "exclusion_reason", "cell_count", "fragments_total", "umi_count"])
            for (ds, cluster) in sorted(k for k in all_clusters if k[0] == dataset):
                row = stats[(ds, cluster)]
                w.writerow([ds, cluster,
                            "y" if (ds, cluster) in included else "n",
                            "y" if (ds, cluster) in upload_eligible else "n",
                            "y" if (ds, cluster) in manifest_eligible else "n",
                            row["reason"], row["cell_count"], row["fragments_total"], row["umi_count"]])
        log(f"wrote {path}")

    return {
        "config": config, "mode": mode, "output_dir": output_dir, "datasets": datasets,
        "included": included, "upload_eligible": upload_eligible,
        "manifest_eligible": manifest_eligible, "excluded": excluded,
        "all_clusters": all_clusters,
    }


# ---------------------------------------------------------------------------
# [1] warm
# ---------------------------------------------------------------------------
def warm(ctx, dry_run):
    """Returns (ok, statuses). ok=False means the portal was unreachable -- the
    caller continues to stage 2 regardless, because losing hours of scE2G compute
    to a portal hiccup is strictly worse than a degraded run."""
    config, output_dir = ctx["config"], ctx["output_dir"]
    state_db = config.get("igvf", {}).get("state_db_path")
    if not state_db:
        raise SystemExit("igvf.state_db_path is required in default mode")

    from igvf_metadata import cell_metadata, state

    cluster_configs = {
        (d, c): cfg for d, clusters in config["clusters"].items() for c, cfg in clusters.items()
    }
    cluster_keys = set(cluster_configs)
    digest = cas.cluster_set_digest(cluster_keys)

    ok = True
    with state_db_lock(state_db, "warm"):
        conn = state.connect(state_db)
        try:
            last_fetch = state.latest_cell_annotation_fetch(conn)
            n_primary = len(state.all_primary_pseudobulks(conn))
            log(f"raw primary-pseudobulk cache: {n_primary} row(s), last portal fetch {last_fetch}")

            if dry_run:
                log("dry run: NOT fetching and NOT deriving. Snapshotting whatever is already cached, "
                    "so the snakemake preview reflects the current cache rather than a hypothetical one.")
                statuses = [
                    {
                        "dataset": d, "cluster": c,
                        "resolved": state.get_cell_annotation(
                            conn, *annotation_lookup_key(d, c, cluster_configs[(d, c)])
                        ) is not None,
                        "reason": "dry-run: not derived", "cell_annotation": None,
                    }
                    for d, c in sorted(cluster_keys)
                ]
            else:
                try:
                    fetched = cell_metadata.fetch_if_stale(conn, _portal_reader(config))
                    log(f"portal fetch {'performed' if fetched else 'skipped (cache fresh)'}")
                except Exception as e:
                    # Credentials are never echoed; only the exception type/message,
                    # which igvf_utils does not populate with key material.
                    ok = False
                    log(f"PORTAL FETCH FAILED ({type(e).__name__}: {e})")
                    log("DEGRADED: continuing with the existing raw cache. Stages 3-4 will be skipped "
                        "and this run will exit non-zero.")
                statuses = cell_metadata.derive_scopes(conn, cluster_keys, cluster_configs)

            # Snapshot: what Snakemake will actually read. Written per dataset,
            # from state.db's rows, immediately before stage 2.
            for dataset in ctx["datasets"]:
                rows = [
                    dict(r) for r in state.all_cell_annotations(conn) if r["dataset"] == dataset
                ]
                path = cas.write_snapshot(
                    cas.snapshot_path(output_dir, dataset),
                    rows,
                    fetched_at=state.latest_cell_annotation_fetch(conn) or "",
                    digest=digest,
                )
                log(f"wrote snapshot {path} ({len(rows)} annotation row(s))")
        finally:
            conn.close()

    # Per-cluster status: the answer to "why isn't this cluster reformatted".
    for dataset in ctx["datasets"]:
        path = os.path.join(output_dir, cas.SNAPSHOT_DIR_NAME, WARM_STATUS_NAME.format(dataset=dataset))
        with open(path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["dataset", "cluster", "resolved", "reason", "cell_annotation"])
            for s in sorted((s for s in statuses if s["dataset"] == dataset),
                            key=lambda s: s["cluster"]):
                w.writerow([s["dataset"], s["cluster"], "y" if s["resolved"] else "n",
                            s["reason"], s["cell_annotation"] or ""])
        log(f"wrote {path}")

    unresolved = [s for s in statuses if not s["resolved"]
                  and (s["dataset"], s["cluster"]) in ctx["manifest_eligible"]]
    for s in unresolved:
        log(f"UNRESOLVED {s['dataset']}/{s['cluster']}: {s['reason']}")
    return ok, statuses


def _portal_reader(config):
    from igvf_metadata import portal_client

    return portal_client.PortalReader(igvf_mode=config.get("igvf", {}).get("mode", "prod"))


# ---------------------------------------------------------------------------
# [2] snakemake
# ---------------------------------------------------------------------------
def run_snakemake(ctx, args):
    conda_prefix = (
        args.conda_prefix
        or ctx["config"].get("snakemake_conda_prefix")
        or DEFAULT_CONDA_PREFIX
    )
    if not os.path.isdir(conda_prefix):
        raise SystemExit(
            f"--conda-prefix {conda_prefix} does not exist. Never run without it: Snakemake would "
            "rebuild every conda env against a different hash, and concurrent runs race."
        )
    cmd = [
        "snakemake",
        "-s", os.path.join(REPO_ROOT, "workflow", "Snakefile"),
        "--configfile", os.path.abspath(args.pipeline_config),
        "--use-conda",
        "--conda-prefix", conda_prefix,
        # Without this an explicit --conda-prefix reads as "software environment
        # changed" and queues a full rebuild of already-correct outputs.
        "--rerun-triggers", "mtime",
        "--keep-going",
        "--jobs", str(args.jobs),
        "--config", f"pipeline_mode={ctx['mode']}", f"sce2g_modules={args.sce2g_modules}",
    ]
    if args.dry_run:
        cmd.append("-n")
    else:
        cmd += ["--executor", "slurm", "--profile", "slurm.smk9", "-p"]
    cmd += args.snakemake_arg
    log("running: " + " ".join(cmd))
    return subprocess.run(cmd, cwd=REPO_ROOT).returncode


# ---------------------------------------------------------------------------
# [3] manifest
# ---------------------------------------------------------------------------
def run_manifest(ctx, args, manifest_dir):
    ready = sorted(ctx["manifest_ready"])
    if not ready:
        log("no manifest-ready clusters -- nothing to generate")
        return 0
    keys = ",".join(f"{d}/{c}" for d, c in ready)
    excluded_keys = ",".join(f"{d}/{c}" for d, c in sorted(ctx["excluded"]))
    cmd = [
        sys.executable,
        os.path.join(SCRIPT_DIR, "manage_igvf_metadata.py"),
        "--pipeline-config", os.path.abspath(args.pipeline_config),
        "--cluster-keys", keys,
        "--excluded-cluster-keys", excluded_keys,
        "--state-db", ctx["config"]["igvf"]["state_db_path"],
        "--manifest-dir", manifest_dir,
        "--mode", args.manifest_mode,
    ]
    log(f"running manage_igvf_metadata.py --mode {args.manifest_mode} for {len(ready)} cluster(s)")
    # orchestrator.run() calls refresh_if_stale at the top, so this subprocess
    # WRITES state.db. Hold the same lock the warm stage uses.
    with state_db_lock(ctx["config"]["igvf"]["state_db_path"], "manifest"):
        return subprocess.run(cmd, cwd=REPO_ROOT).returncode


# ---------------------------------------------------------------------------
# [4] audit
# ---------------------------------------------------------------------------
def audit(ctx, args, manifest_dir):
    """Returns a list of problem strings. Empty means this run is complete."""
    problems = []
    for dataset in ctx["datasets"]:
        path = os.path.join(manifest_dir, dataset, orchestrator.MANIFEST_COVERAGE_NAME)
        if not os.path.exists(path):
            problems.append(f"{dataset}: no {orchestrator.MANIFEST_COVERAGE_NAME} was written")
            continue
        with open(path) as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        gaps = [r for r in rows if r["outcome"] in orchestrator.MANIFEST_GAP_OUTCOMES]
        by_outcome = {}
        for r in gaps:
            by_outcome.setdefault(r["outcome"], []).append(r)
        log(f"{dataset}: {len(rows)} coverage row(s), {len(gaps)} gap(s)")
        for outcome, group in sorted(by_outcome.items()):
            for r in group[:10]:
                log(f"  GAP {outcome}: {r['cluster']} {r['table']}/{r['variant'] or '(default)'} :: {r['reason']}")
            if len(group) > 10:
                log(f"  ... and {len(group) - 10} more {outcome}")
            problems.append(f"{dataset}: {len(group)} {outcome} row(s)")

    for dataset, cluster in sorted(ctx["manifest_eligible"] - ctx["manifest_ready"]):
        problems.append(f"{dataset}/{cluster}: manifest-eligible but has no CellAnnotation")

    cmd = [
        sys.executable, os.path.join(SCRIPT_DIR, "generate_report.py"),
        "--output-dir", ctx["output_dir"],
        "--state-db", ctx["config"]["igvf"]["state_db_path"],
        "--configs-dir", os.path.join(REPO_ROOT, "workflow", "config"),
        "--manifest-dir", manifest_dir,
    ]
    if subprocess.run(cmd, cwd=REPO_ROOT).returncode != 0:
        problems.append("generate_report.py failed")
    return problems


def main():
    args = parse_args()
    started = datetime.now(timezone.utc).isoformat()

    stage_header(0, "preflight")
    ctx = preflight(args)
    mode = ctx["mode"]
    manifest_dir = ctx["config"].get("igvf", {}).get("manifest_dir") or os.path.join(
        ctx["output_dir"], "igvf_manifests"
    )

    degraded = []
    reformat_eligible = set()

    if mode == "default":
        stage_header(1, "warm the CellAnnotation cache")
        fetch_ok, statuses = warm(ctx, args.dry_run)
        if not fetch_ok:
            degraded.append("portal fetch failed; manifest stages skipped")
        resolved = {(s["dataset"], s["cluster"]) for s in statuses if s["resolved"]}
        reformat_eligible = ctx["upload_eligible"] & resolved
    else:
        stage_header(1, "warm", "SKIPPED (local_only: no portal contact)")
        fetch_ok = True

    ctx["manifest_ready"] = ctx["manifest_eligible"] & reformat_eligible
    if mode == "default":
        log(f"manifest-ready      : {len(ctx['manifest_ready'])} of "
            f"{len(ctx['manifest_eligible'])} manifest-eligible")

    rc_snakemake = 0
    snapshot_paths = [cas.snapshot_path(ctx["output_dir"], d) for d in ctx["datasets"]]
    try:
        if args.skip_snakemake:
            stage_header(2, "snakemake", "SKIPPED (--skip-snakemake)")
        else:
            stage_header(2, "snakemake", f"mode={mode}{' (dry run)' if args.dry_run else ''}")
            rc_snakemake = run_snakemake(ctx, args)
            if rc_snakemake != 0:
                degraded.append(f"snakemake exited {rc_snakemake}")
    finally:
        # The snapshot is a temp artifact: removing it means a later bare
        # `snakemake` in default mode aborts pointing at this driver rather than
        # silently dropping the reformat targets, and no stale values can ever
        # carry into a future run.
        #
        # Removed whenever the warm stage wrote one -- NOT conditioned on snakemake
        # having run. --skip-snakemake used to leave it behind, which is precisely
        # the state this deletion exists to prevent: a snapshot on disk that a later
        # bare `snakemake` would happily consume while still inside its freshness
        # window, with nobody having warmed the cache for that invocation.
        if mode == "default":
            for path in snapshot_paths:
                if cas.remove_snapshot(path):
                    log(f"removed temp snapshot {path}")

    problems = []
    if mode != "default":
        stage_header(3, "manifest", "SKIPPED (local_only)")
        stage_header(4, "audit", "SKIPPED (local_only)")
    elif args.dry_run:
        stage_header(3, "manifest", "SKIPPED (dry run)")
        stage_header(4, "audit", "SKIPPED (dry run)")
    elif not fetch_ok:
        stage_header(3, "manifest", "SKIPPED (portal fetch failed)")
        stage_header(4, "audit", "SKIPPED (portal fetch failed)")
    else:
        stage_header(3, "manifest", args.manifest_mode)
        rc = run_manifest(ctx, args, manifest_dir)
        if rc != 0:
            degraded.append(f"manage_igvf_metadata.py exited {rc}")
        stage_header(4, "audit")
        problems = audit(ctx, args, manifest_dir)

    log("")
    log(f"===== summary (started {started})")
    log(f"mode={mode} dry_run={args.dry_run} manifest_dir={manifest_dir}")
    for item in degraded:
        log(f"DEGRADED: {item}")
    for item in problems:
        log(f"PROBLEM : {item}")
    if degraded or problems:
        log("RESULT: INCOMPLETE -- see the itemised list above. Exit 1.")
        return 1
    if args.dry_run:
        log("RESULT: dry run completed. No portal GET, no derive, no rule executed.")
        return 0
    log("RESULT: complete and verified. Exit 0.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
