#!/usr/bin/env python3
"""Download pseudobulk component files from the IGVF Portal.

Discovery lives in igvf_metadata/portal_files.py, the transfer loop in
igvf_metadata/downloader.py, the ledger in igvf_metadata/state.py. This file is
just the CLI that wires them together.

TWO FAMILIES, selected by --family (they select disjoint sets on the portal and
cannot be combined in one run):

  --family kundaje (default)  primary-pseudobulk ATAC/RNA/per-cell-QC files.

    python download_portal_pseudobulks.py \
        --download-root /scratch/users/$USER/pseudobulk_portal \
        --state-db resources/igvf_metadata_state.db \
        --datasets igvf10 \
        --dry-run

  --family catlas             the Yang Li lab CATlas curated-set pseudobulks,
                              fragments only. DEFAULT IS ALL 627 CLUSTERS
                              (both species); --species narrows it.

    # every Yang Li cluster, both species
    python download_portal_pseudobulks.py --family catlas \
        --download-root /scratch/users/$USER/catlas --state-db resources/igvf_metadata_state.db

    # ONLY human (258 clusters, ~0.76 TiB)
    python download_portal_pseudobulks.py --family catlas --species human \
        --download-root /scratch/users/$USER/catlas --state-db resources/igvf_metadata_state.db

Run it under Slurm, never on a login node -- a full run is hundreds to thousands
of files and hours of network I/O. See resources/run_portal_download.sbatch.

OUTPUT LAYOUT, per family. Both put `dataset` first, which is what keeps species
from ever mixing: a CATlas dataset label is "catlas-{species}", so one directory
can only hold one species.

  kundaje: {download-root}/{dataset}/pseudobulks/annotation-{annotation}-{subsample}/
           This MIRRORS what the pipeline already expects, so repointing
           `pseudobulks_root` in a *_pipeline_config.yaml is the only change
           needed to run on downloaded rather than archive-extracted data -- no
           change to filter_atac_fragments.py or filter_rna_counts.py, which
           glob that layout.
  catlas:  {download-root}/catlas-{species}/fragments/{species}-{cluster}_fragments.tsv.gz
           There is no annotation/subsample level: a CATlas set is one cluster
           with one fragments file. These are consumed by exact path via each
           cluster's `atac_frag_file` config key, not by globbing.

Resumability is a property of the ledger, not of this process: every file's
outcome is committed as it completes, so an interrupted run (walltime, node
failure, Ctrl-C) is resumed by re-invoking with the same arguments. Files whose
portal md5sum still matches what was observed at download time are skipped;
files upstream has since replaced are re-fetched (see downloader.needs_download).

$SCRATCH note: the 90-day purge timer is reset only by real content writes. A
resumed run that SKIPS an unchanged file does not refresh that file's timer, so
a long-lived tree still needs periodic real activity to survive.
"""

import argparse
import os
import sys
import threading
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from igvf_metadata import downloader, portal_client, portal_files, state  # noqa: E402


def log(msg):
    print(f"[download_portal_pseudobulks] {msg}", file=sys.stderr)


def _now():
    return datetime.now(timezone.utc).isoformat()


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--download-root",
        required=True,
        help="destination root, e.g. /scratch/users/$USER/pseudobulk_portal. NOT defaulted on "
        "purpose: this writes many GB and must be a path you chose deliberately.",
    )
    p.add_argument("--state-db", required=True, help="the same SQLite ledger manage_igvf_metadata.py uses")
    p.add_argument(
        "--family",
        choices=("kundaje", "catlas"),
        default="kundaje",
        help="which pseudobulk family to fetch. 'kundaje' (default) = the primary "
        "ATAC/RNA/per-cell-QC pseudobulks this pipeline has always consumed. 'catlas' = the "
        "Yang Li lab curated-set CATlas pseudobulks (fragments only). The two select disjoint "
        "sets and cannot be combined in one run.",
    )
    p.add_argument(
        "--species",
        default="",
        help="CATlas only (--family catlas). Comma-separated species to restrict to, e.g. "
        "'human'. DEFAULT (empty) retrieves ALL Yang Li clusters, both species. Species can "
        f"never mix within a dataset regardless: the dataset label is 'catlas-{{species}}'. "
        f"Known species: {','.join(portal_files.CATLAS_SPECIES)}.",
    )
    p.add_argument(
        "--datasets",
        default="",
        help="comma-separated dataset labels to restrict to (e.g. igvf10,igvf9). Default: every "
        "dataset discovered on the portal. --family kundaje only; use --species for CATlas, "
        "whose dataset labels are derived from the species.",
    )
    p.add_argument(
        "--lab",
        default=portal_files.DEFAULT_LAB,
        help=f"only download primary pseudobulks from this lab @id (default {portal_files.DEFAULT_LAB}). "
        "Pass an empty string to disable the filter -- but note other labs' pseudobulks do not "
        "follow this directory convention. --family kundaje only; --family catlas is pinned to "
        f"{portal_files.CATLAS_LAB} by construction.",
    )
    p.add_argument(
        "--content-types",
        default=None,
        help="comma-separated portal content_type values to fetch. Defaults to every target "
        "content_type of the selected --family.",
    )
    p.add_argument("--max-workers", type=int, default=4, help="concurrent transfers (default 4)")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="discover and update the ledger, but transfer nothing. Safe to run anywhere.",
    )
    p.add_argument("--igvf-mode", default="prod", help="passed to igvf_utils (default prod)")
    p.add_argument(
        "--verify-archive",
        default=None,
        help="optional: an existing archive root (e.g. /scratch/users/kaybrand/pseudobulk/igvf). "
        "Reports how many resolved paths already exist there -- a free check that path "
        "resolution is correct, for datasets that should already be present.",
    )
    return p.parse_args(argv)


def record_discovery(conn, records, download_root, now):
    """Persist what the portal says, and compute each file's local path. Returns
    the records annotated with local_path."""
    for rec in records:
        rel = rec.get("rel_path")
        rec["local_path"] = os.path.join(download_root, rel) if rel else None
        reasons = rec.get("review_reasons") or []
        state.upsert_portal_file(
            conn,
            rec["accession"],
            now,
            file_set=rec["file_set"],
            principal_analysis_set=rec["principal_analysis_set"],
            lab=rec["lab"],
            content_type=rec["content_type"],
            file_format=rec["file_format"],
            href=rec["href"],
            s3_uri=rec["s3_uri"],
            alias=rec["alias"],
            md5sum=rec["md5sum"],
            file_size=rec["file_size"],
            portal_status=rec["portal_status"],
            upload_status=rec["upload_status"],
            submitted_file_name=rec["submitted_file_name"],
            dataset=rec["dataset"],
            annotation=rec["annotation"],
            subsample=rec["subsample"],
            local_path=rec["local_path"],
            review_reasons=",".join(reasons) or None,
        )
    return records


def verify_against_archive(records, archive_root, skip_datasets=()):
    """How many resolved paths already exist in an existing archive. Accepts a
    decompressed counterpart too: the portal ships per_cell_qc.tsv.gz while the
    archive holds per_cell_qc.tsv, so a strict existence check would report
    every per-cell-QC file as absent."""
    present = absent = skipped = 0
    missing = []
    for rec in records:
        rel = rec.get("rel_path")
        if not rel:
            continue
        if rec.get("dataset") in skip_datasets:
            skipped += 1
            continue
        path = os.path.join(archive_root, rel)
        if os.path.exists(path) or (path.endswith(".gz") and os.path.exists(path[:-3])):
            present += 1
        else:
            absent += 1
            missing.append(rel)
    return present, absent, skipped, missing


def main(argv=None):
    args = parse_args(argv)
    datasets = {d.strip() for d in args.datasets.split(",") if d.strip()} or None
    species = {s.strip() for s in args.species.split(",") if s.strip()} or None
    content_types = (
        {c.strip() for c in args.content_types.split(",") if c.strip()}
        if args.content_types
        else None
    )

    # Reject flags belonging to the other family up front rather than ignoring
    # them: a run that silently drops "--species human" would download all 627
    # CATlas clusters, and a run that silently drops "--datasets" would download
    # far more than asked. Both are expensive mistakes to discover afterwards.
    if args.family == "catlas":
        if datasets:
            raise SystemExit(
                "--datasets does not apply to --family catlas: CATlas dataset labels are derived "
                "from the species (catlas-human/catlas-mouse). Use --species instead."
            )
        if args.lab != portal_files.DEFAULT_LAB:
            raise SystemExit(
                f"--lab does not apply to --family catlas, which is pinned to {portal_files.CATLAS_LAB}."
            )
        unknown_species = (species or set()) - set(portal_files.CATLAS_SPECIES)
        if unknown_species:
            # Checked here as well as in discover_catlas so a typo is a one-line
            # CLI error, not a traceback -- and so it fails before the Portal GET.
            raise SystemExit(
                f"unknown --species {sorted(unknown_species)}; "
                f"known: {','.join(portal_files.CATLAS_SPECIES)}"
            )
    elif species:
        raise SystemExit("--species only applies to --family catlas.")

    conn = state.connect(args.state_db)
    reader = portal_client.PortalReader(igvf_mode=args.igvf_mode)

    if args.family == "catlas":
        log(
            "family=catlas, species=%s" % (",".join(sorted(species)) if species else "ALL (default)")
        )
        records, report = portal_files.discover_catlas(
            reader, species=species, content_types=content_types
        )
    else:
        records, report = portal_files.discover(
            reader, lab=args.lab or None, content_types=content_types, datasets=datasets
        )
    portal_files.log_report(report)

    now = _now()
    records = record_discovery(conn, records, args.download_root, now)
    log(f"ledger updated for {len(records)} file(s)")

    if args.verify_archive:
        # igvf9/igvf13 are legitimately absent from the archive (see the plan's
        # comparison table), so they'd be false alarms in this particular check.
        present, absent, skipped, missing = verify_against_archive(records, args.verify_archive)
        log(f"archive check vs {args.verify_archive}: {present} present, {absent} absent, {skipped} skipped")
        for m in missing[:20]:
            log(f"  absent: {m}")
        if len(missing) > 20:
            log(f"  ... and {len(missing) - 20} more")

    unresolved = [r for r in records if not r["local_path"]]
    for rec in unresolved:
        state.record_download_result(
            conn,
            rec["accession"],
            "needs_review",
            _now(),
            error="could not resolve a local path: " + ",".join(rec.get("review_reasons") or []),
            bump_attempt=False,
        )
    if unresolved:
        log(f"{len(unresolved)} file(s) have no resolvable local path -> needs_review, not downloaded")

    todo = []
    skipped_unchanged = 0
    for rec in records:
        if not rec["local_path"]:
            continue
        row = state.get_portal_file(conn, rec["accession"])
        should, reason = downloader.needs_download(row, rec["local_path"])
        if should:
            todo.append((rec, reason))
        else:
            skipped_unchanged += 1
    log(f"{len(todo)} file(s) to download; {skipped_unchanged} already present and unchanged")

    if args.dry_run:
        log("--dry-run: no transfers performed")
        by_reason = Counter(reason.split(":")[0].split(" ")[0] for _, reason in todo)
        log(f"  would download, by reason: {dict(by_reason)}")
        conn.close()
        return 0

    if not todo:
        log("nothing to do")
        conn.close()
        return 0

    auth = reader.auth
    if not auth:
        log(
            "ERROR: no IGVF credentials in the environment (IGVF_API_KEY / IGVF_SECRET_KEY). "
            "igvf_utils falls back to anonymous access, which would 403 on every file. "
            "Export both and retry; never pass them on the command line."
        )
        conn.close()
        return 2

    base_url = reader.base_url
    total_bytes = sum(r["file_size"] or 0 for r, _ in todo)
    log(f"downloading {len(todo)} file(s), ~{total_bytes / 1e9:.2f} GB, {args.max_workers} worker(s)")

    # One Session per worker thread. Sessions are not documented as thread-safe,
    # and this sidesteps the question entirely.
    local = threading.local()

    def session_for_thread():
        if not hasattr(local, "session"):
            local.session = downloader.build_session(auth, pool_maxsize=2)
        return local.session

    def fetch(rec):
        return rec, downloader.download(
            session_for_thread(),
            base_url,
            rec["href"],
            rec["local_path"],
            expected_md5=rec["md5sum"],
            expected_size=rec["file_size"],
        )

    outcomes = Counter()
    done_bytes = 0
    # Results are recorded HERE, in the main thread: the sqlite3 connection is
    # bound to its creating thread (check_same_thread), and committing per file
    # as it lands is what makes an interrupted run resumable.
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        futures = {pool.submit(fetch, rec): rec for rec, _ in todo}
        for i, fut in enumerate(as_completed(futures), start=1):
            rec = futures[fut]
            try:
                rec, result = fut.result()
            except Exception as exc:  # noqa: BLE001 -- a worker crash must not lose the run
                state.record_download_result(
                    conn, rec["accession"], "failed", _now(), error=f"worker crashed: {exc!r}"
                )
                outcomes["failed"] += 1
                log(f"  CRASH {rec['accession']}: {exc!r}")
                continue
            state.record_download_result(
                conn,
                rec["accession"],
                result.state,
                _now(),
                bytes_written=result.bytes_written,
                md5_observed=result.md5_observed,
                error=result.error,
            )
            outcomes[result.state] += 1
            done_bytes += result.bytes_written or 0
            if result.state != "done":
                log(f"  {result.state.upper()} {rec['accession']} ({rec['content_type']}): {result.error}")
            elif not result.verified:
                log(f"  UNVERIFIED {rec['accession']}: portal reported no md5sum")
            if i % 50 == 0 or i == len(futures):
                log(f"  progress {i}/{len(futures)} files, {done_bytes / 1e9:.2f} GB")

    log(f"outcomes: {dict(outcomes)}")
    log(f"ledger state totals: {state.portal_file_state_counts(conn)}")
    conn.close()

    bad = outcomes["failed"] + outcomes["md5_mismatch"]
    if bad:
        log(f"{bad} file(s) did not complete cleanly -- re-run to retry; md5 mismatches kept as .part")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
