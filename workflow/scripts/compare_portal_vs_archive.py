#!/usr/bin/env python3
"""Compare portal-downloaded pseudobulk files against the existing archive.

Answers the two questions that matter for release: which pseudobulks can be
shared as they are, and which have to go back through human QC filtering.

    python compare_portal_vs_archive.py \
        --state-db resources/igvf_metadata_state.db \
        --download-root $SCRATCH/IGVF_Data_Portal_download_08-2026 \
        --archive-root /scratch/users/kaybrand/pseudobulk/igvf \
        --out-tsv results/portal_vs_archive.tsv

THE JOIN KEY IS THE DIRECTORY IDENTITY: (dataset, annotation, subsample).
NOT the subsample alone. One MULTI-seq subsample contributes to MANY clusters
(igvf13's IGVFSM2403HODV appears under both "Lymphocyte" and "SMC"), so a
subsample-keyed map silently keeps only the last one seen -- the same trap
cell_metadata.py's docstring records, and one I walked into once while planning
this: it produced impossible symmetric "renames" (igvf10 jurkat ->
jurkat_pma_cd3_4hr AND the exact reverse) before I noticed the dict was
overwriting.

WHY A NAIVE md5 COMPARISON IS WRONG HERE
----------------------------------------
The portal and the archive store the same content with different compression.
Measured on a real downloaded file:

    raw md5, portal  per_cell_qc.tsv.gz   1e622f770a270793273f58e12b82d8ba
    raw md5, archive per_cell_qc.tsv      b78fcbc5e932e844a719647bcf0b46b3
    decompressed, both                    b78fcbc5e932e844a719647bcf0b46b3

Identical content. Comparing raw md5s would report all 1028 per-cell quality
reports as "changed". So every comparison falls back to a CANONICAL CONTENT
hash -- the md5 of the decompressed bytes -- whenever the raw hashes differ.
gzip also embeds an mtime, so even two .gz files of identical content can differ
byte for byte; the same fallback covers that.

FALSE "UNCHANGED" IS THE EXPENSIVE ERROR
----------------------------------------
Per the task spec, wrongly calling a changed file "unchanged" is worse than the
reverse, because it means shipping stale data as current. So:
  - "unchanged" is only ever claimed on a positive hash match, never inferred;
  - anything that differs and cannot be positively explained is reported
    MEANINGFUL or UNCERTAIN, never quietly absorbed;
  - h5ad files whose bytes differ need h5py to introspect. Without it the
    verdict is UNCERTAIN, not "probably fine". Run under an interpreter with
    h5py for full fidelity (the filter_multiome env has h5py + anndata;
    igvf_utils_env has neither).
"""

import argparse
import csv
import gzip
import hashlib
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict

CHUNK = 1 << 20
DIRNAME_RE = re.compile(r"^annotation-(?P<annotation>.+)-(?P<subsample>IGVFSM\w+)$")

# Verdicts, ordered from best to worst news.
IDENTICAL = "identical"            # same bytes
COMPRESSION_ONLY = "compression_only"  # same content, different container
MEANINGFUL = "meaningful"          # content genuinely differs
UNCERTAIN = "uncertain"            # differs, could not characterise -- treat as meaningful


def log(msg):
    print(f"[compare] {msg}", file=sys.stderr)


def _open_maybe_gz(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def raw_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def content_md5(path):
    """md5 of the DEcompressed bytes -- the canonical form, so a .tsv.gz and a
    .tsv of the same content hash equal. Streamed: fragments files reach several
    GB uncompressed and must never be read whole."""
    h = hashlib.md5()
    with _open_maybe_gz(path) as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def counterpart(archive_root, dataset, dirname, filename):
    """The archive path for a portal file, tolerating the compression mismatch:
    the portal ships per_cell_qc.tsv.gz where the archive holds per_cell_qc.tsv
    (and likewise pseudobulk_expression / peaks). Returns None if neither form
    is present."""
    base = os.path.join(archive_root, dataset, "pseudobulks", dirname)
    candidates = [filename]
    if filename.endswith(".gz"):
        candidates.append(filename[:-3])
    else:
        candidates.append(filename + ".gz")
    for cand in candidates:
        path = os.path.join(base, cand)
        if os.path.exists(path):
            return path
    return None


def _tsv_detail(portal_path, archive_path, barcode_col=None):
    """Characterise a difference between two tab-separated files: header change,
    row-count change, the first differing line, and (when a barcode column is
    given) how the barcode sets differ. One streaming pass each."""
    detail = {}
    with _open_maybe_gz(portal_path) as a, _open_maybe_gz(archive_path) as b:
        ha = a.readline().decode("utf-8", "replace").rstrip("\n")
        hb = b.readline().decode("utf-8", "replace").rstrip("\n")
        detail["header_changed"] = ha != hb
        if ha != hb:
            detail["header_portal"] = ha[:300]
            detail["header_archive"] = hb[:300]
        na = nb = 0
        first_diff = None
        bar_a, bar_b = set(), set()
        idx = None
        if barcode_col and not detail["header_changed"]:
            cols = ha.split("\t")
            idx = cols.index(barcode_col) if barcode_col in cols else None
        while True:
            la, lb = a.readline(), b.readline()
            if not la and not lb:
                break
            if la:
                na += 1
                if idx is not None:
                    f = la.decode("utf-8", "replace").rstrip("\n").split("\t")
                    if len(f) > idx:
                        bar_a.add(f[idx])
            if lb:
                nb += 1
                if idx is not None:
                    f = lb.decode("utf-8", "replace").rstrip("\n").split("\t")
                    if len(f) > idx:
                        bar_b.add(f[idx])
            if first_diff is None and la != lb:
                first_diff = max(na, nb)
    detail["rows_portal"] = na
    detail["rows_archive"] = nb
    detail["first_differing_line"] = first_diff
    if idx is not None:
        detail["barcodes_portal"] = len(bar_a)
        detail["barcodes_archive"] = len(bar_b)
        detail["barcodes_only_portal"] = len(bar_a - bar_b)
        detail["barcodes_only_archive"] = len(bar_b - bar_a)
    return detail


def _h5ad_detail(portal_path, archive_path):
    """Compare two .h5ad files structurally. HDF5 embeds creation metadata, so
    identical data can differ byte for byte; this decides whether the DATA is
    the same. Returns (verdict, detail)."""
    try:
        import h5py
        import numpy as np
    except ImportError:
        return UNCERTAIN, {"note": "h5py unavailable -- cannot introspect h5ad; run under filter_multiome env"}

    def summarise(path):
        out = {}
        with h5py.File(path, "r") as f:
            for key in ("obs", "var"):
                if key in f:
                    idx_name = f[key].attrs.get("_index")
                    if isinstance(idx_name, bytes):
                        idx_name = idx_name.decode()
                    if idx_name and idx_name in f[key]:
                        vals = f[key][idx_name][:]
                        out[f"n_{key}"] = len(vals)
                        h = hashlib.md5()
                        for v in vals:
                            h.update(v if isinstance(v, bytes) else str(v).encode())
                        out[f"{key}_names_md5"] = h.hexdigest()
            if "X" in f:
                x = f["X"]
                if isinstance(x, h5py.Group):  # sparse
                    out["X_kind"] = "sparse"
                    for part in ("data", "indices", "indptr"):
                        if part in x:
                            h = hashlib.md5()
                            arr = x[part]
                            for i in range(0, arr.shape[0], 1 << 20):
                                h.update(np.ascontiguousarray(arr[i : i + (1 << 20)]).tobytes())
                            out[f"X_{part}_md5"] = h.hexdigest()
                else:
                    out["X_kind"] = "dense"
                    out["X_shape"] = str(x.shape)
                    h = hashlib.md5()
                    for i in range(0, x.shape[0], 4096):
                        h.update(np.ascontiguousarray(x[i : i + 4096]).tobytes())
                    out["X_md5"] = h.hexdigest()
        return out

    try:
        pa, ar = summarise(portal_path), summarise(archive_path)
    except Exception as exc:  # noqa: BLE001
        return UNCERTAIN, {"note": f"h5ad introspection failed: {type(exc).__name__}: {exc}"}
    diffs = {k: f"portal={pa.get(k)!r} archive={ar.get(k)!r}" for k in set(pa) | set(ar) if pa.get(k) != ar.get(k)}
    if not diffs:
        return COMPRESSION_ONLY, {"note": "HDF5 bytes differ but obs/var names and X payload are identical"}
    return MEANINGFUL, {"h5ad_differences": "; ".join(f"{k}: {v}" for k, v in sorted(diffs.items()))}


def compare_file(portal_path, archive_path, content_type, characterise=True):
    """Returns (status, verdict, detail_dict). status is unchanged/changed."""
    r_portal, r_archive = raw_md5(portal_path), raw_md5(archive_path)
    if r_portal == r_archive:
        return "unchanged", IDENTICAL, {"raw_md5": r_portal}

    c_portal, c_archive = content_md5(portal_path), content_md5(archive_path)
    if c_portal == c_archive:
        return (
            "unchanged",
            COMPRESSION_ONLY,
            {
                "raw_md5_portal": r_portal,
                "raw_md5_archive": r_archive,
                "content_md5": c_portal,
                "note": "decompressed content identical; only the compression container differs",
            },
        )

    detail = {"raw_md5_portal": r_portal, "raw_md5_archive": r_archive,
              "content_md5_portal": c_portal, "content_md5_archive": c_archive}
    if not characterise:
        return "changed", UNCERTAIN, {**detail, "note": "characterisation skipped (--no-characterize)"}

    if content_type == "cell by gene matrix":
        verdict, extra = _h5ad_detail(portal_path, archive_path)
        # A structural match here means the payload is equal, so it is NOT a
        # content change -- but the file is still not byte-identical.
        status = "unchanged" if verdict == COMPRESSION_ONLY else "changed"
        return status, verdict, {**detail, **extra}

    barcode_col = "barcode" if content_type == "per-cell quality report" else None
    try:
        extra = _tsv_detail(portal_path, archive_path, barcode_col=barcode_col)
    except Exception as exc:  # noqa: BLE001
        return "changed", UNCERTAIN, {**detail, "note": f"tsv characterisation failed: {exc}"}
    return "changed", MEANINGFUL, {**detail, **extra}


def archive_directories(archive_root, datasets):
    """Every annotation-*-IGVFSM* directory in the archive, as
    {(dataset, dirname)} -- for spotting archive dirs the portal no longer has."""
    found = set()
    for dataset in datasets:
        base = os.path.join(archive_root, dataset, "pseudobulks")
        if not os.path.isdir(base):
            continue
        for entry in os.listdir(base):
            if DIRNAME_RE.match(entry):
                found.add((dataset, entry))
    return found


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--state-db", required=True)
    p.add_argument("--download-root", required=True)
    p.add_argument("--archive-root", required=True)
    p.add_argument("--out-tsv", required=True, help="per-file comparison table")
    p.add_argument("--out-cluster-tsv", default=None, help="per-(dataset,annotation) rollup (default: alongside --out-tsv)")
    p.add_argument("--datasets", default="", help="comma-separated datasets to restrict to")
    p.add_argument("--no-characterize", action="store_true", help="skip deep diffing of changed files")
    args = p.parse_args(argv)

    only = {d.strip() for d in args.datasets.split(",") if d.strip()} or None

    conn = sqlite3.connect(args.state_db)
    conn.row_factory = sqlite3.Row
    rows = [dict(r) for r in conn.execute("SELECT * FROM portal_files").fetchall()]
    conn.close()
    if only:
        rows = [r for r in rows if r["dataset"] in only]
    log(f"{len(rows)} ledger row(s) in scope")

    results = []
    counts = Counter()
    for r in rows:
        dataset, dirname = r["dataset"], None
        if r["local_path"]:
            dirname = os.path.basename(os.path.dirname(r["local_path"]))
        filename = os.path.basename(r["local_path"]) if r["local_path"] else None
        rec = {
            "dataset": dataset,
            "annotation": r["annotation"],
            "subsample": r["subsample"],
            "dirname": dirname,
            "content_type": r["content_type"],
            "accession": r["accession"],
            "portal_path": r["local_path"],
            "archive_path": "",
            "status": "",
            "verdict": "",
            "portal_md5_reported": r["md5sum"],
            "detail": "",
            "review_reasons": r["review_reasons"] or "",
        }

        if r["download_state"] != "done" or not r["local_path"] or not os.path.exists(r["local_path"]):
            rec["status"] = "not_downloaded"
            rec["verdict"] = ""
            rec["detail"] = f"download_state={r['download_state']}"
            counts["not_downloaded"] += 1
            results.append(rec)
            continue

        arch = counterpart(args.archive_root, dataset, dirname, filename) if dataset and dirname else None
        if not arch:
            rec["status"] = "novel"
            rec["verdict"] = ""
            rec["detail"] = "no counterpart in archive"
            counts["novel"] += 1
            results.append(rec)
            continue

        status, verdict, detail = compare_file(
            r["local_path"], arch, r["content_type"], characterise=not args.no_characterize
        )
        rec["archive_path"] = arch
        rec["status"] = status
        rec["verdict"] = verdict
        rec["detail"] = "; ".join(f"{k}={v}" for k, v in detail.items())
        counts[f"{status}/{verdict}"] += 1
        results.append(rec)

    # Per-file table
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)) or ".", exist_ok=True)
    cols = ["dataset", "annotation", "subsample", "dirname", "content_type", "accession",
            "status", "verdict", "portal_md5_reported", "portal_path", "archive_path",
            "review_reasons", "detail"]
    tmp = args.out_tsv + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for rec in results:
            w.writerow(rec)
    os.replace(tmp, args.out_tsv)
    log(f"wrote {args.out_tsv} ({len(results)} rows)")

    # Per-(dataset, annotation) rollup: the actual release decision.
    by_dir = defaultdict(list)
    for rec in results:
        by_dir[(rec["dataset"], rec["annotation"] or rec["dirname"])].append(rec)
    cluster_rows = []
    for (dataset, annotation), recs in sorted(by_dir.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]))):
        statuses = {rec["status"] for rec in recs}
        verdicts = {rec["verdict"] for rec in recs if rec["verdict"]}
        if "not_downloaded" in statuses:
            decision = "INCOMPLETE_DOWNLOAD"
        elif statuses == {"novel"}:
            decision = "NOVEL_needs_QC_filtering"
        elif "novel" in statuses:
            decision = "PARTIALLY_NOVEL_needs_QC_filtering"
        elif "changed" in statuses:
            decision = "CHANGED_needs_QC_filtering"
        elif statuses == {"unchanged"}:
            decision = "SHAREABLE_unchanged"
        else:
            decision = "REVIEW"
        cluster_rows.append(
            {
                "dataset": dataset,
                "annotation": annotation,
                "n_files": len(recs),
                "n_subsamples": len({rec["subsample"] for rec in recs}),
                "statuses": ",".join(sorted(s for s in statuses if s)),
                "verdicts": ",".join(sorted(verdicts)),
                "decision": decision,
            }
        )
    out_cluster = args.out_cluster_tsv or args.out_tsv.replace(".tsv", "") + "_by_cluster.tsv"
    tmp = out_cluster + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["dataset", "annotation", "n_files", "n_subsamples", "statuses", "verdicts", "decision"],
            delimiter="\t", lineterminator="\n",
        )
        w.writeheader()
        w.writerows(cluster_rows)
    os.replace(tmp, out_cluster)
    log(f"wrote {out_cluster} ({len(cluster_rows)} rows)")

    # Archive directories the portal no longer offers.
    datasets = {r["dataset"] for r in rows if r["dataset"]}
    portal_dirs = {(rec["dataset"], rec["dirname"]) for rec in results if rec["dirname"]}
    orphaned = sorted(archive_directories(args.archive_root, datasets) - portal_dirs)

    log("")
    log("================ SUMMARY ================")
    for k in sorted(counts):
        log(f"  {counts[k]:>6}  {k}")
    log("")
    log("  per-annotation decisions:")
    for decision, n in Counter(c["decision"] for c in cluster_rows).most_common():
        log(f"    {n:>5}  {decision}")
    if orphaned:
        log("")
        log(f"  archive directories with NO portal counterpart: {len(orphaned)}")
        for dataset, dirname in orphaned[:15]:
            log(f"    {dataset}/{dirname}")
        if len(orphaned) > 15:
            log(f"    ... and {len(orphaned) - 15} more")
        log("    (their derived products are superseded -- see the rerun plan)")
    changed = [r for r in results if r["status"] == "changed"]
    if changed:
        log("")
        log(f"  CHANGED files needing a human look: {len(changed)}")
        for rec in changed[:20]:
            log(f"    [{rec['verdict']}] {rec['dataset']}/{rec['dirname']}/{rec['content_type']}")
            log(f"        {rec['detail'][:400]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
