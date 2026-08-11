#!/usr/bin/env python3
"""Download IGVF fragments files listed in a pseudobulk report TSV (e.g.
washu_pseudobulk_report.tsv, produced by get_washu_pseudobulk_report.py) and
verify each download against its portal-reported md5sum.

Built to scale from the 10-file WashU test set up to ~500 pseudobulks: it
streams each download (never holds a whole file in memory), skips any file
already downloaded with a matching md5 (safe to re-run/resume a long batch),
and never fails the whole run over one bad file -- failures are collected and
reported at the end, and the process exits non-zero only if any remain.

Auth mirrors get_igvf_cell_annotation_report.py / get_washu_pseudobulk_report.py:
igvf_utils.Connection reads IGVF_API_KEY/IGVF_SECRET_KEY from the environment,
never from arguments or a file.
"""

import argparse
import csv
import hashlib
import os
import sys
from urllib.parse import urljoin

REQUIRED_COLUMNS = ["FragmentsFile_Href", "FragmentsFile_MD5", "FragmentsFile_Aliases"]


def log(msg):
    print(f"[download_and_verify_fragments] {msg}", file=sys.stderr)


def read_report(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(f"{path}: missing required column(s): {missing}")
        return list(reader)


def accession_from_href(href):
    """"/tabular-files/IGVFFI3117BGUY/@@download/IGVFFI3117BGUY.tsv.gz" ->
    ("IGVFFI3117BGUY", "IGVFFI3117BGUY.tsv.gz") -- the download filename is
    always the last path segment, and the accession is its stem before the
    first "." (every fragments file observed so far is "{accession}.tsv.gz")."""
    filename = href.rstrip("/").rsplit("/", 1)[-1]
    accession = filename.split(".", 1)[0]
    return accession, filename


def md5_of_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_one(conn, base_url, href, dest_path):
    """Streams href to dest_path via a temp file, so a failed/interrupted
    download never leaves a file that looks complete at dest_path (which
    would otherwise be mistaken for a verified download on the next run).

    Two requests, not one auto-following GET: the portal's @@download URl
    302s to a pre-signed S3 URL that carries its own signature in the query
    string. If our IGVF Basic Auth header rides along on that second
    request (which `requests` will do across the redirect in this setup),
    S3 sees two conflicting auth mechanisms and returns 403 -- confirmed
    against the real portal. So the redirect is followed by hand: the first
    (portal) request is authenticated and redirects disabled, the second
    (S3) request uses the Location header with no auth at all."""
    import requests
    import igvf_utils as iu
    from igvf_utils.utils import REQUEST_HEADERS_JSON

    url = urljoin(base_url, href)
    tmp_path = dest_path + ".partial"
    redirect_resp = requests.get(
        url, auth=conn.auth, headers=REQUEST_HEADERS_JSON, allow_redirects=False, timeout=iu.TIMEOUT
    )
    if redirect_resp.status_code in (301, 302, 303, 307, 308):
        signed_url = redirect_resp.headers["Location"]
    else:
        redirect_resp.raise_for_status()
        signed_url = url  # small/local files the portal serves directly, no redirect

    with requests.get(signed_url, stream=True, timeout=iu.TIMEOUT) as resp:
        resp.raise_for_status()
        with open(tmp_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    os.replace(tmp_path, dest_path)


def process_row(conn, base_url, outdir, row, force):
    href = row["FragmentsFile_Href"]
    expected_md5 = row["FragmentsFile_MD5"]
    alias = row.get("FragmentsFile_Aliases", "")
    if not href or not expected_md5:
        return {**row, "LocalPath": "", "ActualMD5": "", "Match": "SKIPPED_NO_HREF_OR_MD5"}

    accession, filename = accession_from_href(href)
    dest_path = os.path.join(outdir, filename)

    if not force and os.path.exists(dest_path):
        actual_md5 = md5_of_file(dest_path)
        if actual_md5 == expected_md5:
            log(f"{alias} ({accession}): already downloaded, md5 verified -- skipping")
            return {**row, "LocalPath": dest_path, "ActualMD5": actual_md5, "Match": "PASS_CACHED"}
        log(f"{alias} ({accession}): existing file md5 mismatch -- re-downloading")

    log(f"{alias} ({accession}): downloading...")
    try:
        download_one(conn, base_url, href, dest_path)
    except Exception as exc:  # noqa: BLE001 -- one bad file must not abort a 500-file batch
        log(f"{alias} ({accession}): DOWNLOAD FAILED: {exc}")
        return {**row, "LocalPath": dest_path, "ActualMD5": "", "Match": "DOWNLOAD_FAILED"}

    actual_md5 = md5_of_file(dest_path)
    match = "PASS" if actual_md5 == expected_md5 else "FAIL_MD5_MISMATCH"
    log(f"{alias} ({accession}): {match} (expected {expected_md5}, got {actual_md5})")
    return {**row, "LocalPath": dest_path, "ActualMD5": actual_md5, "Match": match}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("report", help="TSV with FragmentsFile_Href/FragmentsFile_MD5/FragmentsFile_Aliases columns")
    parser.add_argument("-o", "--outdir", required=True, help="destination directory (e.g. $SCRATCH/CATlas)")
    parser.add_argument("--igvf-mode", default="prod")
    parser.add_argument("--force", action="store_true", help="re-download even if a matching local file exists")
    parser.add_argument("--summary", default=None, help="optional path to write a per-file verification TSV")
    args = parser.parse_args()

    from igvf_utils.connection import Connection

    os.makedirs(args.outdir, exist_ok=True)
    rows = read_report(args.report)
    log(f"{len(rows)} row(s) in {args.report}")

    conn = Connection(igvf_mode=args.igvf_mode, no_log_file=True)
    base_url = conn.igvf_mode.url

    results = [process_row(conn, base_url, args.outdir, row, args.force) for row in rows]

    n_pass = sum(1 for r in results if r["Match"] in ("PASS", "PASS_CACHED"))
    n_fail = len(results) - n_pass
    log(f"{n_pass}/{len(results)} verified (md5 match); {n_fail} failed")

    if args.summary:
        cols = list(rows[0].keys()) + ["LocalPath", "ActualMD5", "Match"] if rows else []
        with open(args.summary, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for r in results:
                writer.writerow(r)
        log(f"wrote {args.summary}")

    if n_fail:
        for r in results:
            if r["Match"] not in ("PASS", "PASS_CACHED"):
                log(f"UNRESOLVED: {r.get('FragmentsFile_Aliases')}: {r['Match']}")
        sys.exit(1)


if __name__ == "__main__":
    main()
