#!/usr/bin/env python3
"""Find CATlas fragment files that are duplicates of each other, keyed by portal alias.

    python workflow/scripts/find_duplicate_fragments.py \
        --state-db resources/igvf_metadata_state.db \
        --dataset catlas-human \
        --out resources/catlas-human_duplicate_fragments.tsv

WHY THIS EXISTS. The IGVF Portal's own duplicate check is `content_md5sum`,
computed over the decompressed bytes IN FILE ORDER. Two files holding the
identical fragments in a different row order therefore get different
content_md5sums and BOTH validate -- the Portal cannot see them. Confirmed pair
(both `validated`, neither flagged):

    yang-li:catlas-human-Chrp_77_fragments                     826713 bytes
    yang-li:catlas-human-supercluster-Choroid_plexus_fragments  838889 bytes
    -> different raw md5, different size, different content_md5sum
    -> IDENTICAL once sorted: 926dcc6b00703c94
    -> and the pipeline's own QC agrees independently: both report
       cell_count 17 / fragments_total 85543 in all_qc_stats.tsv

So "257 content-unique clusters" is only true up to row order. Predicting on both
members of such a pair yields two identical Prediction Sets.

DECISION (2026-09-04): duplicates are ACCEPTED, not excluded. The submitting lab
was asked and said both versions may be used, and which is preferred is still
undecided. This script therefore only RECORDS them -- it changes nothing and
excludes nothing. The record is what a later dedup decision would be made from.

WHAT IT HASHES, and why not the raw downloads. It hashes the SORTED,
decompressed content, preferring sources that already exist so the census costs
reads rather than a second 0.76 TiB of sorting:

  1. `{results}/{dataset}/{cluster}/fragments_filtered.tsv.gz` -- scE2G's own
     output, which is the sorted file with chromosomes restricted to chr_sizes.
     Persists after a run, so this is the cheap path.
  2. `{sorted_fragments_dir}/{dataset}/{cluster}_fragments.sorted.tsv.gz` -- the
     pipeline's sorted copy, if it still exists (it is temp(), so usually gone).
  3. the raw download, sorted on the fly -- correct but expensive, so it is
     opt-in via --allow-resort.

CAVEAT on source 1: it is chromosome-filtered as well as sorted. A pair that
differs ONLY on chromosomes absent from chr_sizes would look identical here. That
is a false positive for "the raw files are identical" but a TRUE positive for
"these two clusters feed the pipeline identical input", which is the thing that
produces duplicate predictions. The `source` column records which path was used
per cluster so a claim can always be traced.
"""

import argparse
import csv
import gzip
import hashlib
import os
import sqlite3
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pipeline_paths import repo_root_from_script, results_dir_base, sorted_fragments_dir  # noqa: E402

CHUNK = 1 << 20


def log(msg):
    print(f"[find_duplicate_fragments] {msg}", file=sys.stderr)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--state-db", required=True)
    p.add_argument("--dataset", default="catlas-human")
    p.add_argument("--out", required=True, help="TSV record of every cluster and its sorted-content hash")
    p.add_argument("--pipeline-config", default=None,
                   help="config to read output_dir/sorted_fragments_dir from; defaults to "
                        "workflow/config/{dataset}_pipeline_config.yaml")
    p.add_argument("--allow-resort", action="store_true",
                   help="if neither pre-sorted source exists for a cluster, sort the raw download "
                        "on the fly. Correct but expensive (0.76 TiB for catlas-human), so opt-in.")
    return p.parse_args(argv)


def _md5_stream(fh):
    h = hashlib.md5()
    for chunk in iter(lambda: fh.read(CHUNK), b""):
        h.update(chunk)
    return h.hexdigest()


def hash_sorted_content(cluster, results_base, dataset, sorted_dir, raw_path, allow_resort):
    """(md5_of_sorted_decompressed_content, source_label) or (None, reason)."""
    filtered = os.path.join(results_base, dataset, cluster, "fragments_filtered.tsv.gz")
    if os.path.exists(filtered):
        with gzip.open(filtered, "rb") as f:
            return _md5_stream(f), "fragments_filtered"

    if sorted_dir:
        pre = os.path.join(sorted_dir, dataset, f"{cluster}_fragments.sorted.tsv.gz")
        if os.path.exists(pre):
            with gzip.open(pre, "rb") as f:
                return _md5_stream(f), "sorted_copy"

    if not allow_resort:
        return None, "no_presorted_source"
    if not raw_path or not os.path.exists(raw_path):
        return None, "raw_missing"
    # Same sort the pipeline applies, so the hash is comparable with the other
    # two sources. Streamed through a pipe -- never materialised.
    cmd = f"zcat {raw_path!r} | LC_ALL=C sort -k1,1 -k2,2n"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    digest = _md5_stream(proc.stdout)
    if proc.wait() != 0:
        return None, "resort_failed"
    return digest, "resorted_raw"


def main(argv=None):
    args = parse_args(argv)
    repo = repo_root_from_script(os.path.abspath(__file__))
    cfg_path = args.pipeline_config or os.path.join(
        repo, "workflow", "config", f"{args.dataset}_pipeline_config.yaml"
    )
    import yaml

    with open(cfg_path) as f:
        config = yaml.safe_load(f) or {}
    results_base = results_dir_base(config, repo)
    sorted_dir = sorted_fragments_dir(config)

    conn = sqlite3.connect(args.state_db)
    rows = conn.execute(
        "SELECT annotation, alias, accession, local_path, file_size, review_reasons "
        "FROM portal_files WHERE dataset = ? AND content_type = 'fragments' ORDER BY annotation",
        (args.dataset,),
    ).fetchall()
    conn.close()
    if not rows:
        raise SystemExit(f"no fragments rows for dataset {args.dataset!r} in {args.state_db}")

    records, skipped = [], []
    for cluster, alias, accession, raw_path, size, reasons in rows:
        digest, source = hash_sorted_content(
            cluster, results_base, args.dataset, sorted_dir, raw_path, args.allow_resort
        )
        if digest is None:
            skipped.append((alias, source))
            continue
        records.append({
            "portal_alias": alias,
            "cluster": cluster,
            "file_accession": accession,
            "sorted_content_md5": digest,
            "file_size": size,
            "portal_flagged": "yes" if "upload_status:" in (reasons or "") else "no",
            "hash_source": source,
        })
        log(f"{alias}: {digest} ({source})")

    by_hash = {}
    for r in records:
        by_hash.setdefault(r["sorted_content_md5"], []).append(r)
    dupes = {h: v for h, v in by_hash.items() if len(v) > 1}

    for r in records:
        group = by_hash[r["sorted_content_md5"]]
        r["duplicate_group_size"] = len(group)
        r["duplicates_of"] = ",".join(
            o["portal_alias"] for o in group if o["portal_alias"] != r["portal_alias"]
        )

    header = ["portal_alias", "cluster", "file_accession", "sorted_content_md5", "file_size",
              "duplicate_group_size", "duplicates_of", "portal_flagged", "hash_source"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in sorted(records, key=lambda x: (-x["duplicate_group_size"], x["portal_alias"])):
            w.writerow(r)

    log("")
    log(f"hashed {len(records)} of {len(rows)} cluster(s); {len(skipped)} not yet hashable")
    if skipped:
        by_reason = {}
        for alias, reason in skipped:
            by_reason.setdefault(reason, []).append(alias)
        for reason, aliases in by_reason.items():
            log(f"  {reason}: {len(aliases)} (e.g. {aliases[0]})")
        log("  re-run once the pipeline has produced their fragments_filtered.tsv.gz,")
        log("  or pass --allow-resort to sort the raw downloads on the fly.")
    log("")
    log(f"DUPLICATE GROUPS FOUND: {len(dupes)}")
    for h, group in sorted(dupes.items()):
        flagged = [g for g in group if g["portal_flagged"] == "yes"]
        log(f"  {h}  ({len(group)} files, {len(flagged)} portal-flagged)")
        for g in sorted(group, key=lambda x: x["portal_alias"]):
            log(f"     {g['portal_alias']}  [{g['file_accession']}, {g['file_size']} bytes,"
                f" portal_flagged={g['portal_flagged']}]")
    log("")
    log(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
