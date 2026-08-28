#!/usr/bin/env python
"""Find portal-format reformatted files whose embedded Cell Annotation no longer
matches the CellAnnotation cache, and optionally delete them so Snakemake
rebuilds exactly those.

WHY THIS IS NEEDED AT ALL
-------------------------
rules/reformat.smk passes the Portal's cell metadata to its rules as `params:`
values, not as input files:

    summary=lambda wildcards: portal_cell_metadata(...)["cell_annotation"]

run_pipeline.py always passes `--rerun-triggers mtime`, deliberately (without it
an explicit --conda-prefix reads as "software environment changed" and queues a
full rebuild). mtime does NOT include params. So when the Portal corrects a Cell
Annotation and the cache is refreshed, every already-written reformatted file
still looks current: `snakemake portal_reformat` reports "Nothing to be done" and
the correction never reaches the files -- while the IGVF manifests DO pick it up,
because manage_igvf_metadata.py rebuilds those from scratch every run. The result
is manifests and data files that disagree about what cell type a prediction
describes. Observed for real on igvf2 and igvf3 (2026-08-28).

WHY DELETION RATHER THAN --forcerun OR A REAL input:
----------------------------------------------------
- `--forcerun` rebuilds every reformat output (828 files across 14 datasets), and
  its nargs='+' swallows a following positional target, which is the same footgun
  class as the --omit-from incident.
- Declaring the snapshot as a real `input:` looks like the Snakemake-native fix
  but is worse: cell_annotation_snapshot.write_snapshot is a full unconditional
  overwrite that stamps a fresh `derived_at` every run, so its mtime is always
  new and every reformat rule would rebuild on every invocation.
- Deleting only the mismatched outputs makes them genuinely missing, which
  `--rerun-triggers mtime` handles correctly and minimally.

Usable two ways: imported by run_pipeline.py (stage 1, straight after the snapshot
is written, so no extra Portal contact and well before the long compute), or as a
standalone CLI against state.db.
"""

import argparse
import gzip
import os
import sqlite3
import sys

# Header field -> the cell_annotations column it must agree with. All three are
# written by the reformat rules from the same portal_cell_metadata() call, but they
# do not always go stale together: a Portal term correction can change
# SampleTermID while leaving the annotation string identical, which a
# CellAnnotation-only check would pass silently.
HEADER_FIELDS = {
    "# CellAnnotation:": "cell_annotation",
    "# SampleTermName:": "term_name",
    "# SampleTermID:": "term_id",
}

# Reformatted outputs carry the header; scE2G's own outputs do not. Matched by
# suffix against the `{dataset}_{cluster}_` prefix the reformat rules use.
REFORMAT_SUFFIXES = (".e2g.tsv.gz", "_element_list.bed.gz", "_gene_list.tsv.gz")


def embedded_header(path):
    """{header field: value} for the HEADER_FIELDS this file carries. Reads only
    the leading comment block -- these files are large. Returns None if the file
    could not be read at all."""
    found = {}
    try:
        with gzip.open(path, "rt") as f:
            for line in f:
                if not line.startswith("#"):
                    break
                for key in HEADER_FIELDS:
                    if line.startswith(key):
                        found[key] = line[len(key):].strip()
    except OSError:
        return None
    return found


def _reformat_files(cluster_dir, dataset, cluster):
    if not os.path.isdir(cluster_dir):
        return []
    prefix = f"{dataset}_{cluster}_"
    return [
        os.path.join(cluster_dir, n)
        for n in sorted(os.listdir(cluster_dir))
        if n.startswith(prefix) and n.endswith(REFORMAT_SUFFIXES)
    ]


def find_stale(results_dir, dataset, annotation_rows):
    """annotation_rows: iterable of mappings with dataset/cluster plus the columns
    named in HEADER_FIELDS -- state.all_cell_annotations() rows, or the snapshot's
    own rows, interchangeably.

    Returns [(path, {header_field: (file_value, expected_value)})] for files whose
    header disagrees. Only fields a file actually carries are compared: the
    element BED omits some, and an older file may predate a header addition, and a
    field we cannot see cannot be judged stale.
    """
    stale = []
    for row in annotation_rows:
        if row.get("dataset") != dataset:
            continue
        cluster = row.get("cluster")
        if not row.get("cell_annotation"):
            continue  # never portal-annotated, therefore never reformatted
        want = {k: (row.get(col) or "") for k, col in HEADER_FIELDS.items()}
        for path in _reformat_files(os.path.join(results_dir, dataset, cluster), dataset, cluster):
            got = embedded_header(path)
            if not got:
                continue
            diffs = {k: (v, want[k]) for k, v in got.items() if v != want[k]}
            if diffs:
                stale.append((path, diffs))
    return stale


def remove_stale(stale, log=print):
    """Delete each stale file and its .tbi sibling if present. Returns the paths
    removed. The .tbi must go too: it indexes byte offsets into the file being
    rebuilt, so keeping it would leave a stale index beside fresh data."""
    removed = []
    for path, _ in stale:
        try:
            os.remove(path)
            removed.append(path)
        except OSError as e:
            log(f"could not remove {path}: {e}")
            continue
        tbi = path + ".tbi"
        if os.path.exists(tbi):
            try:
                os.remove(tbi)
                removed.append(tbi)
            except OSError as e:
                log(f"could not remove {tbi}: {e}")
    return removed


def describe(stale, limit=8):
    """Compact, log-friendly lines: one per stale file, capped."""
    out = []
    for path, diffs in stale[:limit]:
        fields = ", ".join(
            f"{k.strip('# :')}: {fv!r} -> {wv!r}" for k, (fv, wv) in sorted(diffs.items())
        )
        out.append(f"{os.path.basename(path)}  [{fields}]")
    if len(stale) > limit:
        out.append(f"... and {len(stale) - limit} more")
    return out


# --------------------------------------------------------------------------- CLI
def _rows_from_state_db(state_db, datasets):
    conn = sqlite3.connect(f"file:{state_db}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        cols = "dataset, cluster, cell_annotation, term_name, term_id"
        return [dict(r) for r in conn.execute(f"SELECT {cols} FROM cell_annotations")]
    finally:
        conn.close()


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("datasets", nargs="+")
    ap.add_argument("--state-db", default=os.path.join(repo, "resources", "igvf_metadata_state.db"))
    ap.add_argument("--results-dir", default=os.path.join(repo, "results", "uniformly_processed"))
    ap.add_argument("--delete", action="store_true",
                    help="remove the stale files so a portal_reformat pass rebuilds them")
    args = ap.parse_args()

    rows = _rows_from_state_db(args.state_db, args.datasets)
    total = []
    for dataset in args.datasets:
        stale = find_stale(args.results_dir, dataset, rows)
        checked = sum(len(_reformat_files(os.path.join(args.results_dir, dataset, r["cluster"]),
                                          dataset, r["cluster"]))
                      for r in rows if r["dataset"] == dataset and r.get("cell_annotation"))
        print(f"{dataset:<8} {checked - len(stale):>4} current, {len(stale):>3} stale")
        for line in describe(stale, limit=100):
            print(f"    {line}")
        total += stale

    print(f"\nTOTAL STALE: {len(total)} file(s)")
    if not total:
        print("nothing to do -- every reformatted file agrees with the cache")
        return 0
    if args.delete:
        removed = remove_stale(total)
        print(f"removed {len(removed)} file(s) (including .tbi siblings)")
        print("now run the reformat pass: --sce2g-modules false --snakemake-arg=portal_reformat")
    else:
        print("re-run with --delete to remove them, or let run_pipeline.py stage 1 do it")
    return 0


if __name__ == "__main__":
    sys.exit(main())
