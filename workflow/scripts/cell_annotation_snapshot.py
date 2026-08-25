"""The per-dataset CellAnnotation snapshot: a read-only snapshot of state.db's
`cell_annotations` rows that Snakemake reads INSTEAD of opening state.db.

Why this exists at all -- it is not a convenience, it's a concurrency fix.
state.db lives on Lustre in WAL mode. WAL's index is an mmap'd -shm file that
every connection must share, which is only supported when all connections are on
ONE host. But common.smk opened state.db at PARSE time, and every Slurm worker
re-parses the whole Snakefile (see common.smk's KNOWN LIMITATION comment) -- so
with several dataset drivers running at once, driver N's writes overlap driver
N-1's worker-node reads, across hosts. Readers never conflict with readers; reads
overlapping a write across hosts is the unsupported case.

With this snapshot, no worker node ever opens state.db. The only processes that
touch it are the drivers, which serialise on an flock, so state.db has at most one
accessor at any instant -- sidestepping WAL/shm/Lustre coherency entirely rather
than reasoning about which overlaps happen to be safe.

TEMP ARTIFACT, deliberately. An earlier version of this pipeline gated portal
reformatting on a static, manually-refreshed preview TSV and crashed for real when
that TSV claimed a cluster was annotated while state.db was actually cold (commit
b6e62e0). The rules that make this different from that TSV, and they are
load-bearing:

  - Written unconditionally fresh by the driver immediately before Snakemake
    starts. Never read-if-exists, never merged with a previous run's content.
  - Carries the portal fetch timestamp and a digest of the cluster set it was
    built for, and read_snapshot() ENFORCES both. A leftover or foreign file
    fails loudly instead of silently answering the wrong question.
  - Deleted by the driver when the Snakemake stage exits. A later bare
    `snakemake` in default mode therefore finds nothing and aborts pointing at
    the driver, rather than silently dropping reformat targets.
  - Declared as neither input nor output of any rule, so Snakemake's
    --rerun-triggers machinery never sees it.

state.db's cell_annotations table remains the authority: manifest generation reads
it directly (prediction_set and principal_pseudobulk_set, via
cell_metadata.get_metadata_for). This is a snapshot of that, never a substitute.
"""

import csv
import hashlib
import os
from datetime import datetime, timezone

SNAPSHOT_DIR_NAME = "igvf_metadata"
DEFAULT_MAX_AGE_HOURS = 24

# Mirrors state.db's cell_annotations columns. term_id/term_name/cell_annotation are
# what reformat.smk actually puts in the file headers; the rest travel along so the
# snapshot is a faithful snapshot rather than a lossy subset.
COLUMNS = [
    "dataset",
    "cluster",
    "cell_annotation",
    "cl_id",
    "term_id",
    "term_name",
    "cell_qualifier",
    "portal_samples",
    "all_primary_released",
    "principal_uploaded",
    "principal_alias",
]

_PROVENANCE_PREFIX = "#"


def snapshot_path(output_dir, dataset):
    return os.path.join(output_dir, SNAPSHOT_DIR_NAME, f"{dataset}_cell_annotations.tsv")


def cluster_set_digest(cluster_keys):
    """Stable digest of a (dataset, cluster) set. Both the driver (writing) and
    Snakemake's parse (reading) compute this from the same pipeline config, so a
    snapshot built for a different cluster set can never be silently accepted --
    e.g. after someone edits the config between the driver run and a later bare
    `snakemake` invocation."""
    joined = ";".join(sorted(f"{dataset}/{cluster}" for dataset, cluster in cluster_keys))
    return hashlib.sha256(joined.encode("utf-8")).hexdigest()[:16]


class SnapshotError(RuntimeError):
    """Raised with an actionable message -- always naming the command that fixes
    it, since the whole point is that a cold/stale/foreign cache must never
    degrade quietly into "no reformat targets"."""


def write_snapshot(path, rows, fetched_at, digest):
    """Full, unconditional overwrite. rows: dicts with at least COLUMNS' keys
    (extra keys ignored, missing ones written blank)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    derived_at = datetime.now(timezone.utc).isoformat()
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w", newline="") as f:
        f.write(f"{_PROVENANCE_PREFIX}portal_fetched_at\t{fetched_at}\n")
        f.write(f"{_PROVENANCE_PREFIX}derived_at\t{derived_at}\n")
        f.write(f"{_PROVENANCE_PREFIX}cluster_set_digest\t{digest}\n")
        writer = csv.DictWriter(f, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["dataset"], r["cluster"])):
            writer.writerow({c: ("" if row.get(c) is None else row.get(c)) for c in COLUMNS})
    os.replace(tmp, path)  # atomic: worker nodes may be reading a previous run's file
    return path


def _parse(path):
    provenance = {}
    data_lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(_PROVENANCE_PREFIX):
                key, _, value = line[len(_PROVENANCE_PREFIX):].rstrip("\n").partition("\t")
                provenance[key] = value
                continue
            data_lines.append(line)
    rows = list(csv.DictReader(data_lines, delimiter="\t"))
    return provenance, rows


def read_snapshot(path, expected_digest=None, max_age_hours=DEFAULT_MAX_AGE_HOURS, fix_hint=""):
    """{(dataset, cluster): row} for a snapshot that passes every check.

    Raises SnapshotError -- never returns an empty dict as a stand-in for
    "couldn't read it". Silently returning nothing is exactly the failure mode
    this whole mechanism exists to remove: it looks identical to "no cluster is
    annotated yet" and makes the reformat targets vanish without a word.
    """
    hint = f"\n  Fix: {fix_hint}" if fix_hint else ""
    if not os.path.exists(path):
        raise SnapshotError(f"no CellAnnotation snapshot at {path}{hint}")

    provenance, rows = _parse(path)

    if expected_digest is not None and provenance.get("cluster_set_digest") != expected_digest:
        raise SnapshotError(
            f"{path} was built for a different cluster set "
            f"(snapshot {provenance.get('cluster_set_digest')!r} != config {expected_digest!r}) -- "
            f"the pipeline config changed since it was written{hint}"
        )

    fetched_at = provenance.get("portal_fetched_at")
    if max_age_hours is not None:
        if not fetched_at:
            raise SnapshotError(f"{path} records no portal_fetched_at, so its age can't be checked{hint}")
        age_hours = (datetime.now(timezone.utc) - datetime.fromisoformat(fetched_at)).total_seconds() / 3600
        if age_hours > max_age_hours:
            raise SnapshotError(
                f"{path} is built on a portal fetch from {fetched_at} "
                f"({age_hours:.1f}h old, limit {max_age_hours}h){hint}"
            )

    return {(row["dataset"], row["cluster"]): row for row in rows}


def remove_snapshot(path):
    """Best-effort delete -- the driver calls this in a finally block, where a
    missing file is a normal outcome (the run may have failed before writing it)."""
    try:
        os.remove(path)
        return True
    except FileNotFoundError:
        return False
