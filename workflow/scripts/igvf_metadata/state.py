"""SQLite ledger for the IGVF metadata uploader.

Row granularity: one state row per (dataset, cluster, model, table_name,
variant) -- a single failed/pending row never blocks or duplicates its
siblings. `model` and `variant` are normalized to "" rather than left NULL
when not applicable: SQLite's UNIQUE constraint treats NULL as distinct from
every other NULL, so two genuinely-duplicate rows with model=NULL would NOT
violate the UNIQUE(...) index -- "" avoids that trap.

The portal itself -- not this DB -- is the final authority on whether a row
really exists (see orchestrator.process_variant): this ledger is a fast-path
cache plus the audit trail, not a substitute for checking the destination.

cell_annotations (2026-07-21) is a different kind of table -- not an upload
ledger, but a 24h-TTL cache of the one PseudobulkSet-multireport GET per
pipeline trigger (see cell_metadata.py), one row per (dataset, cluster).
"""

import hashlib
import json
import os
import sqlite3
from contextlib import contextmanager

SCHEMA = """
CREATE TABLE IF NOT EXISTS clusters (
    dataset TEXT NOT NULL,
    cluster TEXT NOT NULL,
    excluded INTEGER NOT NULL DEFAULT 0,
    exclusion_reason TEXT,
    PRIMARY KEY (dataset, cluster)
);

CREATE TABLE IF NOT EXISTS uploads (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset TEXT NOT NULL,
    cluster TEXT NOT NULL,
    model TEXT NOT NULL DEFAULT '',
    table_name TEXT NOT NULL,
    variant TEXT NOT NULL DEFAULT '',
    alias TEXT NOT NULL UNIQUE,
    payload_hash TEXT NOT NULL,
    portal_id TEXT,
    status TEXT NOT NULL DEFAULT 'pending',
    attempt_count INTEGER NOT NULL DEFAULT 0,
    last_attempt_at TEXT,
    last_error TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE (dataset, cluster, model, table_name, variant)
);
CREATE INDEX IF NOT EXISTS idx_uploads_status  ON uploads(status);
CREATE INDEX IF NOT EXISTS idx_uploads_cluster ON uploads(dataset, cluster);
CREATE INDEX IF NOT EXISTS idx_uploads_table   ON uploads(table_name);

CREATE TABLE IF NOT EXISTS cell_annotations (
    dataset TEXT NOT NULL,
    cluster TEXT NOT NULL,
    cell_annotation TEXT NOT NULL,
    cl_id TEXT NOT NULL,
    cell_qualifier TEXT,
    portal_samples TEXT NOT NULL,
    all_primary_released INTEGER NOT NULL,
    principal_uploaded INTEGER NOT NULL,
    principal_alias TEXT,
    fetched_at TEXT NOT NULL,
    UNIQUE (dataset, cluster)
);

-- Tracks "did we issue the multireport GET recently" independent of whether
-- any (dataset, cluster) scope actually validated into a cell_annotations
-- row that round -- a run whose only scopes fail local-subset-of-portal
-- validation would otherwise leave cell_annotations empty forever, making
-- MAX(fetched_at) over that table read as "cache never populated" and
-- re-fetch on every single invocation. Single-row table (id always 1).
CREATE TABLE IF NOT EXISTS cell_annotations_fetch_log (
    id INTEGER PRIMARY KEY CHECK (id = 1),
    fetched_at TEXT NOT NULL
);

-- Every primary-pseudobulk row the multireport GET ever returned (2026-07-22),
-- saved unconditionally -- NOT gated on whether its subsample happens to match a
-- currently-configured local cluster, or whether that cluster's whole group
-- validates cleanly. cell_annotations (above) is a stricter, derived, per-cluster
-- view built FROM this table; this table is the actual "did we save what the
-- portal gave us" record, and what a downstream consumer (e.g. the E2G tabular
-- file reformatting step) should read for a given subsample's raw SampleSummaryShort/
-- CL Term ID/Cell Qualifier, independent of any local cluster-grouping concerns.
--
-- Keyed by alias, NOT subsample (2026-07-22 correction): a subsample (an
-- In-Vitro-System, MULTI-seq-tagged) is NOT a unique pseudobulk key -- the final
-- cell-type/cluster annotation is made by downstream human analysis, so one
-- subsample routinely has MANY distinct primary pseudobulks (one per resulting
-- cluster). The unique key is (subsample, cluster); today the alias
-- ("{lab}:{dataset}-{cluster}-{subsample}") is the only field that encodes which
-- cluster a given pseudobulk represents, so it doubles as that unique key. A
-- subsample-keyed table would silently keep only the last-seen pseudobulk per
-- subsample and discard the rest.
CREATE TABLE IF NOT EXISTS cell_metadata_primary_pseudobulks (
    alias TEXT PRIMARY KEY,
    subsample TEXT,
    cell_annotation TEXT,
    cl_id TEXT,
    cell_qualifier TEXT,
    status TEXT,
    fetched_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_cell_metadata_primary_subsample ON cell_metadata_primary_pseudobulks(subsample);

-- Every principal-pseudobulk row the multireport GET ever returned, saved
-- unconditionally -- the "is this Cell Annotation already locked in by an
-- uploaded principal" evidence, independent of local cluster-grouping.
CREATE TABLE IF NOT EXISTS cell_metadata_principal_pseudobulks (
    alias TEXT PRIMARY KEY,
    cell_annotation TEXT,
    fetched_at TEXT NOT NULL
);
"""


def connect(db_path: str) -> sqlite3.Connection:
    dirname = os.path.dirname(db_path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    conn = sqlite3.connect(db_path, timeout=30)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.executescript(SCHEMA)
    conn.commit()
    return conn


def payload_hash(payload: dict) -> str:
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


@contextmanager
def transaction(conn: sqlite3.Connection):
    conn.execute("BEGIN IMMEDIATE")
    try:
        yield conn
        conn.commit()
    except Exception:
        conn.rollback()
        raise


def get_upload(conn, dataset, cluster, model, table_name, variant):
    row = conn.execute(
        "SELECT * FROM uploads WHERE dataset=? AND cluster=? AND model=? AND table_name=? AND variant=?",
        (dataset, cluster, model or "", table_name, variant or ""),
    ).fetchone()
    return dict(row) if row else None


def claim_pending(conn, dataset, cluster, model, table_name, variant, alias, payload_hash_value, now):
    """Upsert this row to a fresh 'pending' claim ahead of any network call,
    inside its own short transaction (kept separate from the network I/O
    that follows, so we never hold the write lock during a slow HTTP call).
    A row already 'uploaded' with a matching hash is left untouched -- this
    only resets status when something's actually changed."""
    with transaction(conn):
        conn.execute(
            """INSERT INTO uploads
               (dataset, cluster, model, table_name, variant, alias, payload_hash,
                status, attempt_count, created_at, updated_at)
               VALUES (?,?,?,?,?,?,?, 'pending', 0, ?, ?)
               ON CONFLICT(dataset, cluster, model, table_name, variant) DO UPDATE SET
                   payload_hash=excluded.payload_hash,
                   updated_at=excluded.updated_at,
                   status = CASE WHEN status='uploaded' AND payload_hash=excluded.payload_hash
                                 THEN status ELSE 'pending' END""",
            (dataset, cluster, model or "", table_name, variant or "", alias, payload_hash_value, now, now),
        )
        row = conn.execute(
            "SELECT * FROM uploads WHERE dataset=? AND cluster=? AND model=? AND table_name=? AND variant=?",
            (dataset, cluster, model or "", table_name, variant or ""),
        ).fetchone()
        return dict(row)


def record_result(conn, row_id, status, portal_id=None, error=None, now=None):
    with transaction(conn):
        conn.execute(
            """UPDATE uploads SET status=?, portal_id=COALESCE(?, portal_id),
               last_error=?, attempt_count=attempt_count+1, last_attempt_at=?, updated_at=?
               WHERE id=?""",
            (status, portal_id, error, now, now, row_id),
        )


def latest_cell_annotation_fetch(conn):
    """When the multireport GET last ran -- cell_metadata.refresh_if_stale's
    24h TTL check. Reads cell_annotations_fetch_log, NOT MAX(fetched_at) over
    cell_annotations: a round where every scope fails validation caches zero
    rows, but the GET still happened and the TTL still applies."""
    row = conn.execute("SELECT fetched_at FROM cell_annotations_fetch_log WHERE id=1").fetchone()
    return row["fetched_at"] if row else None


def record_cell_annotations_fetch(conn, now):
    with transaction(conn):
        conn.execute(
            """INSERT INTO cell_annotations_fetch_log (id, fetched_at) VALUES (1, ?)
               ON CONFLICT(id) DO UPDATE SET fetched_at=excluded.fetched_at""",
            (now,),
        )


def get_cell_annotation(conn, dataset, cluster):
    row = conn.execute(
        "SELECT * FROM cell_annotations WHERE dataset=? AND cluster=?", (dataset, cluster)
    ).fetchone()
    return dict(row) if row else None


def all_cell_annotations(conn):
    return [dict(row) for row in conn.execute("SELECT * FROM cell_annotations").fetchall()]


def upsert_cell_annotation(
    conn,
    dataset,
    cluster,
    cell_annotation,
    cl_id,
    cell_qualifier,
    portal_samples,
    all_primary_released,
    principal_uploaded,
    principal_alias,
    now,
):
    with transaction(conn):
        conn.execute(
            """INSERT INTO cell_annotations
               (dataset, cluster, cell_annotation, cl_id, cell_qualifier, portal_samples,
                all_primary_released, principal_uploaded, principal_alias, fetched_at)
               VALUES (?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(dataset, cluster) DO UPDATE SET
                   cell_annotation=excluded.cell_annotation,
                   cl_id=excluded.cl_id,
                   cell_qualifier=excluded.cell_qualifier,
                   portal_samples=excluded.portal_samples,
                   all_primary_released=excluded.all_primary_released,
                   principal_uploaded=excluded.principal_uploaded,
                   principal_alias=excluded.principal_alias,
                   fetched_at=excluded.fetched_at""",
            (
                dataset,
                cluster,
                cell_annotation,
                cl_id,
                cell_qualifier,
                portal_samples,
                int(all_primary_released),
                int(principal_uploaded),
                principal_alias,
                now,
            ),
        )


def upsert_primary_pseudobulk(conn, alias, subsample, cell_annotation, cl_id, cell_qualifier, status, now):
    with transaction(conn):
        conn.execute(
            """INSERT INTO cell_metadata_primary_pseudobulks
               (alias, subsample, cell_annotation, cl_id, cell_qualifier, status, fetched_at)
               VALUES (?,?,?,?,?,?,?)
               ON CONFLICT(alias) DO UPDATE SET
                   subsample=excluded.subsample,
                   cell_annotation=excluded.cell_annotation,
                   cl_id=excluded.cl_id,
                   cell_qualifier=excluded.cell_qualifier,
                   status=excluded.status,
                   fetched_at=excluded.fetched_at""",
            (alias, subsample, cell_annotation, cl_id, cell_qualifier, status, now),
        )


def get_primary_pseudobulks_by_subsample(conn, subsample):
    """A subsample routinely has MANY distinct primary pseudobulks (one per
    resulting cluster) -- always a list, never a single row."""
    return [
        dict(row)
        for row in conn.execute(
            "SELECT * FROM cell_metadata_primary_pseudobulks WHERE subsample=?", (subsample,)
        ).fetchall()
    ]


def all_primary_pseudobulks(conn):
    return [dict(row) for row in conn.execute("SELECT * FROM cell_metadata_primary_pseudobulks").fetchall()]


def upsert_principal_pseudobulk(conn, alias, cell_annotation, now):
    with transaction(conn):
        conn.execute(
            """INSERT INTO cell_metadata_principal_pseudobulks (alias, cell_annotation, fetched_at)
               VALUES (?,?,?)
               ON CONFLICT(alias) DO UPDATE SET
                   cell_annotation=excluded.cell_annotation,
                   fetched_at=excluded.fetched_at""",
            (alias, cell_annotation, now),
        )


def all_principal_pseudobulks(conn):
    return [dict(row) for row in conn.execute("SELECT * FROM cell_metadata_principal_pseudobulks").fetchall()]


def mark_excluded(conn, dataset, cluster, reason, now):
    # Avoid referencing SQLite's "excluded" upsert pseudo-table in the DO UPDATE
    # clause -- our own column is also named "excluded", which is an
    # unnecessary ambiguity risk. Passing `reason` twice sidesteps it entirely.
    with transaction(conn):
        conn.execute(
            """INSERT INTO clusters (dataset, cluster, excluded, exclusion_reason)
               VALUES (?,?,1,?)
               ON CONFLICT(dataset, cluster) DO UPDATE SET excluded=1, exclusion_reason=?""",
            (dataset, cluster, reason, reason),
        )
