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
