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
    term_id TEXT NOT NULL,
    term_name TEXT NOT NULL,
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
    term_id TEXT,
    term_name TEXT,
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

-- Portal File objects belonging to PRIMARY pseudobulk sets, and the state of
-- downloading each one (2026-08-17). See portal_files.py for the discovery
-- query and download_portal_pseudobulks.py for the fetch loop.
--
-- Keyed on the portal's own file `accession`, NOT on (dataset, cluster) like
-- every other table here: `dataset` is an informal local label (igvf0..igvf18)
-- slated to be replaced by each dataset's principal analysis set accession
-- (see context.make_alias's UNRESOLVED note). A (dataset, ...) key would
-- orphan every row the day that rename lands; an accession never moves.
-- principal_analysis_set is captured now for the same reason -- it is the
-- future identity, available for free today from input_file_sets.
--
-- dataset/annotation/subsample are RESOLVED values, parsed from the file's
-- own submitted_file_name (the only field that carries all three), and are
-- nullable: a set whose submitted_file_name doesn't follow the convention is
-- recorded with them NULL and reported, never guessed at and never dropped.
--
-- href/file_size are nullable too -- href is a calculated portal field and a
-- not-yet-uploaded file may have neither.
CREATE TABLE IF NOT EXISTS portal_files (
    accession TEXT PRIMARY KEY,
    file_set TEXT,
    principal_analysis_set TEXT,
    lab TEXT,
    content_type TEXT NOT NULL,
    file_format TEXT,
    href TEXT,
    alias TEXT,
    md5sum TEXT,
    file_size INTEGER,
    portal_status TEXT,
    upload_status TEXT,
    submitted_file_name TEXT,
    dataset TEXT,
    annotation TEXT,
    subsample TEXT,
    local_path TEXT,
    -- Comma-joined tags from portal_files.resolve_scope: why a human should
    -- look at this row (irregular path shape, unexpected filename, incomplete
    -- parent set, ...). Non-empty does NOT mean unusable -- most such files
    -- still download fine; it means "surfaced, not silently normalised".
    review_reasons TEXT,
    -- pending | done | md5_mismatch | failed | needs_review | skipped
    download_state TEXT NOT NULL DEFAULT 'pending',
    bytes_written INTEGER,
    md5_observed TEXT,
    attempt_count INTEGER NOT NULL DEFAULT 0,
    last_error TEXT,
    first_seen_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_portal_files_state   ON portal_files(download_state);
CREATE INDEX IF NOT EXISTS idx_portal_files_scope   ON portal_files(dataset, annotation, subsample);
CREATE INDEX IF NOT EXISTS idx_portal_files_content ON portal_files(content_type);
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
    term_id,
    term_name,
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
               (dataset, cluster, cell_annotation, cl_id, term_id, term_name, cell_qualifier, portal_samples,
                all_primary_released, principal_uploaded, principal_alias, fetched_at)
               VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
               ON CONFLICT(dataset, cluster) DO UPDATE SET
                   cell_annotation=excluded.cell_annotation,
                   cl_id=excluded.cl_id,
                   term_id=excluded.term_id,
                   term_name=excluded.term_name,
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
                term_id,
                term_name,
                cell_qualifier,
                portal_samples,
                int(all_primary_released),
                int(principal_uploaded),
                principal_alias,
                now,
            ),
        )


def upsert_primary_pseudobulk(
    conn, alias, subsample, cell_annotation, cl_id, term_id, term_name, cell_qualifier, status, now
):
    with transaction(conn):
        conn.execute(
            """INSERT INTO cell_metadata_primary_pseudobulks
               (alias, subsample, cell_annotation, cl_id, term_id, term_name, cell_qualifier, status, fetched_at)
               VALUES (?,?,?,?,?,?,?,?,?)
               ON CONFLICT(alias) DO UPDATE SET
                   subsample=excluded.subsample,
                   cell_annotation=excluded.cell_annotation,
                   cl_id=excluded.cl_id,
                   term_id=excluded.term_id,
                   term_name=excluded.term_name,
                   cell_qualifier=excluded.cell_qualifier,
                   status=excluded.status,
                   fetched_at=excluded.fetched_at""",
            (alias, subsample, cell_annotation, cl_id, term_id, term_name, cell_qualifier, status, now),
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


_PORTAL_FILE_FIELDS = (
    "file_set",
    "principal_analysis_set",
    "lab",
    "content_type",
    "file_format",
    "href",
    "alias",
    "md5sum",
    "file_size",
    "portal_status",
    "upload_status",
    "submitted_file_name",
    "dataset",
    "annotation",
    "subsample",
    "local_path",
    "review_reasons",
)


def upsert_portal_file(conn, accession, now, **fields):
    """Records what the portal currently says about one File, WITHOUT touching
    any download-progress column (download_state/bytes_written/md5_observed/
    attempt_count/last_error). Discovery and downloading are separate passes:
    re-running discovery must refresh portal metadata -- including a changed
    md5sum, which is exactly how an updated upstream file gets noticed -- while
    leaving an in-flight or completed download's own bookkeeping intact.

    Deciding whether a refreshed md5sum invalidates an existing download is the
    downloader's job (see download_portal_pseudobulks.needs_download), not this
    function's: doing it here would mean a plain --dry-run discovery pass
    silently reset completed rows."""
    unknown = set(fields) - set(_PORTAL_FILE_FIELDS)
    if unknown:
        raise ValueError(f"unknown portal_files field(s): {sorted(unknown)}")
    cols = [c for c in _PORTAL_FILE_FIELDS if c in fields]
    assignments = ", ".join(f"{c}=excluded.{c}" for c in cols)
    placeholders = ", ".join("?" for _ in cols)
    with transaction(conn):
        conn.execute(
            f"""INSERT INTO portal_files (accession, {', '.join(cols)}, first_seen_at, updated_at)
                VALUES (?, {placeholders}, ?, ?)
                ON CONFLICT(accession) DO UPDATE SET {assignments}, updated_at=excluded.updated_at""",
            (accession, *[fields[c] for c in cols], now, now),
        )


def record_download_result(
    conn, accession, download_state, now, bytes_written=None, md5_observed=None, error=None, bump_attempt=True
):
    """The download-progress counterpart to upsert_portal_file. attempt_count is
    incremented by default; pass bump_attempt=False for a state change that
    wasn't an actual transfer attempt (e.g. marking a row 'skipped')."""
    with transaction(conn):
        conn.execute(
            f"""UPDATE portal_files
                SET download_state=?, bytes_written=COALESCE(?, bytes_written),
                    md5_observed=COALESCE(?, md5_observed), last_error=?,
                    attempt_count=attempt_count+{1 if bump_attempt else 0}, updated_at=?
                WHERE accession=?""",
            (download_state, bytes_written, md5_observed, error, now, accession),
        )


def get_portal_file(conn, accession):
    row = conn.execute("SELECT * FROM portal_files WHERE accession=?", (accession,)).fetchone()
    return dict(row) if row else None


def all_portal_files(conn, dataset=None, content_type=None, download_state=None):
    where, params = [], []
    for col, val in (("dataset", dataset), ("content_type", content_type), ("download_state", download_state)):
        if val is not None:
            where.append(f"{col}=?")
            params.append(val)
    sql = "SELECT * FROM portal_files"
    if where:
        sql += " WHERE " + " AND ".join(where)
    return [dict(r) for r in conn.execute(sql + " ORDER BY accession", params).fetchall()]


def portal_file_state_counts(conn):
    return {
        r["download_state"]: r["n"]
        for r in conn.execute(
            "SELECT download_state, COUNT(*) AS n FROM portal_files GROUP BY download_state"
        ).fetchall()
    }


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
