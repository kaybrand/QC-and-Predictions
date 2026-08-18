"""Streaming, md5-verifying, resumable file download from the IGVF Portal.

Deliberately does NOT use igvf_utils' Connection.download(), which is unusable
at this scale (~3100 files, some multi-GB). Confirmed by reading
igvf_utils/connection.py at v3.1.1:
  - chunk_size=512 bytes (connection.py:1943) -- 2M write syscalls per GB;
  - no md5 verification of the downloaded bytes anywhere;
  - filename parsed from Content-Disposition via a bare
    r.headers["Content-Disposition"] (connection.py:1936), a KeyError if the
    header is absent, and unquoted;
  - verify=False (connection.py:1930), silently disabling TLS verification;
  - get_stream=True returns AFTER open(filename,"wb") (connection.py:1940-41),
    leaving a 0-byte file behind;
  - and Connection has no retry/backoff of any kind -- no Session, no
    HTTPAdapter, no urllib3 Retry, no 429/503 handling.

This module keeps TLS verification ON, matching portal_client.get_multireport's
own deliberate choice not to copy that library-internal skip.

Crash safety: bytes are streamed to "{dest}.part" and only os.replace()d into
place after the md5 matches. os.replace is atomic within a filesystem, so an
interrupted or corrupt transfer can never leave a truncated file that a later
run would mistake for complete. A failed md5 leaves the .part in place for
inspection and does NOT publish.

Two layers of retry, because they cover different failures:
  1. urllib3 Retry on the adapter -- retries the REQUEST when the response
     itself fails to establish (429/5xx, connection refused).
  2. _stream_to_part's own loop with byte-range resume -- retries the BODY.
     This layer is the important one here and is easy to miss: `href` answers
     with a 307 cross-host redirect to a presigned S3 URL, and requests passes
     redirect=False down to urllib3, so adapter-level Retry only ever covers
     one hop's response establishment. A connection dropped partway through
     iter_content() of a multi-GB body raises ChunkedEncodingError/
     ProtocolError and is NOT retried by layer 1 at all. For thousands of large
     files that is the single most likely failure mode.

Auth across that redirect is safe: requests' Session.rebuild_auth() strips the
Authorization header when a redirect changes host, so Basic credentials are
never sent to S3. Hence allow_redirects stays True and auth lives on the
Session, not on a manually-built header.
"""

import hashlib
import os
import re
import socket
import sys
import time
from http.client import IncompleteRead

CHUNK_SIZE = 1 << 20  # 1 MiB
# (connect, read). The read timeout is per-socket-read INACTIVITY, not a cap on
# total transfer duration, so a large but healthy body is never cut short.
DEFAULT_TIMEOUT = (10, 300)
MAX_STREAM_ATTEMPTS = 5


def log(msg):
    print(f"[downloader] {msg}", file=sys.stderr)


def redact(text):
    """Strip query strings out of any URL in a message before it is logged or
    stored.

    `href` redirects to a PRESIGNED S3 URL whose query string carries
    AWSAccessKeyId, Signature and x-amz-security-token. requests puts the full
    final URL into its exception text, so an unredacted error ends up in the
    state.db ledger and the job log -- which is how a temporary AWS credential
    got persisted to disk on the first full run. The tokens are short-lived and
    object-scoped, and are not the IGVF key pair, but they have no business being
    written down. The path is kept, since that is the part with diagnostic value.
    """
    if not text:
        return text
    return re.sub(r"(https?://[^\s?]+)\?[^\s]*", r"\1?<redacted>", str(text))


class Result:
    """state: "done" | "md5_mismatch" | "failed". `verified` is True only when an
    expected md5 was supplied AND matched -- a portal file with no md5sum can be
    fetched but not proven correct, and the caller decides what to do about that
    rather than having it silently count as success."""

    __slots__ = ("state", "bytes_written", "md5_observed", "error", "verified")

    def __init__(self, state, bytes_written=None, md5_observed=None, error=None, verified=False):
        self.state = state
        self.bytes_written = bytes_written
        self.md5_observed = md5_observed
        self.error = error
        self.verified = verified

    def __repr__(self):
        return (
            f"Result(state={self.state!r}, bytes={self.bytes_written}, "
            f"md5={self.md5_observed!r}, verified={self.verified}, error={self.error!r})"
        )


def build_session(auth, total_retries=5, backoff_factor=1.0, pool_maxsize=8):
    """A Session with adapter-level retry for response establishment. See the
    module docstring for why this is necessary but not sufficient."""
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    retry = Retry(
        total=total_retries,
        connect=total_retries,
        read=total_retries,
        status=total_retries,
        backoff_factor=backoff_factor,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET", "HEAD"}),
        respect_retry_after_header=True,
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=pool_maxsize, pool_maxsize=pool_maxsize)
    session = requests.Session()
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    session.auth = auth
    return session


def resolve_url(base_url, href):
    """href is a portal PATH, e.g.
    "/matrix-files/IGVFFI0694XKXE/@@download/IGVFFI0694XKXE.h5ad".
    Joined onto the same api.data.igvf.org base the metadata came from."""
    if not href:
        raise ValueError("file has no href")
    if href.startswith("http://") or href.startswith("https://"):
        return href
    return base_url.rstrip("/") + "/" + href.lstrip("/")


def _hash_existing(path, hasher, chunk_size=CHUNK_SIZE):
    """Feed an existing partial file into `hasher` and return its length. Called
    once per resume attempt: the md5 must be rebuilt from the bytes actually on
    disk, since a previous attempt's in-memory hasher is long gone."""
    size = 0
    with open(path, "rb") as fh:
        while True:
            block = fh.read(chunk_size)
            if not block:
                break
            hasher.update(block)
            size += len(block)
    return size


def _retryable(exc):
    import requests

    return isinstance(
        exc,
        (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ConnectionError,
            requests.exceptions.Timeout,
            IncompleteRead,
            socket.timeout,
        ),
    )


def _stream_to_part(session, url, part_path, expected_size, chunk_size, timeout, max_attempts, sleep):
    """Stream `url` into `part_path`, resuming by byte range across attempts.
    Returns (bytes_on_disk, md5_hexdigest). Raises on unrecoverable failure."""
    import requests

    last_exc = None
    for attempt in range(1, max_attempts + 1):
        hasher = hashlib.md5()
        offset = _hash_existing(part_path, hasher) if os.path.exists(part_path) else 0
        if expected_size and offset == expected_size:
            return offset, hasher.hexdigest()  # already fully on disk
        if expected_size and offset > expected_size:
            # Overlong partial: can't be a prefix of the real object. Start over
            # rather than trying to reason about it.
            log(f"  partial is longer than expected ({offset} > {expected_size}); restarting from 0")
            os.remove(part_path)
            offset, hasher = 0, hashlib.md5()

        headers = {"Range": f"bytes={offset}-"} if offset else {}
        try:
            with session.get(
                url, stream=True, timeout=timeout, headers=headers, allow_redirects=True
            ) as resp:
                resp.raise_for_status()
                mode = "ab"
                if offset and resp.status_code != 206:
                    # Range ignored (200 with the whole body). Anything already
                    # written would be duplicated, so discard and restart.
                    log(f"  server ignored Range (HTTP {resp.status_code}); restarting from 0")
                    offset, hasher, mode = 0, hashlib.md5(), "wb"
                written = offset
                with open(part_path, mode) as out:
                    for chunk in resp.iter_content(chunk_size=chunk_size):
                        if not chunk:
                            continue
                        out.write(chunk)
                        hasher.update(chunk)
                        written += len(chunk)
                    out.flush()
                    os.fsync(out.fileno())
            if expected_size and written != expected_size:
                raise requests.exceptions.ChunkedEncodingError(
                    f"short body: got {written} of {expected_size} bytes"
                )
            return written, hasher.hexdigest()
        except Exception as exc:  # noqa: BLE001 -- re-raised below unless retryable
            last_exc = exc
            if not _retryable(exc) or attempt == max_attempts:
                raise
            delay = min(2 ** (attempt - 1), 30)
            have = os.path.getsize(part_path) if os.path.exists(part_path) else 0
            log(
                f"  attempt {attempt}/{max_attempts} failed ({type(exc).__name__}: {redact(exc)}); "
                f"{have} bytes on disk, resuming in {delay}s"
            )
            sleep(delay)
    raise last_exc  # unreachable; kept so a future edit can't fall through silently


def download(
    session,
    base_url,
    href,
    dest_path,
    expected_md5=None,
    expected_size=None,
    chunk_size=CHUNK_SIZE,
    timeout=DEFAULT_TIMEOUT,
    max_attempts=MAX_STREAM_ATTEMPTS,
    sleep=time.sleep,
):
    """Fetch one file to dest_path, atomically and md5-verified.

    Never raises for an expected failure -- returns a Result whose state the
    caller records. A wrong md5 leaves "{dest_path}.part" on disk and does not
    publish, because a false "this file is fine" is the costliest outcome here.
    """
    part_path = dest_path + ".part"
    try:
        url = resolve_url(base_url, href)
    except ValueError as exc:
        return Result("failed", error=str(exc))

    parent = os.path.dirname(dest_path)
    if parent:
        os.makedirs(parent, exist_ok=True)

    try:
        written, md5_observed = _stream_to_part(
            session, url, part_path, expected_size, chunk_size, timeout, max_attempts, sleep
        )
    except Exception as exc:  # noqa: BLE001 -- surfaced as a recorded failure
        return Result("failed", error=redact(f"{type(exc).__name__}: {exc}"))

    if expected_md5 and md5_observed != expected_md5:
        return Result(
            "md5_mismatch",
            bytes_written=written,
            md5_observed=md5_observed,
            error=(
                f"md5 mismatch: portal={expected_md5} observed={md5_observed} "
                f"({written} bytes); left {os.path.basename(part_path)} in place"
            ),
        )

    os.replace(part_path, dest_path)
    return Result("done", bytes_written=written, md5_observed=md5_observed, verified=bool(expected_md5))


def needs_download(row, dest_path):
    """Should this file be fetched? Returns (bool, reason).

    The md5 clause is the point of this function: a row already 'done' is only
    skipped if the md5 we OBSERVED when we downloaded it still equals the md5
    the portal reports NOW. When upstream replaces a file, discovery refreshes
    portal_files.md5sum, the two stop agreeing, and the file is re-fetched. A
    plain "state == done" check would skip it forever -- and a false "unchanged"
    is worse than a false "changed", since it means shipping stale data as
    current.
    """
    if not os.path.exists(dest_path):
        return True, "absent"
    if row.get("download_state") != "done":
        return True, f"state={row.get('download_state')}"
    expected_size = row.get("file_size")
    if expected_size is not None and os.path.getsize(dest_path) != expected_size:
        return True, f"size {os.path.getsize(dest_path)} != portal {expected_size}"
    portal_md5, observed = row.get("md5sum"), row.get("md5_observed")
    if portal_md5 and observed and portal_md5 != observed:
        return True, "portal md5 changed since download"
    if portal_md5 and not observed:
        return True, "no recorded md5 to compare"
    return False, "unchanged"
