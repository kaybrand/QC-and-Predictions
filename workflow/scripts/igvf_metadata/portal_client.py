"""Portal access for the IGVF metadata uploader: a read-only existence
check via igvf_utils.connection.Connection, plus TSV generation and a
subprocess wrapper around the real registration script --

    /oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py

Confirmed by reading that script and igvf_utils/connection.py directly
(2026-07-13):
  - PROFILE_KEY = "_profile", IGVFID_KEY = "_igvf_id" -- NOT "@type"/"uuid"
    as an earlier draft of this file guessed.
  - Connection.get(rec_ids, ignore404=True) returns {} (falsy) when not
    found, not an exception -- matches get_by_alias's `or None` below.
  - Real dry-run behavior lives inside Connection itself (constructed with
    dry_run=...); iu_register.py never fakes it at the script level, so we
    don't either.
  - iu_register.py POSTs/PATCHes an entire TSV/JSON/JSONL file against ONE
    profile at a time; PATCH rows are identified by a literal `record_id`
    column in the input file (RECORD_ID_FIELD below), which iu_register.py
    itself translates to IGVFID_KEY internally -- our TSV writer must use
    "record_id", not "_igvf_id".
  - iu_register.py's own TSV reader (create_payloads_from_tsv) parses each
    line with a plain `line.split("\t")` -- no csv module, no quote/escape
    handling at all -- then calls json.loads directly on the raw field text
    for object/array-of-object fields. write_tsv below must therefore never
    csv-quote a field: Python's csv.writer's default QUOTE_MINIMAL wraps any
    field containing a `"` (every JSON-dumped dict/list value) in an outer
    quote pair and doubles internal quotes, which iu_register.py never
    unescapes -- json.loads then fails immediately on the literal leading
    `"` (confirmed 2026-08-05: a real --mode upload attempt failed with
    "Extra data: line 1 column 4", i.e. json.loads successfully parsed the
    3-character string `"{"` from `"{""path""...` and choked on what
    followed). Plain str.join, no quoting, is correct here precisely because
    none of our field values ever contain a literal tab.
  - It POSTs/PATCHes rows in a single Python loop with no per-row
    try/except beyond a JSONDecodeError check -- an error partway through
    one file can abort the remaining rows in that file. Callers of
    invoke_register in "upload" mode should re-verify each intended row via
    get_by_alias afterward rather than trusting the subprocess's exit code.

get_multireport() is the one other read path here (2026-07-21): a raw GET
against /multireport/ for cell_metadata.py's PseudobulkSet lookup --
Connection.search() can't be reused for this, it hardcodes "search/" as
the path.
"""

import json
import os
import subprocess
import sys

IU_REGISTER_DEFAULT_PATH = (
    "/oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/igvf_utils/igvf_utils/MetaDataRegistration/iu_register.py"
)

RECORD_ID_FIELD = "record_id"  # must match iu_register.py's RECORD_ID_FIELD exactly


class PortalReader:
    """Read-only: the idempotency check of record. Always trust the
    portal's own answer over the local state ledger, since the ledger can
    lag behind a crash between "portal accepted the write" and "ledger
    commit". Cheap because aliases are directly addressable -- no
    full-tree listing needed the way Synapse's getChildren requires.
    Never performs a write; iu_register.py (via invoke_register) is the
    only thing in this package that does."""

    def __init__(self, igvf_mode=None):
        self.igvf_mode = igvf_mode
        self._conn = None

    def _connection(self):
        if self._conn is None:
            from igvf_utils.connection import Connection  # deferred: only needed once we actually query

            self._conn = Connection(igvf_mode=self.igvf_mode)
        return self._conn

    def get_by_alias(self, alias: str):
        """Returns the portal record dict, or None if it doesn't exist yet."""
        return self._connection().get(alias, ignore404=True) or None

    def get_multireport(self, query_string: str):
        """One raw GET against the /multireport/ endpoint -- e.g. cell_metadata.py's
        PseudobulkSet lookup. Not reusable via Connection.search(): that method
        hardcodes its URL to "search/?..." via make_search_url() and can't be
        pointed at a different report path. Mirrors search()'s own request shape
        (auth/timeout/headers), but does NOT copy its verify=False -- that
        library-internal TLS-verification skip is the vendored code's own
        choice, not something to reintroduce here without asking. If a real
        run against data.igvf.org hits a cert-verification error, surface it
        rather than silently disabling verification.
        Returns the parsed "@graph" list, same shape search() returns."""
        import requests
        import igvf_utils as iu
        import igvf_utils.utils as iuu

        conn = self._connection()
        url = iuu.url_join([conn.igvf_mode.url, "multireport/?"]) + query_string
        response = requests.get(url, auth=conn.auth, timeout=iu.TIMEOUT, headers=iuu.REQUEST_HEADERS_JSON)
        response.raise_for_status()
        return response.json()["@graph"]


def _tsv_cell(value):
    """Follows iu_register.py's create_payloads_from_tsv conventions: array
    of scalars -> comma-joined (brackets optional, so we just omit them);
    array of objects / bare object -> JSON; everything else -> str()."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple)):
        if value and isinstance(value[0], dict):
            return json.dumps(list(value))
        return ",".join(str(v) for v in value)
    if isinstance(value, dict):
        return json.dumps(value)
    if isinstance(value, bool):
        return str(value).lower()
    return str(value)


def write_tsv(path, rows, record_ids=None):
    """rows: list[dict] of schema-property -> value, NOT including the
    _profile key -- iu_register.py takes that as its --profile_id CLI arg,
    one profile per file, and injects it into every payload itself.

    record_ids: optional list, parallel to rows, of portal identifiers to
    PATCH. When given, a `record_id` column is added -- iu_register.py's
    own PATCH-row convention, not a name we're free to change.

    Always written (never gated on upload permission): this is exactly the
    file "for the user to peruse" in preview mode, and the same file gets
    handed to iu_register.py verbatim when upload is actually permitted.
    """
    if not rows:
        return None
    columns = sorted({k for row in rows for k in row})
    if record_ids is not None:
        columns = [RECORD_ID_FIELD] + columns
    dirname = os.path.dirname(path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    # Plain str.join, NOT csv.writer -- see module docstring for why
    # csv-quoting corrupts iu_register.py's own naive line.split("\t") parse.
    with open(path, "w", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for i, row in enumerate(rows):
            values = dict(row)
            if record_ids is not None:
                values[RECORD_ID_FIELD] = record_ids[i]
            f.write("\t".join(_tsv_cell(values.get(c)) for c in columns) + "\n")
    return path


def invoke_register(
    infile, profile_id, patch=False, dry_run=True, igvf_mode=None, iu_register_path=IU_REGISTER_DEFAULT_PATH
):
    """Shells out to the real igvf_utils registration script -- this is the
    ONLY code path in this package that can ever perform a live write to
    the portal, and only when dry_run=False. dry_run=True still invokes
    the script (exercising its real schema validation/type-casting) but
    passes --dry-run, so igvf_utils itself guarantees no write happens.

    Returns the completed subprocess.CompletedProcess. See this module's
    docstring: iu_register.py doesn't isolate one bad row from the rest of
    its input file, so callers running for real should re-verify each row
    via PortalReader.get_by_alias afterward rather than trusting returncode
    alone.
    """
    cmd = [sys.executable, iu_register_path, "--profile_id", profile_id, "--infile", infile]
    if patch:
        cmd.append("--patch")
    if dry_run:
        cmd.append("--dry-run")
    if igvf_mode:
        cmd += ["--igvf-mode", igvf_mode]
    return subprocess.run(cmd, capture_output=True, text=True)
