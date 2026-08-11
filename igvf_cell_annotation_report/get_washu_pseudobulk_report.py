#!/usr/bin/env python3
"""One-off IGVF Data Portal query for the Yang Li (WashU) PseudobulkSets.

Mirrors the GET pattern in get_igvf_cell_annotation_report.py (same auth via
igvf_utils.Connection, same multireport endpoint, limit=all) but scoped to
lab.title=Yang Li, WashU and reporting a different set of fields: full
alias, cell_annotation, cell_qualifier, cell_type (term_id/term_name),
input_file_sets (to confirm they're curated sets, not analysis sets), and
the files list (to find the "fragments" file's aliases/status/href/md5sum).
"""

import json
import sys
from urllib.parse import quote

QUERY_FIELDS = [
    "@id",
    "aliases",
    "cell_annotation",
    "cell_qualifier",
    "cell_type",
    "cell_type.term_name",
    "cell_type.term_id",
    "status",
    "input_file_sets",
    "files",
    "files.content_type",
    "files.aliases",
    "files.status",
    "files.href",
    "files.md5sum",
]


def fetch_multireport(igvf_mode="prod"):
    import requests
    import igvf_utils as iu
    import igvf_utils.utils as iuu
    from igvf_utils.connection import Connection

    conn = Connection(igvf_mode=igvf_mode, no_log_file=True)
    lab_title = quote("Yang Li, WashU")
    query = (
        f"type=PseudobulkSet&status%21=deleted&lab.title={lab_title}&limit=all"
        + "".join(f"&field={quote(f)}" for f in QUERY_FIELDS)
    )
    url = iuu.url_join([conn.igvf_mode.url, "multireport/?"]) + query
    response = requests.get(url, auth=conn.auth, timeout=iu.TIMEOUT, headers=iuu.REQUEST_HEADERS_JSON)
    response.raise_for_status()
    return response.json()["@graph"]


def main():
    rows = fetch_multireport()
    print(f"# {len(rows)} PseudobulkSet(s) returned", file=sys.stderr)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
