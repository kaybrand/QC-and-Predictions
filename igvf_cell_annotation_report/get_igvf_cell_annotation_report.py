#!/usr/bin/env python3
"""Standalone IGVF Data Portal cell-annotation report.

Produces one TSV row per {dataset}-{cluster}, surveying every primary
pseudobulk PseudobulkSet currently on the portal: CellAnnotation,
SampleTermID, SampleTermName, CellQualifier, the list of contributing
subsamples, a released/in-progress status, and (see MODES below) which QC
guide file, if any, was used to resolve that row.

Deliberately standalone -- no dependency on this repo's own
workflow/scripts/igvf_metadata/ package, so this one file (plus the
IGVF_GET_REQUEST_USER_GUIDE.md next to it) can be copied out and sent on its
own. It mirrors (does not import) three pieces of that package, credited
here so the two can be kept in sync if the portal's own schema changes:
  - the multireport GET (query string + auth pattern):
    workflow/scripts/igvf_metadata/portal_client.py's PortalReader.get_multireport
    and workflow/scripts/igvf_metadata/cell_metadata.py's _MULTIREPORT_QUERY.
  - primary/principal classification: workflow/scripts/igvf_metadata/
    pseudobulk_sets.py's classify() (was cell_metadata._classify until
    2026-08-17). Cell_type field extraction: cell_metadata.py's
    _term_id_from_cell_type/_term_name_from_cell_type.
  - the most-contributing-subsample tie-break and its subsample-frequency
    source (only used in QC-guide mode, see below):
    subsamples.py's subsamples_by_frequency and cell_metadata.refresh_if_stale.

Never touches credentials beyond reading IGVF_API_KEY/IGVF_SECRET_KEY/
IGVF_MODE from the environment (via igvf_utils) -- it never writes them
anywhere, and this repo is public, so nothing in this folder should ever
contain a real key/secret value. See IGVF_GET_REQUEST_USER_GUIDE.md.

MODES -- one script, two behaviors, chosen by whether --qc-guide-dir is given:

  1. No --qc-guide-dir (works everywhere, for anyone with an IGVF key pair):
     each dataset-cluster's CellAnnotation/SampleTermID/SampleTermName/
     CellQualifier is every DISTINCT value seen across its contributing
     primary pseudobulks, " | "-joined -- a lone value the overwhelming
     majority of the time; more than one means those primaries disagree.
     Subsamples are listed sorted alphabetically (no contribution-order
     information available without a QC guide). QCGuideFile is blank on
     every row.

  2. --qc-guide-dir given (needs this lab's own per-cluster filtered-barcode
     QC guides, so only useful to whoever has them -- but resolves a single
     value, not a disagreement report): for each dataset-cluster, look for a
     QC guide under {qc-guide-dir}/{dataset}/{cluster}/ and, if found (and
     that cluster's contributing primaries agree on term_id/term_name --
     they're supposed to always agree, so disagreement there is treated the
     same as no guide), resolve CellAnnotation/SampleTermID/SampleTermName/
     CellQualifier to the single value from the subsample contributing the
     most cells in that guide, and order Subsamples by descending
     contribution (only subsamples contributing >=1 cell). QCGuideFile is
     that guide file's name. Any dataset-cluster with no matching guide (or
     with contributing primaries too inconsistent to trust) falls back to
     mode 1's behavior for that one row instead of failing the whole run --
     QCGuideFile stays blank for those rows, same as mode 1, so which mode
     actually applied to a given row is always visible at a glance.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import Counter

# Same as cell_metadata._MULTIREPORT_QUERY (workflow/scripts/igvf_metadata/cell_metadata.py):
# limit=all is required -- the multireport endpoint otherwise silently caps at a default
# page size (25) instead of returning every PseudobulkSet.
MULTIREPORT_QUERY = (
    "type=PseudobulkSet&status%21=deleted&limit=all"
    "&field=%40id&field=cell_annotation&field=aliases&field=cell_type"
    "&field=cell_type.term_name&field=cell_type.term_id&field=summary&field=cell_qualifier"
    "&field=input_file_sets&field=lab&field=samples&field=status"
)

# Joins multi-value columns -- NOT a comma: cell_annotation/term_name are free-text
# portal fields that routinely contain a literal comma themselves (confirmed against
# real production data, e.g. cell_type.term_name "B cell, CD19-positive") -- joining
# distinct values with "," would make a single, consistent value indistinguishable
# from two disagreeing ones. " | " is not expected in any of these fields.
JOIN = " | "

OUTPUT_COLUMNS = [
    "Dataset_Cluster",
    "Subsamples",
    "N_Subsamples",
    "CellAnnotation",
    "SampleTermID",
    "SampleTermName",
    "CellQualifier",
    "Status",
    "QCGuideFile",
]


def log(msg):
    print(f"[igvf_cell_annotation_report] {msg}", file=sys.stderr)


def fetch_multireport(igvf_mode):
    """One raw GET against /multireport/ -- mirrors
    portal_client.PortalReader.get_multireport exactly (auth/timeout/headers),
    without copying that module's verify=False-avoidance comment out of
    context: this also never disables TLS verification."""
    import requests
    import igvf_utils as iu
    import igvf_utils.utils as iuu
    from igvf_utils.connection import Connection

    # no_log_file=True: avoids igvf_utils dropping an IU_Logs/ directory wherever
    # a recipient of this standalone script happens to run it.
    conn = Connection(igvf_mode=igvf_mode, no_log_file=True)
    url = iuu.url_join([conn.igvf_mode.url, "multireport/?"]) + MULTIREPORT_QUERY
    response = requests.get(url, auth=conn.auth, timeout=iu.TIMEOUT, headers=iuu.REQUEST_HEADERS_JSON)
    response.raise_for_status()
    return response.json()["@graph"]


def _classify(row):
    """"primary" (input_file_sets all AnalysisSets), "principal" (all
    PseudobulkSets), or None (ambiguous/unparseable). Mirrors
    igvf_metadata/pseudobulk_sets.py's classify() -- input_file_sets entries never carry an
    "@type" sub-field, confirmed against a real production call, so
    classification is by "@id" path prefix instead."""
    ids = [entry.get("@id", "") for entry in row.get("input_file_sets") or [] if isinstance(entry, dict)]
    if not ids:
        return None
    if all(i.startswith("/analysis-sets/") for i in ids):
        return "primary"
    if all(i.startswith("/pseudobulk-sets/") for i in ids):
        return "principal"
    return None


def _term_id_from_cell_type(cell_type):
    """The portal's own CURIE-style term_id (e.g. "CL:0000235"). Mirrors
    cell_metadata._term_id_from_cell_type."""
    if not isinstance(cell_type, dict):
        return None
    return cell_type.get("term_id") or None


def _term_name_from_cell_type(cell_type):
    """The human-readable cell_type name (e.g. "macrophage"). Mirrors
    cell_metadata._term_name_from_cell_type."""
    if not isinstance(cell_type, dict):
        return None
    return cell_type.get("term_name") or None


def _dataset_cluster(alias, subsample):
    """"{lab}:{dataset}-{cluster}-{subsample}" -> "{dataset}-{cluster}" --
    same derivation as the ad hoc SQL used earlier this session
    (substr(alias, 1, length(alias) - length(subsample) - 1)), computed
    per-row here instead of via a SQL table, so no local cluster config is
    needed to know which (dataset, cluster) pairs exist."""
    suffix = alias.split(":", 1)[-1] if ":" in alias else alias
    tail = f"-{subsample}"
    if not suffix.endswith(tail):
        return None
    return suffix[: -len(tail)]


def split_dataset_cluster(dataset_cluster):
    """"{dataset}-{cluster}" -> (dataset, cluster), splitting on the FIRST
    "-". Safe because every real dataset name observed on the portal is a
    plain "igvfN" token with no internal hyphen (e.g. "igvf10") -- the
    hyphen immediately after it is the one this splits on. Returns (None,
    None) if there's no hyphen at all (shouldn't happen for anything
    collect_primary_rows produced, but callers should treat that as "no
    QC guide lookup possible" rather than raising)."""
    if "-" not in dataset_cluster:
        return None, None
    dataset, cluster = dataset_cluster.split("-", 1)
    return dataset, cluster


def collect_primary_rows(rows):
    """Every primary pseudobulk row with exactly one contributing sample and
    a usable alias -- same "exactly 1 contributing sample" filter as
    cell_metadata.py's module docstring: "not a data-quality heuristic --
    it's the actual filter for pipeline membership." Returns
    {dataset_cluster: [row, ...]}."""
    by_group = {}
    skipped_multi_sample = 0
    skipped_no_alias = 0
    skipped_unparseable_alias = 0
    for row in rows:
        if _classify(row) != "primary":
            continue
        sample_accessions = [s.get("accession") for s in row.get("samples") or [] if isinstance(s, dict)]
        if len(sample_accessions) != 1:
            skipped_multi_sample += 1
            continue
        subsample = sample_accessions[0]
        aliases = row.get("aliases") or []
        alias = aliases[0] if aliases else None
        if not alias:
            skipped_no_alias += 1
            continue
        dataset_cluster = _dataset_cluster(alias, subsample)
        if dataset_cluster is None:
            skipped_unparseable_alias += 1
            continue
        by_group.setdefault(dataset_cluster, []).append(
            {
                "subsample": subsample,
                "cell_annotation": row.get("cell_annotation"),
                "term_id": _term_id_from_cell_type(row.get("cell_type")),
                "term_name": _term_name_from_cell_type(row.get("cell_type")),
                "cell_qualifier": row.get("cell_qualifier"),
                "status": row.get("status"),
            }
        )
    log(
        f"{len(rows)} total portal row(s) -> {sum(len(v) for v in by_group.values())} primary pseudobulk(s) "
        f"across {len(by_group)} dataset-cluster group(s) ({skipped_multi_sample} skipped for having != 1 "
        f"contributing sample, {skipped_no_alias} for lacking an alias, {skipped_unparseable_alias} for an "
        "alias not ending in \"-{subsample}\")"
    )
    return by_group


def distinct_values(values):
    """The distinct non-empty value set -- disagreement is judged on this
    set's size, never by scanning a joined string (see JOIN)."""
    return sorted({v for v in values if v})


def _status_of(members):
    return "released" if all(m["status"] == "released" for m in members) else "in progress"


def build_fallback_row(dataset_cluster, members):
    """Mode 1's behavior for one dataset-cluster: every distinct value,
    " | "-joined, no contribution ordering, QCGuideFile blank. Used for
    every row when --qc-guide-dir isn't given at all, and per-row in
    QC-guide mode whenever that cluster has no usable guide."""
    return {
        "Dataset_Cluster": dataset_cluster,
        "Subsamples": JOIN.join(sorted({m["subsample"] for m in members})),
        "N_Subsamples": len({m["subsample"] for m in members}),
        "CellAnnotation": JOIN.join(distinct_values(m["cell_annotation"] for m in members)),
        "SampleTermID": JOIN.join(distinct_values(m["term_id"] for m in members)),
        "SampleTermName": JOIN.join(distinct_values(m["term_name"] for m in members)),
        "CellQualifier": JOIN.join(distinct_values(m["cell_qualifier"] for m in members)),
        "Status": _status_of(members),
        "QCGuideFile": "",
    }


# A real per-cluster QC output directory (e.g.
# QC_pseudobulks/plots/{dataset}/{cluster}/) routinely has multiple .tsv/.tsv.gz
# files that are NOT the per-barcode guide -- filtered_cell_subsample_metrics.tsv
# (a per-subsample summary) and qc_thresholds.tsv (threshold config) sit right next
# to the real guide, confirmed against every real cluster directory on Oak. Only
# this exact filename has the one-row-per-barcode "subsample" column
# subsample_counts_from_guide needs -- picking anything else would silently
# undercount or KeyError. Always this exact name, always gzipped (per user, 2026-07-31).
_DEFAULT_QC_GUIDE_NAME = "filtered_barcodes_with_subsamples.tsv.gz"


def find_qc_guide(qc_guide_dir, dataset, cluster):
    """Looks in {qc_guide_dir}/{dataset}/{cluster}/ for the guide file.
    Returns (path, filename), or (None, None) if the directory doesn't
    exist or no guide can be identified unambiguously -- callers must treat
    that as "fall back to mode 1 for this row," never as an error (the
    whole point of this mode is that the lab has only shared guides for
    some clusters, under a directory naming convention that isn't
    guaranteed to match the portal's own cluster name for every cluster).

    Preference order: the well-known default filename if present (this is
    every real cluster observed so far); else the sole .tsv/.tsv.gz file in
    the directory, if there's exactly one (unambiguous even under an
    unfamiliar naming convention); else this is genuinely ambiguous --
    don't guess, fall back instead."""
    if dataset is None or cluster is None:
        return None, None
    candidate_dir = os.path.join(qc_guide_dir, dataset, cluster)
    if not os.path.isdir(candidate_dir):
        return None, None
    entries = set(os.listdir(candidate_dir))
    if _DEFAULT_QC_GUIDE_NAME in entries:
        return os.path.join(candidate_dir, _DEFAULT_QC_GUIDE_NAME), _DEFAULT_QC_GUIDE_NAME
    candidates = sorted(f for f in entries if f.endswith(".tsv") or f.endswith(".tsv.gz"))
    if len(candidates) == 1:
        return os.path.join(candidate_dir, candidates[0]), candidates[0]
    if len(candidates) > 1:
        log(
            f"{dataset}/{cluster}: {len(candidates)} .tsv/.tsv.gz file(s) in {candidate_dir}, none named "
            f"{_DEFAULT_QC_GUIDE_NAME} -- ambiguous which is the real guide, falling back instead of guessing"
        )
    return None, None


def subsample_counts_from_guide(path):
    """{subsample: contributing cell count}, from the QC guide's own
    "subsample" column (one row per barcode/cell) -- mirrors
    subsamples.py's subsamples_by_frequency exactly, just inlined here so
    this file has no import on that package."""
    opener = gzip.open if path.endswith(".gz") else open
    counts = Counter()
    with opener(path, "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            counts[row["subsample"]] += 1
    return counts


def resolve_group(dataset_cluster, members, qc_guide_dir):
    """One dataset-cluster's output row. Falls back to
    build_fallback_row (mode 1) for anything that would otherwise require
    guessing: no guide found, an empty/unreadable guide, term_id/term_name
    disagreement across contributing primaries (should never happen; treated
    like a missing guide rather than trusted), or the guide's own
    top-contributing subsample having no matching primary pseudobulk on the
    portal yet."""
    dataset, cluster = split_dataset_cluster(dataset_cluster)
    guide_path, guide_name = find_qc_guide(qc_guide_dir, dataset, cluster)
    if guide_path is None:
        return build_fallback_row(dataset_cluster, members)

    counts = subsample_counts_from_guide(guide_path)
    # "only the ones that contribute at least 1 cell" -- trivially true for anything
    # appearing as a row in the guide at all, kept explicit for defensiveness.
    ordered = [s for s, c in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])) if c >= 1]
    if not ordered:
        log(f"{dataset}/{cluster}: QC guide {guide_name} has no contributing subsamples -- falling back")
        return build_fallback_row(dataset_cluster, members)

    term_ids = distinct_values(m["term_id"] for m in members)
    term_names = distinct_values(m["term_name"] for m in members)
    if len(term_ids) > 1 or len(term_names) > 1:
        log(
            f"WARNING {dataset}/{cluster}: {len(term_ids)} distinct SampleTermID / {len(term_names)} distinct "
            "SampleTermName across its primary pseudobulks -- this should never happen -- not trusting the "
            f"QC guide {guide_name} for this row, falling back"
        )
        return build_fallback_row(dataset_cluster, members)

    winning_subsample = ordered[0]
    winning_member = next((m for m in members if m["subsample"] == winning_subsample), None)
    if winning_member is None:
        log(
            f"{dataset}/{cluster}: QC guide {guide_name}'s top-contributing subsample ({winning_subsample}) has "
            "no matching primary pseudobulk on the portal yet -- falling back"
        )
        return build_fallback_row(dataset_cluster, members)

    return {
        "Dataset_Cluster": dataset_cluster,
        "Subsamples": JOIN.join(ordered),
        "N_Subsamples": len(ordered),
        "CellAnnotation": winning_member["cell_annotation"] or "",
        "SampleTermID": term_ids[0] if term_ids else "",
        "SampleTermName": term_names[0] if term_names else "",
        "CellQualifier": winning_member["cell_qualifier"] or "",
        "Status": _status_of(members),
        "QCGuideFile": guide_name,
    }


def build_report_rows(by_group, qc_guide_dir):
    report_rows = []
    resolved_count = 0
    for dataset_cluster, members in sorted(by_group.items()):
        if qc_guide_dir:
            row = resolve_group(dataset_cluster, members, qc_guide_dir)
        else:
            row = build_fallback_row(dataset_cluster, members)
        if row["QCGuideFile"]:
            resolved_count += 1
        report_rows.append(row)
    inconsistent = sum(1 for r in report_rows if JOIN in r["CellAnnotation"])
    log(
        f"{len(report_rows)} dataset-cluster row(s) written; {resolved_count} resolved via a QC guide "
        f"(QCGuideFile set); {inconsistent} have a multi-value CellAnnotation (disagreement, not resolved)"
    )
    return report_rows


def write_tsv(rows, path):
    # lineterminator="\n": csv's default "excel" dialect writes "\r\n" regardless of
    # platform, which is valid TSV but trips up naive `awk`/`grep -c`/etc. on the last
    # column (confirmed while verifying this script -- a genuinely blank last-column
    # value looked non-empty to awk because of the trailing "\r").
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(OUTPUT_COLUMNS)
        for row in rows:
            writer.writerow([row.get(c, "") for c in OUTPUT_COLUMNS])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "-o", "--output", default="cell_annotations_by_dataset_cluster.tsv", help="output TSV path"
    )
    parser.add_argument(
        "--igvf-mode",
        default="prod",
        help="prod/staging/sandbox -- defaults to prod (the real production portal)",
    )
    parser.add_argument(
        "--qc-guide-dir",
        default=None,
        help=(
            "Optional. Root of a directory structured {dataset}/{cluster}/{qc_guide_file}.tsv[.gz] "
            "(this lab's per-cluster filtered-barcode QC guides). If given, resolves a single "
            "CellAnnotation/SampleTermID/SampleTermName/CellQualifier per dataset-cluster and orders "
            "Subsamples by contribution wherever a matching guide is found; every other dataset-cluster "
            "(and every row at all, if this isn't given) falls back to listing every distinct value "
            "instead. Never required, never causes a failure by itself."
        ),
    )
    args = parser.parse_args()

    log(
        "mode: "
        + (
            f"QC-guide-resolved where available (--qc-guide-dir={args.qc_guide_dir})"
            if args.qc_guide_dir
            else "report-all (no --qc-guide-dir given)"
        )
    )
    log("issuing multireport GET against the IGVF Data Portal...")
    rows = fetch_multireport(args.igvf_mode)
    by_group = collect_primary_rows(rows)
    report_rows = build_report_rows(by_group, args.qc_guide_dir)
    write_tsv(report_rows, args.output)
    log(f"wrote {args.output}")


if __name__ == "__main__":
    main()
