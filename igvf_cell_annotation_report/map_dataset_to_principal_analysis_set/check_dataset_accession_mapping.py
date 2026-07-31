#!/usr/bin/env python3
"""Exploratory, one-off: check whether "dataset" (parsed from a primary
PseudobulkSet's alias, "{lab}:{dataset}-{cluster}-{subsample}") maps 1:1 to
input_file_sets accession, or whether one dataset spans multiple distinct
input_file_sets accessions (relevant to the {dataset}-{cluster}-{subsample}
-> {principal analysis set accession}-{term name}-{subsample} aliasing
migration).

Auth/GET pattern mirrors workflow/scripts/igvf_metadata/portal_client.py's
PortalReader.get_multireport + cell_metadata.py's _MULTIREPORT_QUERY -- never
touches credentials beyond letting igvf_utils.connection.Connection read
IGVF_API_KEY/IGVF_SECRET_KEY/IGVF_MODE from the environment itself.
"""

import argparse
import json
import sys
from collections import defaultdict

QUERY = (
    "type=PseudobulkSet&status%21=deleted&limit=all"
    "&field=%40id&field=aliases&field=input_file_sets"
    "&field=input_file_sets.%40id&field=input_file_sets.accession&field=samples"
)


def log(msg):
    print(f"[check_dataset_accession_mapping] {msg}", file=sys.stderr)


def fetch(igvf_mode=None):
    import requests
    import igvf_utils as iu
    import igvf_utils.utils as iuu
    from igvf_utils.connection import Connection

    conn = Connection(igvf_mode=igvf_mode, no_log_file=True)
    url = iuu.url_join([conn.igvf_mode.url, "multireport/?"]) + QUERY
    response = requests.get(url, auth=conn.auth, timeout=iu.TIMEOUT, headers=iuu.REQUEST_HEADERS_JSON)
    response.raise_for_status()
    return response.json()["@graph"]


def classify(row):
    ids = [e.get("@id", "") for e in row.get("input_file_sets") or [] if isinstance(e, dict)]
    if not ids:
        return None
    if all(i.startswith("/analysis-sets/") for i in ids):
        return "primary"
    if all(i.startswith("/pseudobulk-sets/") for i in ids):
        return "principal"
    return None


def dataset_from_alias(alias):
    """"{lab}:{dataset}-{cluster}-{subsample}" -> dataset. Every real
    dataset name on the portal is a plain "igvfN" token with no internal
    hyphen (confirmed in cell_metadata.py's split_dataset_cluster comment),
    so the first "-"-delimited token of the alias suffix is the dataset."""
    suffix = alias.split(":", 1)[-1] if ":" in alias else alias
    return suffix.split("-", 1)[0] if "-" in suffix else None


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--json-out",
        default=None,
        help=(
            "Optional path to write the dataset -> input_file_sets.accession (i.e. principal analysis set "
            "accession(s)) mapping as JSON, one entry per dataset with data available. Value is a sorted list "
            "of every distinct accession seen for that dataset -- a single-element list for the (expected) "
            "1:1 case, more than one element if that dataset's primaries disagree."
        ),
    )
    args = parser.parse_args()

    log("issuing multireport GET against the IGVF Data Portal...")
    rows = fetch()
    log(f"{len(rows)} total PseudobulkSet row(s) returned")

    dataset_to_accessions = defaultdict(set)
    dataset_to_dataset_clusters = defaultdict(set)
    skipped_not_primary = 0
    skipped_no_alias = 0
    skipped_multi_input_file_set = 0
    skipped_unparseable_dataset = 0

    for row in rows:
        if classify(row) != "primary":
            skipped_not_primary += 1
            continue
        aliases = row.get("aliases") or []
        alias = aliases[0] if aliases else None
        if not alias:
            skipped_no_alias += 1
            continue
        input_file_sets = [e for e in row.get("input_file_sets") or [] if isinstance(e, dict)]
        accessions = {e.get("accession") for e in input_file_sets if e.get("accession")}
        if len(accessions) != 1:
            skipped_multi_input_file_set += 1
            continue
        accession = next(iter(accessions))
        dataset = dataset_from_alias(alias)
        if dataset is None:
            skipped_unparseable_dataset += 1
            continue
        dataset_to_accessions[dataset].add(accession)
        suffix = alias.split(":", 1)[-1] if ":" in alias else alias
        sample_accessions = [s.get("accession") for s in row.get("samples") or [] if isinstance(s, dict)]
        subsample = sample_accessions[0] if len(sample_accessions) == 1 else None
        dataset_cluster = None
        if subsample and suffix.endswith(f"-{subsample}"):
            dataset_cluster = suffix[: -len(f"-{subsample}")]
        dataset_to_dataset_clusters[dataset].add(dataset_cluster or suffix)

    log(
        f"{skipped_not_primary} skipped (not classified as primary), {skipped_no_alias} skipped (no alias), "
        f"{skipped_multi_input_file_set} skipped (!=1 distinct input_file_sets accession), "
        f"{skipped_unparseable_dataset} skipped (unparseable dataset from alias)"
    )

    n_datasets = len(dataset_to_accessions)
    total_dataset_clusters = sum(len(v) for v in dataset_to_dataset_clusters.values())
    log(f"{n_datasets} distinct dataset(s) across {total_dataset_clusters} dataset-cluster group(s)")

    counts = {}
    for dataset, accessions in sorted(dataset_to_accessions.items()):
        counts[len(accessions)] = counts.get(len(accessions), 0) + 1

    print("\n=== Distribution: number of datasets by count of distinct input_file_sets accessions ===")
    for n_acc, n_ds in sorted(counts.items()):
        print(f"  {n_acc} distinct accession(s): {n_ds} dataset(s)")

    multi = {d: a for d, a in dataset_to_accessions.items() if len(a) > 1}
    print(f"\n=== 1:1 datasets: {n_datasets - len(multi)} / {n_datasets} ===")
    if multi:
        offenders = ", ".join(f"{d} ({sorted(a)})" for d, a in sorted(multi.items()))
        log(
            f"WARNING: {len(multi)}/{n_datasets} dataset(s) break the assumed dataset <-> principal analysis "
            f"set accession 1:1 mapping: {offenders} -- all distinct accessions are still listed for them in "
            "the JSON output below, but downstream consumers expecting a single accession per dataset should "
            "not silently pick one."
        )
        print(f"\n=== Datasets with >1 distinct input_file_sets accession ({len(multi)}) ===")
        for dataset, accessions in sorted(multi.items()):
            n_clusters = len(dataset_to_dataset_clusters[dataset])
            print(f"  {dataset}: {len(accessions)} accessions ({sorted(accessions)}), {n_clusters} dataset-cluster(s)")

    # dataset -> sorted list of every distinct principal analysis set accession seen for it --
    # every dataset with data available, not just the 1:1 ones, so a future disagreement is
    # visible in the JSON itself rather than silently dropping that dataset from it.
    mapping = {d: sorted(a) for d, a in sorted(dataset_to_accessions.items())}
    print("\n=== dataset -> principal analysis set accession(s) (input_file_sets.accession) ===")
    print(json.dumps(mapping, indent=2))
    if args.json_out:
        with open(args.json_out, "w") as f:
            json.dump(mapping, f, indent=2)
        log(f"wrote {args.json_out}")


if __name__ == "__main__":
    main()
