#!/usr/bin/env python
"""Read-only reconnaissance: finds every FILE currently on Synapse under
--parent-id (default syn53469845, "predictions") that was created by the
current user, excluding whole top-level folders named in --exclude-top, and
that would NOT be covered (new or overwrite) by a predictions manifest TSV
in the path/parent format manage_synapse_manifest.py writes.

File-level, unlike list_synapse_orphans.py's cluster-level diff against
filtered_data_parent_id -- a cluster folder can be partially covered (e.g. a
threshold value changed name) without being a whole orphaned cluster, and
this pipeline's "predictions" space is collaborative (~10 other models'
files coexist per cluster folder), so orphan status has to be checked file by
file, not folder by folder.

Two-phase to avoid an expensive syn.get() per file: a cheap getChildren-only
inventory pass first (all files under every non-excluded top-level folder),
then per-file syn.get() (for createdBy) only on the residual not covered by
the manifest.

Purely getChildren/get calls -- never store/delete. The user reviews the
output and decides what to do about each entry; this script never acts on
Synapse.
"""

import argparse
import csv
import os
import sys

import synapseclient

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

DEFAULT_PARENT_ID = "syn53469845"  # predictions_parent_id
DEFAULT_MANIFEST = os.path.join(
    REPO_ROOT, "results", "uniformly_processed", "synapse_manifests", "predictions_manifest.tsv"
)
# Top-level folders under DEFAULT_PARENT_ID that are never part of this
# pipeline's own dataset/cluster naming scheme, or that this checkout simply
# doesn't have local outputs for yet -- see git history/PR discussion for why
# each one is here rather than hardcoding an assumption in the diff logic.
DEFAULT_EXCLUDE_TOP = "(Archived)Y2_versions,GM12878_10XMultiome,K562_ERR9847049_Multiome,igvf7,kasowski"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--parent-id", default=DEFAULT_PARENT_ID)
    p.add_argument("--manifest", default=DEFAULT_MANIFEST, help="path/parent manifest TSV to diff coverage against")
    p.add_argument(
        "--exclude-top",
        default=DEFAULT_EXCLUDE_TOP,
        help="comma-separated top-level folder names under --parent-id to skip entirely",
    )
    p.add_argument(
        "--only",
        default="",
        help="comma-separated 'dataset/cluster' tokens -- if set, report only these regardless of "
        "--exclude-top or manifest coverage (for a targeted synID lookup on specific clusters)",
    )
    return p.parse_args()


def load_covered(manifest_path):
    covered = set()
    with open(manifest_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            parts = row["path"].rstrip("/").split("/")
            filename, cluster, dataset = parts[-1], parts[-2], parts[-3]
            covered.add((dataset, cluster, filename))
    return covered


def build_inventory(syn, parent_id, exclude_top, only):
    """(dataset, cluster, filename) -> {"id", "modifiedBy"} for every file under
    parent_id, skipping excluded top-level folders unless `only` is set (in
    which case exclude_top is ignored and just the requested clusters are scanned)."""
    inventory = {}
    for tf in syn.getChildren(parent_id, includeTypes=["file"]):
        inventory[("<root>", "<root>", tf["name"])] = {"id": tf["id"], "modifiedBy": tf.get("modifiedBy")}
    for tc in syn.getChildren(parent_id, includeTypes=["folder"]):
        dataset = tc["name"]
        if only:
            if not any(d == dataset for d, _ in only):
                continue
        elif dataset in exclude_top:
            continue
        for cc in syn.getChildren(tc["id"], includeTypes=["folder"]):
            cluster = cc["name"]
            if only and (dataset, cluster) not in only:
                continue
            for fc in syn.getChildren(cc["id"], includeTypes=["file"]):
                inventory[(dataset, cluster, fc["name"])] = {"id": fc["id"], "modifiedBy": fc.get("modifiedBy")}
    return inventory


def main():
    args = parse_args()
    exclude_top = {t for t in args.exclude_top.split(",") if t}
    only = set()
    for token in args.only.split(","):
        if not token:
            continue
        dataset, _, cluster = token.partition("/")
        only.add((dataset, cluster))

    syn = synapseclient.login()
    my_id = str(syn.getUserProfile()["ownerId"])

    covered = load_covered(args.manifest) if not only else set()
    inventory = build_inventory(syn, args.parent_id, exclude_top, only)
    not_covered = sorted(set(inventory) - covered)

    mine_not_covered = []
    for key in not_covered:
        meta = inventory[key]
        entity = syn.get(meta["id"], downloadFile=False)
        created_by = str(entity.createdBy)
        modified_by = str(entity.modifiedBy)
        if created_by == my_id:
            mine_not_covered.append((*key, meta["id"], created_by, modified_by))

    print(
        f"[list_synapse_prediction_orphans] user {my_id}; manifest covers {len(covered)} keys; "
        f"{len(inventory)} file(s) in scope under {args.parent_id}"
        + ("" if only else f" (excluding {sorted(exclude_top)})")
        + f"; {len(not_covered)} not covered; {len(mine_not_covered)} of those created by you",
        file=sys.stderr,
    )
    print("dataset\tcluster\tfilename\tsynapse_id\tcreated_by\tmodified_by")
    for dataset, cluster, filename, fid, created_by, modified_by in mine_not_covered:
        print(f"{dataset}\t{cluster}\t{filename}\t{fid}\t{created_by}\t{modified_by}")


if __name__ == "__main__":
    main()
