#!/usr/bin/env python
"""Read-only reconnaissance: walks every dataset/cluster folder currently on
Synapse under filtered_data_parent_id ("Multiome datasets") and diffs against
this round's in-scope cluster set (report.tsv), to surface any
multiome_dataset entry with no replacement this round (igvf7 expected --
anything else here needs a manual second look before assuming the same).

Purely getChildren calls -- never store/delete. The user reviews the output
and decides what to do about each orphan; this script never acts on Synapse.
manage_synapse_manifest.py's own list_owned_entities() can't surface these by
design: it's scoped to exactly the cluster_keys it's told about, so it can
never discover a dataset/cluster it wasn't asked to look at.
"""

import argparse
import csv
import os
import sys

import synapseclient

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))

DEFAULT_PARENT_ID = "syn53469844"  # "Multiome datasets" -- filtered_data_parent_id


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent-id", default=DEFAULT_PARENT_ID)
    p.add_argument("--output-dir", default=os.path.join(REPO_ROOT, "results"))
    p.add_argument("--in-scope-tsv", default=None, help="default: {output_dir}/report.tsv")
    p.add_argument("--out", default=None, help="default: {output_dir}/synapse_orphans.tsv")
    return p.parse_args()


def load_in_scope(path):
    in_scope = set()
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            in_scope.add((row["dataset"], row["cluster"]))
    return in_scope


def list_synapse_dataset_clusters(syn, parent_id):
    """Every (dataset, cluster) folder pair currently on Synapse under
    parent_id -> (folder_id, n_children)."""
    found = {}
    dataset_folders = {c["name"]: c["id"] for c in syn.getChildren(parent_id, includeTypes=["folder"])}
    for dataset, dataset_folder_id in sorted(dataset_folders.items()):
        cluster_folders = {c["name"]: c["id"] for c in syn.getChildren(dataset_folder_id, includeTypes=["folder"])}
        for cluster, cluster_folder_id in sorted(cluster_folders.items()):
            n_children = sum(1 for _ in syn.getChildren(cluster_folder_id))
            found[(dataset, cluster)] = (cluster_folder_id, n_children)
    return found


def describe_orphan_files(syn, cluster_folder_id, my_user_id):
    """Per-file detail for one orphaned cluster folder: name, who last
    modified/created it, and whether that's the current user (my_user_id) --
    deleting someone else's upload is a materially different, more sensitive
    situation than deleting your own stale files, so this is surfaced
    explicitly rather than assumed."""
    files = []
    for child in syn.getChildren(cluster_folder_id, includeTypes=["file"]):
        entity = syn.get(child["id"], downloadFile=False)
        modified_by = str(entity.modifiedBy)
        created_by = str(entity.createdBy)
        files.append({
            "name": child["name"],
            "id": child["id"],
            "created_by": created_by,
            "modified_by": modified_by,
            "is_mine": modified_by == my_user_id and created_by == my_user_id,
        })
    return files


def main():
    args = parse_args()
    in_scope_tsv = args.in_scope_tsv or os.path.join(args.output_dir, "report.tsv")
    out_path = args.out or os.path.join(args.output_dir, "synapse_orphans.tsv")

    in_scope = load_in_scope(in_scope_tsv)

    syn = synapseclient.login()
    my_user_id = str(syn.getUserProfile()["ownerId"])
    on_synapse = list_synapse_dataset_clusters(syn, args.parent_id)

    orphans = sorted(set(on_synapse) - in_scope)
    orphan_files = {
        (dataset, cluster): describe_orphan_files(syn, on_synapse[(dataset, cluster)][0], my_user_id)
        for dataset, cluster in orphans
    }

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "dataset", "cluster", "synapse_folder_id", "n_children",
            "file_name", "file_id", "created_by", "modified_by", "uploaded_by_me",
        ])
        for dataset, cluster in orphans:
            folder_id, n_children = on_synapse[(dataset, cluster)]
            files = orphan_files[(dataset, cluster)]
            if not files:
                writer.writerow([dataset, cluster, folder_id, n_children, "", "", "", "", ""])
                continue
            for f_info in files:
                writer.writerow([
                    dataset, cluster, folder_id, n_children,
                    f_info["name"], f_info["id"], f_info["created_by"], f_info["modified_by"],
                    "y" if f_info["is_mine"] else "n",
                ])

    print(f"[list_synapse_orphans] logged in as user id {my_user_id}", file=sys.stderr)
    print(f"[list_synapse_orphans] {len(on_synapse)} (dataset, cluster) folder(s) found on Synapse under {args.parent_id}", file=sys.stderr)
    print(f"[list_synapse_orphans] {len(orphans)} orphan cluster(s) (no replacement in this round's scope) written to {out_path}", file=sys.stderr)
    if orphans:
        print(f"[list_synapse_orphans] orphans: {orphans}", file=sys.stderr)
        not_mine = sorted({
            (ds, cl) for (ds, cl), files in orphan_files.items()
            if any(not f_info["is_mine"] for f_info in files)
        })
        if not_mine:
            print(
                f"[list_synapse_orphans] WARNING: {len(not_mine)} orphan cluster(s) contain file(s) NOT "
                f"created+modified by you (user {my_user_id}) -- review before any deletion: {not_mine}",
                file=sys.stderr,
            )


if __name__ == "__main__":
    main()
