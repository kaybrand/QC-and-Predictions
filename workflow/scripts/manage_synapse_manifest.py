"""
Live-Synapse-scoped manifest management for one product type (filtered_data /
predictions / candidates / features), across one or more datasets in a
single run.

For each run: query Synapse itself (not a locally cached manifest) for
whatever this pipeline already owns under the target parent, scoped so it
can NEVER see or touch another team's or another model's files in a
collaborative folder. Diff against the fresh local file set for this run's
upload-eligible (dataset, cluster) pairs. Anything missing goes in a
path/parent manifest TSV for synapseutils.syncToSynapse(); anything stale (a
cluster that was reprocessed under a new filename, e.g. a model's
score_threshold changed, or a cluster that just became excluded) is reported
loudly and only actually deleted with --confirm-delete. Overwriting an
EXISTING entity's content is reported loudly and, for the shared "Multiome
datasets" folder specifically, requires --confirm-overwrite before syncing
for real.

Scoping rules (see the design plan's "Synapse manifests" section):
  - nested (predictions/features/filtered_data): only look under {parent}/{dataset}/{cluster}/...
  - flat (candidates): only look at children of {parent} whose name starts with "{dataset}_{cluster}_"
  - collaborative (predictions/features): additionally require the filename to
    contain one of --model-tokens, so files from the ~10 other models sharing
    that folder are invisible to this script entirely.

Cluster identity is always the (dataset, cluster) pair -- passed on the CLI
as "dataset/cluster" tokens, since cluster names are only unique within
their own dataset.
"""

import argparse
import csv
import os
import sys

import synapseclient
import synapseutils


def ensure_folder_path(syn, parent_id, path_parts, dry_run):
    """Walk/create each level of path_parts under parent_id, returning the final folder id.
    In dry-run mode, missing folders are not created -- callers should only need the id
    for entities that will actually be synced for real."""
    current = parent_id
    for part in path_parts:
        children = {c["name"]: c["id"] for c in syn.getChildren(current, includeTypes=["folder"])}
        if part in children:
            current = children[part]
        elif dry_run:
            current = f"<would-create:{part}>"
        else:
            folder = synapseclient.Folder(name=part, parent=current)
            current = syn.store(folder).id
    return current


def list_owned_entities(syn, parent_id, cluster_keys, nested, model_tokens):
    """Live query, scoped to exactly what this pipeline could have produced.
    cluster_keys: set of (dataset, cluster) tuples this run cares about."""
    owned = {}  # (dataset, cluster, name) -> {"id", "modifiedBy", "modifiedOn"}
    datasets = {dataset for dataset, _ in cluster_keys}

    def _matches(name):
        if model_tokens and not any(tok in name for tok in model_tokens):
            return False
        return True

    if nested:
        dataset_children = {c["name"]: c["id"] for c in syn.getChildren(parent_id, includeTypes=["folder"])}
        for dataset in datasets:
            if dataset not in dataset_children:
                continue
            dataset_folder_id = dataset_children[dataset]
            cluster_children = {c["name"]: c["id"] for c in syn.getChildren(dataset_folder_id, includeTypes=["folder"])}
            clusters_for_dataset = {c for d, c in cluster_keys if d == dataset}
            for cluster in clusters_for_dataset:
                if cluster not in cluster_children:
                    continue
                cluster_folder_id = cluster_children[cluster]
                for entity in syn.getChildren(cluster_folder_id):
                    if _matches(entity["name"]):
                        owned[(dataset, cluster, entity["name"])] = {
                            "id": entity["id"],
                            "modifiedBy": entity.get("modifiedBy"), "modifiedOn": entity.get("modifiedOn"),
                        }
    else:
        prefixes = {(dataset, cluster): f"{dataset}_{cluster}_" for dataset, cluster in cluster_keys}
        for entity in syn.getChildren(parent_id):
            for (dataset, cluster), prefix in prefixes.items():
                if entity["name"].startswith(prefix) and _matches(entity["name"]):
                    owned[(dataset, cluster, entity["name"])] = {
                        "id": entity["id"],
                        "modifiedBy": entity.get("modifiedBy"), "modifiedOn": entity.get("modifiedOn"),
                    }
                    break
    return owned


def synapse_key_for_local_path(local_path, dataset, cluster, nested):
    """The (dataset, cluster, name) key this file should have on Synapse, matching list_owned_entities' keying."""
    basename = os.path.basename(local_path)
    if nested:
        return dataset, cluster, basename
    name = basename if basename.startswith(f"{dataset}_{cluster}_") else f"{dataset}_{cluster}_{basename}"
    return dataset, cluster, name


def relative_folder_parts(local_path, dataset, cluster, nested):
    """Extra folder nesting between {dataset}/{cluster}/ and the file itself
    (e.g. an RNA matrix directory's contents), preserved on Synapse to match
    the local layout."""
    if not nested:
        return []
    marker = os.path.join(dataset, cluster) + os.sep
    marker = marker.replace("\\", "/")
    normalized = local_path.replace("\\", "/")
    idx = normalized.find(marker)
    if idx == -1:
        raise ValueError(f"{local_path} does not contain expected {marker}")
    rel_dir = os.path.dirname(normalized[idx + len(marker):])
    return [dataset, cluster] + ([p for p in rel_dir.split("/") if p] if rel_dir else [])


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--product", required=True)
    p.add_argument("--parent-id", required=True)
    p.add_argument("--cluster-keys", required=True, help='comma-separated "dataset/cluster" tokens, upload-eligible this run')
    p.add_argument("--nested", action="store_true")
    p.add_argument("--model-tokens", default="", help="comma-separated; only used for collaborative folders")
    p.add_argument("--manifest-out", required=True)
    p.add_argument("--dry-run", type=lambda s: s.lower() == "true", default=True)
    p.add_argument("--confirm-delete", type=lambda s: s.lower() == "true", default=False)
    p.add_argument("--confirm-overwrite", type=lambda s: s.lower() == "true", default=False)
    p.add_argument("should_exist_files", nargs="*")
    args = p.parse_args()

    cluster_keys = set()
    for token in args.cluster_keys.split(","):
        if not token:
            continue
        dataset, _, cluster = token.partition("/")
        cluster_keys.add((dataset, cluster))
    model_tokens = [t for t in args.model_tokens.split(",") if t]

    syn = synapseclient.login()

    owned = list_owned_entities(syn, args.parent_id, cluster_keys, args.nested, model_tokens)

    should_exist = {}  # (dataset, cluster, name) -> local_path
    for local_path in args.should_exist_files:
        match = next(
            ((d, c) for d, c in cluster_keys if f"/{d}/{c}/" in local_path.replace(os.sep, "/") or os.path.basename(local_path).startswith(f"{d}_{c}_")),
            None,
        )
        if match is None:
            print(f"[manage_synapse_manifest] WARNING: could not determine (dataset, cluster) for {local_path}, skipping", file=sys.stderr)
            continue
        dataset, cluster = match
        key = synapse_key_for_local_path(local_path, dataset, cluster, args.nested)
        should_exist[key] = local_path

    to_delete = {key: meta for key, meta in owned.items() if key not in should_exist}
    to_upload_new = {key: path for key, path in should_exist.items() if key not in owned}
    to_overwrite = {key: path for key, path in should_exist.items() if key in owned}

    print(f"[manage_synapse_manifest:{args.product}] {len(to_upload_new)} new, {len(to_overwrite)} overwrite, {len(to_delete)} stale/excluded")

    if to_delete:
        print(f"[manage_synapse_manifest:{args.product}] STALE/EXCLUDED ENTITIES (would delete):")
        for (dataset, cluster, name), meta in to_delete.items():
            print(f"    {dataset}/{cluster}/{name}  (synapse id {meta['id']}, last modified by {meta['modifiedBy']} on {meta['modifiedOn']})")
        if args.confirm_delete and not args.dry_run:
            for meta in to_delete.values():
                syn.delete(meta["id"])
        elif not args.dry_run:
            print(f"[manage_synapse_manifest:{args.product}] confirm-delete is False -- NOT deleting the above. Review and re-run with --confirm-delete true if correct.")

    if to_overwrite:
        print(f"[manage_synapse_manifest:{args.product}] WOULD OVERWRITE EXISTING CONTENT:")
        for key, local_path in to_overwrite.items():
            dataset, cluster, name = key
            meta = owned[key]
            print(f"    {dataset}/{cluster}/{name}  (existing synapse id {meta['id']}, last modified by {meta['modifiedBy']} on {meta['modifiedOn']}) <- {local_path}")
        if args.product == "filtered_data" and not args.confirm_overwrite and not args.dry_run:
            print(
                f"[manage_synapse_manifest:{args.product}] confirm-overwrite is False for the shared 'Multiome datasets' "
                "folder -- ABORTING real sync. Other team members may already consume these files; confirm the overwrite "
                "list above is correct (and consider emailing prior consumers) before re-running with --confirm-overwrite true."
            )
            sys.exit(1)

    manifest_rows = []
    for key, local_path in {**to_upload_new, **to_overwrite}.items():
        dataset, cluster, name = key
        if args.nested:
            path_parts = relative_folder_parts(local_path, dataset, cluster, args.nested)
            parent_id = ensure_folder_path(syn, args.parent_id, path_parts, args.dry_run)
        else:
            parent_id = args.parent_id
        manifest_rows.append({"path": local_path, "parent": parent_id})

    os.makedirs(os.path.dirname(args.manifest_out), exist_ok=True)
    with open(args.manifest_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["path", "parent"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    if manifest_rows and not args.dry_run:
        synapseutils.syncToSynapse(syn, manifestFile=args.manifest_out, dryRun=False)
    else:
        print(f"[manage_synapse_manifest:{args.product}] dry_run={args.dry_run} ({'no rows to sync' if not manifest_rows else 'not syncing'}) -- manifest written to {args.manifest_out} for review")


if __name__ == "__main__":
    main()
