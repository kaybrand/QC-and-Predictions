"""
Writes/updates the two config artifacts scE2G and the reformat/candidates/
features rules need, ONE DATASET AT A TIME (called once per distinct dataset
among this run's included clusters -- see common.smk):

1. {dataset}_cell_clusters.tsv -- scE2G's `cell_clusters` input table for
   THAT dataset only (cluster, rna_matrix_file, atac_frag_file, HiC_*,
   alt_*, model_dir). Existing rows for other clusters in the same dataset
   are preserved; only rows for this run's clusters in this dataset are
   added/replaced (per the qc-filter-pseudobulks skill's convention of
   adding a row rather than recreating the table). Partitioning per dataset
   is also what makes cluster-name collisions across datasets a non-issue:
   scE2G only ever sees one dataset's table at a time.

2. {dataset}_cluster_metadata.tsv -- ontology_id / cell_type / summary per
   cluster in THAT dataset, joined from the team's lab_annotations_with_cl.tsv.
   `summary` is just the cluster name (the E2G Pillar Project's own
   convention -- no lookup; unique within a dataset, which is the only
   uniqueness this convention has ever required). `cell_type`/`ontology_id`
   are looked up by (dataset, lab_celltype); missing fields are written as
   "TODO: <field>" rather than silently left incomplete. Re-running this
   NEVER overwrites a value that isn't a TODO placeholder, so hand-filled
   gaps survive a later re-join once the lab annotation table catches up.
"""

import csv
import os

LAB_ANNOTATIONS_PATH = "/oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/scE2G_products_table/metadata/lab_annotations_with_cl.tsv"

CELL_CLUSTERS_HEADER = [
    "cluster", "rna_matrix_file", "atac_frag_file",
    "HiC_file", "HiC_type", "HiC_resolution",
    "alt_TSS", "alt_genes", "model_dir",
]
CLUSTER_METADATA_HEADER = ["cluster", "ontology_id", "cell_type", "summary"]


def _load_lab_annotations(path):
    """Returns {(dataset.upper(), lab_celltype.lower()): row_dict}."""
    lookup = {}
    with open(path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            key = (row["dataset"].strip().upper(), row["lab_celltype"].strip().lower())
            lookup[key] = row
    return lookup


def _resolve_cell_type_and_ontology(dataset, pseudobulk_annotation, lab_annotations):
    row = lab_annotations.get((dataset.upper(), pseudobulk_annotation.lower()))
    if row is None:
        return "TODO: cell_type (no lab_annotations_with_cl.tsv row found)", "TODO: ontology_id"

    cl_term = row.get("CL term", "").strip()
    qualifier = row.get("qualifier", "").strip()
    ontology_id = row.get("CL_ID", "").strip()

    cell_type = " ".join(part for part in (cl_term, qualifier) if part) or "TODO: cell_type (blank CL term/qualifier)"
    ontology_id = ontology_id or "TODO: ontology_id (blank CL_ID)"
    return cell_type, ontology_id


def _read_existing_tsv(path, key_col):
    rows = {}
    if os.path.exists(path):
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                rows[row[key_col]] = row
    return rows


def _write_tsv(path, header, rows_by_key):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for key in sorted(rows_by_key):
            writer.writerow(rows_by_key[key])


def write_cell_clusters_table(dataset, dataset_clusters_cfg, included_cluster_names, out_dir, out_dir_data, scE2G_dir):
    """
    dataset_clusters_cfg: config["clusters"][dataset] (bare cluster name -> cluster config)
    included_cluster_names: bare cluster names (within this dataset) included this run
    out_dir_data: the multiome_data/{dataset} directory holding filtered ATAC/RNA outputs

    model_dir entries must be ABSOLUTE paths: scE2G is imported here as a Snakemake
    `module` from a different working directory than scE2G_dir, so the relative
    "models/{name}" form that works when scE2G is run standalone from its own
    directory resolves to nothing (or the wrong thing) here.
    """
    path = os.path.join(out_dir, f"{dataset}_cell_clusters.tsv")
    existing = _read_existing_tsv(path, "cluster")

    for cluster in included_cluster_names:
        cluster_cfg = dataset_clusters_cfg[cluster]
        is_atac_only = cluster_cfg["models"] == ["scATAC_powerlaw_v3"]
        atac_frag_file = os.path.join(out_dir_data, cluster, f"atac_fragments_{dataset}_{cluster}.tsv.gz")
        rna_matrix_file = "" if is_atac_only else os.path.join(out_dir_data, cluster, f"rna_count_matrix_{dataset}_{cluster}")
        existing[cluster] = {
            "cluster": cluster,
            "rna_matrix_file": rna_matrix_file,
            "atac_frag_file": atac_frag_file,
            "HiC_file": "", "HiC_type": "", "HiC_resolution": "",
            "alt_TSS": "", "alt_genes": "",
            "model_dir": ",".join(os.path.join(scE2G_dir, "models", m) for m in cluster_cfg["models"]),
        }

    _write_tsv(path, CELL_CLUSTERS_HEADER, existing)
    return path


def write_cluster_metadata_table(dataset, dataset_clusters_cfg, included_cluster_names, out_dir):
    path = os.path.join(out_dir, f"{dataset}_cluster_metadata.tsv")
    existing = _read_existing_tsv(path, "cluster")
    lab_annotations = _load_lab_annotations(LAB_ANNOTATIONS_PATH)

    for cluster in included_cluster_names:
        cluster_cfg = dataset_clusters_cfg[cluster]
        current = existing.get(cluster, {})
        # Never clobber a manually-filled (non-TODO) value with a fresh join result.
        if current.get("cell_type", "TODO").startswith("TODO") or "cell_type" not in current:
            cell_type, ontology_id = _resolve_cell_type_and_ontology(
                dataset, cluster_cfg["pseudobulk_annotation"], lab_annotations
            )
        else:
            cell_type, ontology_id = current["cell_type"], current["ontology_id"]

        existing[cluster] = {
            "cluster": cluster,
            "ontology_id": ontology_id,
            "cell_type": cell_type,
            "summary": cluster,  # E2G Pillar Project convention: summary == cluster, unique within a dataset
        }

    _write_tsv(path, CLUSTER_METADATA_HEADER, existing)
    return path


def load_cluster_metadata(path):
    """Returns {cluster: {ontology_id, cell_type, summary}} for one dataset's table."""
    return _read_existing_tsv(path, "cluster")
