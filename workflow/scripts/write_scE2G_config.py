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

Concurrency: under --executor slurm, EVERY worker node independently
re-parses this pipeline's whole Snakefile (confirmed empirically -- each
job's own log shows a fresh "Building DAG of jobs..."), and several jobs run
concurrently. That means these two tables get read-merged-written by
multiple processes at once. A plain read-then-truncate-then-write is not
safe under that: one process's write (which truncates the file immediately
on open) can be observed mid-write by another process's read, producing a
torn/incomplete file and a KeyError downstream. _locked_read_merge_write
guards against this with both an advisory flock (primary defense, serializes
the whole read-merge-write) and a write-to-temp-then-os.replace (atomic
rename, so even if a lock is ever not honored -- e.g. an unexpected NFS
mount -- no reader can ever observe a half-written file).
"""

import csv
import fcntl
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


def _read_tsv_rows(path, key_col):
    rows = {}
    if os.path.exists(path):
        with open(path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                rows[row[key_col]] = row
    return rows


def _atomic_write_tsv(path, header, rows_by_key):
    """Write to a same-directory temp file, then atomically replace the target.
    Guarantees any concurrent reader sees either the fully-old or fully-new
    file, never a partial one (os.replace is atomic on the same filesystem)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp_path = f"{path}.tmp.{os.getpid()}"
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for key in sorted(rows_by_key):
            writer.writerow(rows_by_key[key])
    os.replace(tmp_path, path)


def _locked_read_merge_write(path, header, key_col, merge_fn):
    """Read existing rows, let merge_fn add/update rows in place, write back --
    the whole read-merge-write happens under an exclusive advisory lock on a
    sibling `.lock` file, so concurrent worker-node processes serialize rather
    than interleave. merge_fn(existing_rows: dict) is called with the freshly
    read rows and must mutate/return the dict to write.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    lock_path = f"{path}.lock"
    lock_fd = os.open(lock_path, os.O_CREAT | os.O_RDWR)
    try:
        fcntl.flock(lock_fd, fcntl.LOCK_EX)
        existing = _read_tsv_rows(path, key_col)
        merged = merge_fn(existing)
        _atomic_write_tsv(path, header, merged)
    finally:
        fcntl.flock(lock_fd, fcntl.LOCK_UN)
        os.close(lock_fd)


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

    def merge(existing):
        for cluster in included_cluster_names:
            cluster_cfg = dataset_clusters_cfg[cluster]
            is_atac_only = cluster_cfg["models"] == ["scATAC_powerlaw_v3"]
            if "atac_frag_file" in cluster_cfg:
                # Already-filtered/QC'd fragments (e.g. catlas) -- use the literal
                # path directly, bypassing atac_fragment_file/filter_atac_fragments.py
                # entirely (that rule is simply never requested for these clusters).
                atac_frag_file = cluster_cfg["atac_frag_file"]
            else:
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
        return existing

    _locked_read_merge_write(path, CELL_CLUSTERS_HEADER, "cluster", merge)
    return path


def write_cluster_metadata_table(dataset, dataset_clusters_cfg, included_cluster_names, out_dir):
    path = os.path.join(out_dir, f"{dataset}_cluster_metadata.tsv")
    lab_annotations = _load_lab_annotations(LAB_ANNOTATIONS_PATH)

    def merge(existing):
        for cluster in included_cluster_names:
            cluster_cfg = dataset_clusters_cfg[cluster]
            current = existing.get(cluster, {})
            # Never clobber a manually-filled (non-TODO) value with a fresh join result.
            if current.get("cell_type", "TODO").startswith("TODO") or "cell_type" not in current:
                # cluster_cfg.get(..., cluster): clusters with no pseudobulk_annotation
                # (e.g. catlas's atac_frag_file clusters) fall through to the existing
                # "no lab_annotations_with_cl.tsv row found" TODO placeholder below,
                # keyed on the bare cluster name instead of erroring.
                cell_type, ontology_id = _resolve_cell_type_and_ontology(
                    dataset, cluster_cfg.get("pseudobulk_annotation", cluster), lab_annotations
                )
            else:
                cell_type, ontology_id = current["cell_type"], current["ontology_id"]

            existing[cluster] = {
                "cluster": cluster,
                "ontology_id": ontology_id,
                "cell_type": cell_type,
                "summary": cluster,  # E2G Pillar Project convention: summary == cluster, unique within a dataset
            }
        return existing

    _locked_read_merge_write(path, CLUSTER_METADATA_HEADER, "cluster", merge)
    return path


def load_cluster_metadata(path):
    """Returns {cluster: {ontology_id, cell_type, summary}} for one dataset's table.
    Read-only, no lock needed: writers only ever atomically replace this file, so a
    concurrent reader sees either the fully-old or fully-new version, never a torn one."""
    return _read_tsv_rows(path, "cluster")
