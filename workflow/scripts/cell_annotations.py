"""Shared helper for CellAnnotation lookups.

annotation_lookup_key() is pure key-resolution logic, with no I/O of its
own: it takes a cluster's already-loaded config dict and returns the
(dataset, cluster-or-alias) tuple to look up in whichever backing store a
caller uses. Every caller in this pipeline (common.smk's
REFORMAT_ELIGIBLE_CLUSTERS, reformat.smk's portal_cell_metadata(),
generate_report.py) passes the resolved key straight to
igvf_metadata.state.get_cell_annotation() against the live state.db cache --
populated only by cell_metadata.refresh_if_stale's real IGVF Portal GET (see
manage_igvf_metadata.py). This file deliberately contains no code that reads
cell_annotations_by_dataset_cluster.tsv (a manually-refreshed PREVIEW
snapshot) -- gating reformat eligibility on that TSV once caused a real
crash (it said igvf2 had annotations while state.db was still cold, so
rule all requested a reformat rule that couldn't resolve).
"""


def annotation_lookup_key(dataset, cluster, cluster_cfg):
    """The (dataset, cluster) key to use when looking up a cluster's
    CellAnnotation, in whichever backing store -- honors a cluster's
    cell_annotation_key override. Used by ATAC-only variant clusters (e.g.
    pancreatic_delta_d32_ATAC_only), whose real CellAnnotation/SampleTermID/
    SampleTermName lives under the base (non-suffixed) name, not the
    suffixed cluster key our pipeline uses for file names/output paths/
    report rows."""
    return (dataset, cluster_cfg.get("cell_annotation_key", cluster))
