"""Reads the local filtered-barcode TSV (cluster_cfg["qc_guide"] -- a path,
not a portal concept) to list the unique subsample/in-vitro-system values
contributing to a cluster. Purely a local-file helper: makes no claim about
how that data maps onto any IGVF Data Portal object (e.g. a future
Filtered Barcode Set table) -- callers that need a portal-facing value
built from this should do that mapping themselves.

Column is "subsample" (singular) -- confirmed against real QC guide files
(e.g. igvf2's), which also have a third column, "analysis_accession", not
yet used anywhere in this package but worth checking against the still-open
primary_pseudobulk_*_aliases/annotation_table_alias stubs in refs.py.
"""

import csv
import gzip


def unique_subsamples(ctx):
    with gzip.open(ctx.cluster_cfg["qc_guide"], "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return sorted({row["subsample"] for row in reader})
