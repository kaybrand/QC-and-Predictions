"""Reads the local filtered-barcode TSV (cluster_cfg["qc_guide"] -- a path,
not a portal concept) to list the unique subsample/in-vitro-system values
contributing to a cluster. Purely a local-file helper: makes no claim about
how that data maps onto any IGVF Data Portal object (e.g. a future
Filtered Barcode Set table) -- callers that need a portal-facing value
built from this should do that mapping themselves.
"""

import csv
import gzip


def unique_subsamples(ctx):
    with gzip.open(ctx.cluster_cfg["qc_guide"], "rt") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return sorted({row["subsamples"] for row in reader})
