# QC-and-Predictions

End-to-end Snakemake pipeline for the IGVF E2G Pillar Project: takes one or
more QC-filtered clusters from pseudobulk to Synapse-hosted scE2G data
products.

## Where this fits

QC threshold decision-making (per-cell QC plotting, choosing barcode
thresholds, producing the `filtered_barcodes_with_subsamples.tsv.gz` QC
guide per cluster under `plots/{dataset}/{cluster}/`) is a manual,
judgment-driven step and is **not** part of this pipeline — those scripts
remain in `scripts/`, untouched and untracked here.

This pipeline picks up once a QC guide exists for a cluster and takes it all
the way through:

1. ATAC fragment + RNA count matrix filtering (`workflow/scripts/filter_atac_fragments.py`,
   `filter_rna_counts.py`)
2. scE2G config generation and cluster metadata (`workflow/scripts/write_scE2G_config.py`)
3. Running scE2G (imported as a Snakemake module from the external `scE2G` repo)
4. Reformatting predictions/lists into the standard sharing format
5. Isolating candidate E2G pairs
6. Generating scE2G feature tables
7. Generating/updating Synapse manifests for all four shared result spaces
   and (optionally) uploading

It's designed to run identically whether you target one cluster or a
few thousand — every rule is parameterized per-cluster, so scaling up is
just adding rows to a config, not changing the pipeline.

## Layout

```
Snakefile              # legacy bulk ATAC/RNA filtering pipeline (kept working, superseded by workflow/ for new runs)
workflow/
  Snakefile             # main entry point for the end-to-end pipeline
  rules/                # Snakemake rule files
  scripts/              # Python/R scripts called by rules
  envs/                 # conda environment definitions
  config/
    example_pipeline_config.yaml    # template — copy and fill in per dataset
    example_cluster_metadata.tsv    # template — copy and fill in per dataset
    {dataset}_pipeline_config.yaml  # real configs, gitignored (real paths/synIDs)
    {dataset}_cluster_metadata.tsv  # real cluster metadata, gitignored
```

## Running

```bash
mamba activate run_snakemake9
snakemake -s workflow/Snakefile \
  --configfile workflow/config/{dataset}_pipeline_config.yaml \
  --executor slurm --profile slurm.smk9 --use-conda -p
```

Always dry-run first (`-n`) and inspect the plan, especially before the
first run for a new dataset. Synapse uploads and stale-entity deletions
default to dry-run/disabled (`synapse.dry_run: true`,
`synapse.confirm_delete: false`, `synapse.confirm_overwrite: false` in the
config) — flip these deliberately, one at a time, after reviewing what the
pipeline reports it would do.

See `workflow/config/example_pipeline_config.yaml` for the full config
schema and inline documentation of every field.
