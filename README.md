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

A second, separate destination — the IGVF Data Portal — is being built out
on the `igvf-portal-submission` branch/worktree (not yet merged to `main`).
It runs alongside the Synapse upload system, not instead of it: the two
destinations take different subsets of pipeline outputs, packaged as
different metadata objects, so they're developed independently and only
share the underlying cluster/exclusion resolution both consume. See
[IGVF Data Portal uploads](#igvf-data-portal-uploads) below.

## Layout

```
Snakefile              # legacy bulk ATAC/RNA filtering pipeline (kept working, superseded by workflow/ for new runs)
workflow/
  Snakefile             # main entry point for the end-to-end pipeline
  rules/                # Snakemake rule files
  scripts/              # Python/R scripts called by rules
    igvf_metadata/       # IGVF Data Portal metadata uploader (igvf-portal-submission branch only)
      registry.py          # declarative table/variant registration
      state.py             # SQLite upload ledger
      orchestrator.py       # plan-then-execute upload loop (preview/validate/upload modes)
      portal_client.py      # read-only portal checks + igvf_utils' iu_register.py wrapper
      context.py, refs.py, subsamples.py   # shared per-scope context, cross-table alias refs, helpers
      tables/               # one module per IGVF metadata table (Prediction Set, Signal Files, ...)
    manage_igvf_metadata.py  # CLI entrypoint for the above (igvf-portal-submission branch only)
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

## IGVF Data Portal uploads

`igvf-portal-submission` branch only — not yet merged to `main`, and not
yet wired into the Snakemake pipeline above (that's the next step; for now
it's run standalone). Ten metadata tables are registered so far (Prediction
Set, Prediction Tabular Files, Principal Pseudobulk Set, Filtered Barcode
Lists, Filtered ATAC Fragment Files, Filtered Matrix Files, Signal Files,
ATAC Index File, BEDPE Index File, Documents) — several still block on
pieces marked TBD in the code (a couple of alias formulas and one live
portal lookup), which read loudly rather than silently uploading
incomplete metadata.

Real writes always go through `igvf_utils`' own `iu_register.py`, never a
hand-rolled API call, and default to a dry, reviewable mode:

```bash
python workflow/scripts/manage_igvf_metadata.py \
  --pipeline-config workflow/config/{dataset}_pipeline_config.yaml \
  --cluster-keys "{dataset}/{cluster},..." \
  --state-db path/to/igvf_metadata_state.db \
  --manifest-dir path/to/igvf_metadata_manifests \
  --mode preview
```

- `--mode preview` (default): no contact with `iu_register.py` at all —
  just writes the post/patch TSVs for every eligible row to
  `--manifest-dir`, for review. The only network contact is a read-only
  alias-existence check against the portal.
- `--mode validate`: also runs `iu_register.py --dry-run` against each
  generated TSV — exercises real schema validation/type-casting, still
  zero portal writes.
- `--mode upload`: runs `iu_register.py` for real. Only reachable by
  explicitly choosing this mode.

Upload state (what's uploaded/pending/failed per cluster+table+row) lives
in its own SQLite ledger, entirely separate from the Synapse manifest
system — the two destinations never share state, only the
`--cluster-keys` eligibility list computed by the same
`resolve_exclusions.py` both consume.
