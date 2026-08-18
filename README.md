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
   `filter_rna_counts.py`) — a cluster's `pseudobulk_annotation` can be a
   single raw pseudobulk directory name or a comma-separated list (to merge
   several raw annotations into one cluster; each is filtered independently
   and QC-passing rows are concatenated)
2. scE2G config generation and cluster metadata (`workflow/scripts/write_scE2G_config.py`)
3. Running scE2G (imported as a Snakemake module from the external `scE2G` repo)
4. Reformatting predictions/lists into the standard sharing format
5. Isolating candidate E2G pairs
6. Generating scE2G feature tables
7. Generating/updating Synapse manifests for all four shared result spaces
   and (optionally) uploading

Quality gating happens once, at Snakemake parse time
(`workflow/scripts/resolve_exclusions.py`): a cluster's cell/fragment/UMI
counts are checked against `exclusion.auto_thresholds`, and the outcome
(pass, below a specific threshold, or missing QC-guide/per-cell-QC inputs
entirely) is persisted per-dataset to `{output_dir}/cluster_stats/`. A
cluster that fails quality gets no rule instantiated for it at all — not a
failed job.

It's designed to run identically whether you target one cluster or a
few thousand — every rule is parameterized per-cluster, so scaling up is
just adding rows to a config, not changing the pipeline.

### `data_dir` vs `output_dir`

`data_dir` is read-only: pre-existing QC-guide plots and per-cell-QC
datatables the pipeline reads from (`{data_dir}/plots/...`,
`{data_dir}/datatables/...`), never writes to. Everything the pipeline
*writes* — filtered ATAC/RNA outputs, scE2G's own predictions/candidates/
features/QC-plot tree, cluster-stats tables, Synapse manifests, the
coverage report — is consolidated under a separate `output_dir` config key
(default `./results`, resolved against the repo root regardless of current
working directory; see `workflow/scripts/pipeline_paths.py`). Setting these
to different roots is deliberate: the code checkout, the QC-guide input
location, and the pipeline's own output location are three independent
things that don't have to (and often don't) live in the same place.

A second, separate destination — the IGVF Data Portal — is being built out
on this branch (`igvf-portal-submission`), in parallel with the
Synapse-only `synapse-submission` branch. There is no `main` branch — this
repo's branches (`igvf-portal-submission`, `synapse-submission`,
`CATlas-predictions`) act more like per-use-case forks than feature
branches awaiting merge into a trunk; `igvf-portal-submission` is the most
commonly used. (Other, temporary feature branches, e.g.
`download-from-igvf-portal`, still get merged and deleted normally — the
fork-like structure applies to these three, not every branch in the repo.)
It runs alongside the Synapse upload system, not
instead of it: the two destinations take different subsets of pipeline
outputs, packaged as different metadata objects, so they're developed
independently and only share the underlying cluster/exclusion resolution
both consume. See [IGVF Data Portal uploads](#igvf-data-portal-uploads) below.

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
    cell_annotations.py      # annotation_lookup_key() -- shared (dataset, cluster) key resolution for the
                              # live IGVF Portal CellAnnotation cache (state.db); never reads any preview TSV
    generate_report.py       # writes {output_dir}/report.tsv -- one coverage row per cluster (quality/
                              # annotation/prediction/reformat status), safe to re-run at any point
    list_synapse_orphans.py  # read-only: diffs Synapse's actual folder contents against the in-scope
                              # cluster set from report.tsv; flags orphans for manual review, never deletes
  envs/                 # conda environment definitions
  config/
    example_pipeline_config.yaml    # template — copy and fill in per dataset
    example_cluster_metadata.tsv    # template — copy and fill in per dataset
    {dataset}_pipeline_config.yaml  # real configs, gitignored (real paths/synIDs, data_dir, output_dir)
    {dataset}_cluster_metadata.tsv  # real cluster metadata, gitignored
local_scripts/           # local one-off tools, gitignored (not re-included by any .gitignore pattern)
  generate_pipeline_configs.py  # bulk-generates {dataset}_pipeline_config.yaml for datasets without a
                                 # hand-curated one, from the raw pseudobulk directory listing
resources/                # cluster-account-specific operational files, gitignored
  igvf_metadata_state.db  # the state.db this branch's igvf.state_db_path config key points at
  *.sbatch                # sbatch driver scripts for running this pipeline on Sherlock
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

`igvf-portal-submission` branch only (not present on the Synapse-only
`synapse-submission` branch), and not yet wired into the Snakemake pipeline
above (that's the next step; for now it's run standalone). Ten metadata
tables are registered so far (Prediction
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

### Reformat eligibility is gated on the live CellAnnotation cache, not a preview TSV

`reformat_predictions`/`reformat_predictions_thresholded`/`reformat_lists`
(the rules that turn scE2G's own output into portal-format Prediction
Tabular Files) only run for a `(dataset, cluster)` with a real, cached
CellAnnotation row already in `state.db` — computed at Snakemake parse time
as `common.smk`'s `REFORMAT_ELIGIBLE_CLUSTERS`, keyed through
`cell_annotations.py`'s `annotation_lookup_key()` (which resolves an
ATAC-only variant cluster's real CellAnnotation under its base, non-suffixed
name). **A quality-passing cluster with no cached CellAnnotation yet still
gets full scE2G predictions/candidates/features** — it's simply not
reformatted into portal format until `manage_igvf_metadata.py` has warmed
the cache for it. Core file generation is never blocked by Portal metadata
availability, and this is intentional: `state.db`'s cache is populated only
by an explicit, separate `manage_igvf_metadata.py` invocation, never
triggered automatically by a Snakemake run.

### Portal-facing packaging generated directly by this pipeline

A few files exist purely so `manage_igvf_metadata.py`'s uploader tables have
something to find on disk — none of these require Portal contact to
produce:

- `Neighborhoods/EnhancerList.bed.gz` + `.tbi` (`rules/enhancer_list_packaging.smk`) —
  bgzip/tabix of scE2G's own (already coordinate-sorted) `EnhancerList.bed`.
- Thresholded prediction `.bedpe.gz` + `.tbi` (`rules/qc_stats.smk`'s
  `bgzip_index_bedpe`) — Synapse submission uses scE2G's raw `.bedpe`
  as-is; the Portal's BEDPE Index File needs the bgzipped/indexed form.
- `rna_count_matrix_{dataset}_{cluster}.tar.gz` (`rules/filter_pseudobulks.smk`'s
  `package_rna_count_matrix`) — a flat, decompressed tar archive matching
  Filtered Matrix Files' own file-format spec, distinct from the
  gzipped-per-file directory scE2G itself reads as input.

All three are scoped to quality-passing clusters (`UPLOAD_ELIGIBLE_CLUSTERS`),
independent of CellAnnotation availability, same as core prediction
generation.
