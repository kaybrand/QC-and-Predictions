---
name: qc-filter-pseudobulks
description: Filter IGVF multiome pseudobulks to QC-passing ATAC fragment files and RNA count matrices, then generate scE2G config files. Use when the user has custom QC guides ready and wants to produce filtered multiome datasets.
allowed-tools: Read, Write, Edit, Bash, Glob, Grep
---

# QC Filter Pseudobulks → Filtered Multiome Datasets

This skill converts pseudobulks to filtered multiome datasets (ATAC + RNA) using custom QC guides,
then generates scE2G configuration files. It is used after the user has investigated per-cell QC
plots and decided on a filtering strategy for each cluster.

## Gather Inputs

Ask the user (or read from context) for:
- **Dataset name** (e.g. `igvf7`) — used in output paths and config file names
- **Pseudobulks directory** — root dir containing `annotation-{cell_type}` or `annotation-{cell_type}-IGVF*` subdirs
- **Cluster list** — names and any special renaming (e.g. `hcasm_calcification` → `hcasm_calcification_d6`)
- **QC guide for each cluster** — full path to the `.tsv.gz` file (barcode / subsample / analysis_accession)

## Key Paths

```
WDIR=/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks
Scripts:   $WDIR/scripts/filter_atac_fragments.py
           $WDIR/scripts/filter_rna_counts.py
QC guides: $WDIR/plots/{dataset}/{cluster}/
Output:    $WDIR/multiome_data/{dataset}/{cluster}/
Logs:      $WDIR/multiome_data/{dataset}/logs/
Pseudobulks (default):   /oak/stanford/groups/akundaje/sbaek1/web/igvf/pseudobulk/{dataset}/pseudobulks/
           Use this as the default when the user does not specify a path. If the user provides a $SCRATCH path, use that instead.
GTF:       /oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/uniformly_processed/GENCODEv43/IGVFFI9573KOZR.gtf.gz
Chrom sizes (default): $WDIR/reference/IGVF.DACC.GRCh38.chrom.sizes.tsv
Chrom sizes (alt):     $WDIR/reference/GRCh38_EBV.no_alt.chrom.sizes.tsv
scE2G config tables: /oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/config/mine/tables/
scE2G config yamls:  /oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/config/mine/
Template ATAC script: /oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/filter_atac_igvf7.sh
Template RNA script:  /oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/filter_rna_igvf7.sh
```

Always use `IGVF.DACC.GRCh38.chrom.sizes.tsv` as the default chromosome sizes file.
Offer `GRCh38_EBV.no_alt.chrom.sizes.tsv` as an alternative only if the user asks.

## QC Guide Format

QC guides are gzipped TSV files with a header and three columns:
```
barcode\tsubsample\tanalysis_accession
AAACAGCCACATAACT_IGVFSM1016TCPU\tIGVFSM0405BZBU\tIGVFDS0436AMDC
```

An annotation table (`annotation_{dataset}_rna_barcodes_human.tsv`, columns: analysis_accession,
barcode, subsample, annotation) is a useful authoritative source of all barcodes when available.
However, it may not exist for every dataset — it may need to be sourced from another consortium
lab or reconstructed. Do not assume it is available; ask the user.

Verify that barcodes in a new QC guide overlap with those in the existing
`$WDIR/plots/{dataset}/{cluster}/filtered_barcodes_with_subsamples.tsv.gz`
before proceeding. If overlap is negligible or zero, the most likely cause is an
**ATAC/RNA onlist mismatch**: the barcodes may have been drawn from the ATAC 10x Genomics
onlist rather than the RNA onlist (or vice versa). To diagnose this, check which onlist
the barcodes belong to. For IGVF7 the onlists are:

```
RNA onlist:  /oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/rna_whitelist.Robj  # object name: rna_wlist
ATAC onlist: /oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/atac_whitelist.Robj # object name: atac_wlist
```

For other datasets, locate the equivalent onlist files and check barcode membership accordingly.

## Common QC Guide Strategies

| Situation | Guide name pattern | Source |
|---|---|---|
| Standard QC thresholds, no quirks | `filtered_barcodes_with_subsamples.tsv.gz` | Existing file |
| Remove a bad subsample from QC-filtered set | `QC_filtered_excluding_{SUBSAMPLE}.tsv.gz` | Filter existing filtered_barcodes |
| Use annotation table, remove bad subsample | `McGinnis_filtered_excluding_{SUBSAMPLE}.tsv.gz` | Filter annotation table |
| Isolate a specific subsample's passing cells | `QC_filtered_{SUBSAMPLE}.tsv.gz` | Subset existing filtered_barcodes |

**Splitting one pseudobulk into multiple output clusters** (e.g. d6 vs d10): use `CELL_TYPE`
(the directory name) for pseudobulk lookup and `OUT_NAME` for output naming. Both output clusters
point to the same pseudobulk directory but use different QC guides to select disjoint barcode sets.

**Merging multiple pseudobulks into one output cluster**: this is common when pseudobulks were
created at the subsample level and you want a single combined output. List all matching
`annotation-{cell_type}-IGVF*` directories as sources — `filter_atac_fragments.py` and
`filter_rna_counts.py` automatically concatenate all matching directories for a given `--cell-type`.
The QC guide should contain barcodes from all subsamples to be merged.

## Step 1 — Create Output Directories

```bash
mkdir -p $WDIR/multiome_data/{dataset}/{cluster}  # one per output cluster
mkdir -p $WDIR/multiome_data/{dataset}/logs
```

## Step 2 — Filter ATAC Fragments (Slurm array job)

Environment: `mamba activate e2g_jamboree` + `ml biology bedops`
Script: `$WDIR/scripts/filter_atac_fragments.py`

Key arguments:
- `--qc-guide`    path to the QC guide .tsv.gz
- `--pseudobulks` root pseudobulks directory
- `--cell-type`   annotation name used in directory lookup (string after `annotation-`)
- `--chrom-sizes` chromosome sizes file (default: `IGVF.DACC.GRCh38.chrom.sizes.tsv`)
- `--out`         output path: `$OUTDIR/{out_name}/filtered_ATAC_fragments_{out_name}.tsv.gz`

The script automatically falls back to `annotation-{cell_type}` (no `-IGVF*` suffix) when
pseudobulks were not created at the subsample level.

Copy and adapt the template: `/oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/filter_atac_igvf7.sh`

## Step 3 — Filter RNA Count Matrices (Slurm array job)

Environment: `mamba activate filter_multiome`
Script: `$WDIR/scripts/filter_rna_counts.py`

Key arguments:
- `--qc-guide`    path to the QC guide .tsv.gz
- `--pseudobulks` root pseudobulks directory
- `--cell-type`   annotation name used in directory lookup
- `--out`         output path: `$OUTDIR/{out_name}/rna_count_matrix_{dataset}_{out_name}.mtx`
- `--gtf`         GTF file for Ensembl ID → gene symbol conversion
- `--standard-chromosomes-only`  always use this flag
- `--log`         `$OUTDIR/logs/{out_name}_gtf_mapping.txt`

Always use gene symbols (default). Do NOT use `--ensemblIDs-as-genes`.
Per-cluster run logs go to `$OUTDIR/logs/{out_name}_rna_run.log` via `2>> $RUN_LOG`.

Memory: 32G per cluster (each array task is an independent job).
Copy and adapt the template: `/oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/filter_rna_igvf7.sh`

## Step 4 — Generate scE2G Config Table

Before creating either config file, check whether `{dataset}_config.tsv` and `{dataset}_config.yaml` already exist. If the dataset config already exists, **add a new row** to the existing table rather than recreating it, and leave the YAML unchanged.

File: `/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/config/mine/tables/{dataset}_config.tsv`

Columns (tab-separated, leave HiC_file / HiC_type / HiC_resolution / alt_TSS / alt_genes empty):
```
cluster\trna_matrix_file\tatac_frag_file\tHiC_file\tHiC_type\tHiC_resolution\talt_TSS\talt_genes\tmodel_dir
{cluster}\t{OUTDIR}/{cluster}/rna_count_matrix_{dataset}_{cluster}\t{OUTDIR}/{cluster}/filtered_ATAC_fragments_{cluster}.tsv.gz\t\t\t\t\t\tmodels/multiome_powerlaw_v3,models/scATAC_powerlaw_v3
```

Note: `rna_matrix_file` is the MTX **directory** path (no trailing slash, no extension).
Reference example: `/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/config/mine/tables/igvf10_config.tsv`

## Step 5 — Generate scE2G Config YAML

File: `/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/config/mine/{dataset}_config.yaml`

Copy `igvf10_config.yaml` and replace all instances of `igvf10` with `{dataset}`.
All other fields (model paths, reference files, workflow params) remain unchanged.

## Checklist Before Submitting Jobs

- [ ] QC guide barcodes verified against existing `filtered_barcodes_with_subsamples.tsv.gz` (non-zero overlap); if zero, check ATAC/RNA onlist membership
- [ ] Output directories created on Oak
- [ ] Slurm log directory for job stdout/stderr exists on `$SCRATCH`
- [ ] ATAC script uses `IGVF.DACC.GRCh38.chrom.sizes.tsv` (default)
- [ ] ATAC script uses correct `CELL_TYPE` (directory lookup) vs `OUT_NAME` (output name)
- [ ] RNA script includes `--standard-chromosomes-only` and `--log`
- [ ] RNA script memory set to 32G per task
- [ ] Config table `rna_matrix_file` points to directory (no extension)
- [ ] Config yaml has correct dataset name in `results_dir` and `IGV_dir`
