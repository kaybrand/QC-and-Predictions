---
name: qc-filter-pseudobulks
description: Filter IGVF multiome pseudobulks to QC-passing ATAC fragment files and RNA count matrices, then generate scE2G config files. Use when the user has custom QC guides ready and wants to produce filtered multiome datasets.
allowed-tools: Read, Write, Edit, Bash, Glob, Grep
---

# QC Filter Pseudobulks → Filtered Multiome Datasets

Turns downloaded IGVF primary pseudobulks into QC-passing ATAC fragment files and RNA
count matrices, using per-cluster QC guides built from thresholds a human chose. Since
2026-10 this all happens inside the **QC-and-Predictions** repo (`igvf-portal-submission`
branch): QC guides come from one plotting script, and the filtering + scE2G config are
done by the pipeline's own rules. Do not hand-run the old per-dataset shell templates.

**Never use `/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks`**
(retiring) for a new dataset — not for scripts, guides, outputs, `data_dir`, or reference
files. Never write to the shared scE2G checkout (`.../Users/kaybrand/scE2G_preprint/scE2G`).

## Key Paths

```
REPO=/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks-igvf-portal-submission
Pseudobulks root:  /scratch/users/kaybrand/IGVF_Data_Portal_download_08-2026
                   (= --download-root of workflow/scripts/download_portal_pseudobulks.py;
                   the name is historical -- the 2026-10 syncs land here too; do not rename,
                   state.db's portal_files.local_path stores absolute paths into it)
  layout:          {root}/{dataset}/pseudobulks/{annotation}-{IGVFSM...}/          (IGVFDS* datasets)
                   {root}/{dataset}/pseudobulks/annotation-{annotation}-{IGVFSM...}/ (igvfN datasets)
  per set:         fragments.tsv.gz, rna_counts_mtx.h5ad, per_cell_qc.tsv.gz
Per-cell QC tables: $REPO/results/QC_datatables/{dataset}_data/{cluster}_per_cell_qc.tsv
QC plots + guides:  $REPO/results/QC_plots/{dataset}/{cluster}/
Plotting script:    $REPO/local_scripts/plot_per_cell_qc.R          (conda env: qc_per_cell)
Table builder:      $REPO/workflow/scripts/build_qc_datatables.py   (env: igvf_utils_env)
Filter scripts:     $REPO/workflow/scripts/filter_atac_fragments.py, filter_rna_counts.py
                    (run by workflow/rules/filter_pseudobulks.smk; env workflow/envs/filter_multiome_env.yaml)
Filtered outputs:   $REPO/results/multiome_data/{dataset}/{cluster}/atac_fragments_{dataset}_{cluster}.tsv.gz(.tbi)
                    $REPO/results/multiome_data/{dataset}/{cluster}/rna_count_matrix_{dataset}_{cluster}/
Chrom sizes:        $REPO/reference/IGVF.DACC.GRCh38.chrom.sizes.tsv
GTF:                /oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/uniformly_processed/GENCODEv43/IGVFFI9573KOZR.gtf.gz
Pipeline configs:   $REPO/workflow/config/{dataset}_pipeline_config.yaml  (gitignored)
scE2G cell_clusters: $REPO/workflow/config/{dataset}_cell_clusters.tsv   (written by the pipeline)
Envs:               /oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/{qc_per_cell,igvf_utils_env,filter_multiome}
```

`$REPO/results/plots/` holds the human QC guides of igvf9/igvf13 only — never write
placeholder output there. For IGVFDS datasets the human guides live in `results/QC_plots/`.

## Step 1 — Per-cell QC tables (one per cluster)

```bash
cd $REPO
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python \
  workflow/scripts/build_qc_datatables.py \
  --pseudobulks-root /scratch/users/kaybrand/IGVF_Data_Portal_download_08-2026 \
  --out-dir $REPO/results --datasets IGVFDS... [--force]
```

Concatenates every subsample's `per_cell_qc.tsv.gz` (header once, subsample dirs in sorted
order) and normalizes the new Portal column names (`analysis_set_accession` →
`analysis_accession`, `barcode_sample` → `barcode`). Run it as a job, not on a login node.

**Merged cluster** (several annotations → one cluster, e.g. IGVF7's
`hcasm_calcification` = `hcasm_calcification_d10` + `hcasm_calcification_d6`): build one
table spanning both, keyed by the merged cluster name:

```bash
cd $REPO/workflow/scripts && /oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python -c "
import build_qc_datatables as b
print(b.build_one('/scratch/users/kaybrand/IGVF_Data_Portal_download_08-2026', 'IGVFDS3445MTOM',
      'hcasm_calcification', ['hcasm_calcification_d10', 'hcasm_calcification_d6'],
      '$REPO/results/QC_datatables/IGVFDS3445MTOM_data/hcasm_calcification_per_cell_qc.tsv'))"
```

Check the row count equals the sum of the constituents and that no barcode is duplicated.

## Step 2 — QC plots and QC guides (human thresholds)

The human reviews `ATAC_QC_plots.pdf` / `RNA_QC_plots.pdf` and supplies thresholds,
usually as a TSV with columns `Dataset, Cluster, RNA UMIs (min), Gene Count (min),
% mito (max), % ribosomal (max), ATAC fragments (min), nucleosomal signal (max),
TSS enrichment (min)`. Map them onto the script's flags (blank = script default):

| Threshold column | flag | default |
|---|---|---|
| RNA UMIs (min) | `--rna-min` | 1000 |
| Gene Count (min) | `--gene-min` | 0 |
| % mito (max) | `--pct-mt-max` | 30 |
| % ribosomal (max) | `--pct-ribo-max` | 100 (no cutoff) |
| ATAC fragments (min) | `--atac-min` | 1000 |
| nucleosomal signal (max) | `--nuc-max` | 1.5 |
| TSS enrichment (min) | `--tss-min` | 3 |

Also available: `--rna-max --gene-max --atac-max --pct-dup-max --frip-min`.

```bash
source /oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/etc/profile.d/conda.sh
source /oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/etc/profile.d/mamba.sh
mamba activate qc_per_cell
cd $REPO/results
Rscript $REPO/local_scripts/plot_per_cell_qc.R \
  QC_datatables/{dataset}_data/{cluster}_per_cell_qc.tsv QC_plots/{dataset}/{cluster} \
  --rna-min 1000 --gene-min 1000 --pct-mt-max 30 --atac-min 10000 --nuc-max 1 --tss-min 3
```

Writes, per cluster: the plots, `filtered_barcodes_with_subsamples.tsv.gz` (**the QC
guide**), `filtered_cell_subsample_metrics.tsv` (read by the pipeline's quality gate) and
`qc_thresholds.tsv` (attached to the IGVF QC_documents manifest). Keep all four together —
the pipeline derives the other paths from the guide's directory. Save the exact commands
(e.g. a small script under `$SCRATCH/temp`) and verify each `qc_thresholds.tsv` matches
what the human gave. For a merged cluster with the same thresholds as its constituents,
cells kept should equal the constituents' sum.

Exploratory plots made before thresholds exist use placeholders (`--rna-min 1000
--gene-min 1000 --pct-mt-max 30 --atac-min 1000 --nuc-max 1.5 --tss-min 3`); their
guides are **not** human-approved — never point a config at them.

## QC guide format

Gzipped TSV, header + three columns:
```
barcode\tsubsample\tanalysis_accession
AAACAGCCAGGAACAT_IGVFSM1016TCPU\tIGVFSM0405BZBU\tIGVFDS0436AMDC
```
`barcode` is the full `{16-mer}_{IGVFSM...}` string — the same string as column 4 of
`fragments.tsv.gz` and the h5ad `obs_names`, old and new datasets alike. Its IGVFSM suffix
can differ from `subsample` (multiplexed lanes); matching is always on the full string,
never the 16-mer alone.

If a guide's barcodes barely overlap the pseudobulk's, suspect an ATAC/RNA onlist mix-up
(IGVF7 onlists: `/oak/stanford/groups/engreitz/Users/kaybrand/IGVF_Consortium/IGVF7/{rna,atac}_whitelist.Robj`,
objects `rna_wlist` / `atac_wlist`).

## Step 3 — Pipeline config (one per dataset)

```bash
cd $REPO; R=$REPO/results
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python \
  local_scripts/generate_pipeline_configs.py \
  --datasets IGVFDS... --qc-guide-dir $R/QC_plots \
  --data-dir $R --output-dir $R --qc-datatables-dir $R/QC_datatables [--force]
```

- The dataset key must be the IGVFDS accession (download dirs and Portal aliases use it).
- Each cluster: `pseudobulk_annotation` (the dir annotation; becomes `--cell-type`),
  `qc_guide`, `models: [multiome_powerlaw_v3, scATAC_powerlaw_v3]`.
- **Merged cluster** (igvf18 `mcf7` style): `pseudobulk_annotation: a,b`
  (comma-separated) and `submitter_comment: Merged a and b pseudobulks`. Both filter
  scripts split on the comma, find each annotation's dirs, filter them all against the
  one merged guide and concatenate (verified for IGVF7: 9 dirs, all 5,588 guide barcodes
  found exactly once). The generator's `MERGED_CLUSTERS` drops the constituents as
  separate clusters — keep them only if the human wants both.
- `--force` rewrites only the datasets named in `--datasets`.

## Step 4 — Cell Annotation (the Portal is the source of truth)

CellAnnotation / SampleTermName / SampleTermID come from the IGVF Portal through the 24h
cache in `$REPO/resources/igvf_metadata_state.db` (refreshed by the downloader's discovery
GET or the driver's warm stage) — never from `lab_annotations_with_cl.tsv` or
`cell_annotations_by_dataset_cluster.tsv` (snapshots; do not edit them). A cluster missing
any of the three on the Portal is an error, not a blank.

A merged cluster otherwise gets the bare cell-type term. To give it a specific annotation,
or to correct a Portal typo in our own records, set a durable override — its own table,
never written by a cache refresh, re-applied on every derive:

```bash
/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python \
  $REPO/workflow/scripts/set_cell_annotation_override.py \
  --state-db $REPO/resources/igvf_metadata_state.db --dataset IGVFDS3445MTOM --cluster hcasm_calcification \
  --cell-annotation "calcification HCASMC 2105 hTERT dCas9 smooth muscle cell of the coronary artery derived from HCASMC-hTERT" \
  --cell-qualifier "calcification HCASMC 2105 hTERT dCas9" --reason "merged d6+d10; no day"
# --list shows all overrides; --clear removes one
```

## Step 5 — Filtering and scE2G config: the pipeline does it

```bash
cd $REPO
sbatch resources/run_pipeline.sbatch {dataset} default -n                          # ALWAYS dry-run first
sbatch resources/run_pipeline.sbatch {dataset} default --snakemake-arg=--retries=1
```

`atac_fragment_file` / `rna_count_matrix` run the filter scripts with the chrom sizes,
GTF, `--standard-chromosomes-only` and log flags built in; `write_scE2G_config.py` writes
`workflow/config/{dataset}_cell_clusters.tsv` at parse time. Do not hand-write scE2G config
tables or YAMLs. Standalone use of a filter script (debugging only) must point `--out` at
`$SCRATCH/temp`, never at the results tree.

## Checklist

- [ ] Per-cell QC tables built from the CURRENT download (re-issued Portal files replace old ones)
- [ ] Human thresholds mapped onto flags; blanks left at defaults; commands saved
- [ ] Each cluster's `qc_thresholds.tsv` matches the human's values; guide rows = cells kept
- [ ] Merged cluster: table rows = sum of constituents, 0 duplicate barcodes, comma-separated
      `pseudobulk_annotation`, `submitter_comment`, annotation override if the human wants one
- [ ] Config: IGVFDS key, `qc_guide` under `results/QC_plots/`, absolute `data_dir`/`output_dir`,
      no `QC_pseudobulks` anywhere
- [ ] Dry-run the driver before the real run
