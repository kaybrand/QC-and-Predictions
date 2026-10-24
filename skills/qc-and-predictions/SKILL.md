---
name: qc-and-predictions
description: Run the QC-and-Predictions pipeline end to end on a new IGVF primary-pseudobulk dataset (an IGVFDS accession) — Portal download, per-cell QC plots, human thresholds, pipeline config, scE2G predictions, and IGVF Portal manifests ready for a human to validate and upload. Use when asked to process, sync, QC, run predictions for, or prepare Portal manifests for a new IGVFDS dataset in the QC_pseudobulks-igvf-portal-submission repo, or to resume such work.
allowed-tools: Read, Write, Edit, Bash, Glob, Grep
---

# QC-and-Predictions: new IGVFDS dataset, Portal to Portal

Proven end to end on IGVFDS3445MTOM (= igvf7, 2026-10-19 → 10-22): primary pseudobulks
downloaded from the IGVF Portal → filtered multiome data → scE2G predictions → 84
manifest rows → uploaded by the user. Every command below was run for real.

Works in the QC-and-Predictions repo (`$REPO` below, branch `igvf-portal-submission`).
**Read first:** `$REPO/CLAUDE.md` (hard rules) and the newest `$REPO/.baton/` file
(`ls -t $REPO/.baton/ | head`). Batons hold the live status; this skill holds the method.

## Hard rules (from CLAUDE.md and the user — non-negotiable)

- **Never upload.** No `iu_register.py`, no `--mode upload`, from anything you run.
  Manifests stop at `--mode preview`; `--mode validate` only with the user's explicit OK.
  The user types the upload themselves.
- **Never expose IGVF credentials** (`IGVF_API_KEY` / `IGVF_SECRET_KEY`). Never print,
  echo, `env`/`printenv`, or put them in a file, log, commit or command line. Never read
  `~/.bash_profile`, `~/.bashrc`, `~/.netrc` or keypair files. If they leak, say so loudly
  at once — the repo is public. Commands that even *mention* the variable names get blocked;
  check their presence only inside an sbatch script.
- **Writable paths:** this repo and `$SCRATCH/temp`. Never edit the shared scE2G checkout
  (`/oak/.../Users/kaybrand/scE2G_preprint/scE2G`) — Snakemake writing envs into its
  `.snakemake/conda` is the only exception, and only through the pipeline.
- **Commits:** modular, cherry-pickable, succinct, and **no Co-Authored-By trailer**. Commit
  only when asked; never push.
- **Never touch datasets that are live on the Portal:** igvf0–6, 9–13, 17, 18, kasowski,
  and any IGVFDS the user has uploaded. igvf18 and igvf6 manifests were hand-corrected:
  never regenerate them.
- **Sherlock:** no compute or network Python on a login node — `sbatch` (or an `sh_dev`
  allocation). Every job gets an explicit `--time`. Wait ≥60 s between status checks;
  prefer `--dependency` chains to polling. Use env vars (`$SCRATCH`, `$OAK`) in paths you
  write.
- **Scope:** multiome only. A dataset whose sets lack ATAC fragments (RNA-only) is out of
  scope.

## Session hygiene

- `git fetch && git log --oneline origin/igvf-portal-submission..HEAD` first: an unfetched
  remote once made a "complete and verified" output wrong (missing table).
- Other sessions (and humans) share this working tree. Re-check `git status`,
  `squeue -u $USER` and fresh globs right before any real action.
- Partitions: `dev` rejects batch jobs. Short diagnostics: `-p normal,engreitz,owners`.
  Keep multi-hour drivers on `engreitz` (`owners` can preempt them).
- Fresh clone: `resources/` is gitignored. `mkdir -p resources && cp
  workflow/sbatch/run_pipeline.sbatch.example resources/run_pipeline.sbatch`, then set `REPO`/`FORGE`.

## Fixed paths and tools

These are Kayla Brand's installs and download root (2026-10). Another lab member
substitutes their own `PY`/miniforge, `ROOT` and scE2G checkout, keeping the layout.

```
REPO=/oak/stanford/groups/engreitz/Projects/IGVF-E2GPillarProject/QC_pseudobulks-igvf-portal-submission
PY=/oak/stanford/groups/engreitz/Users/kaybrand/software/miniforge3/envs/igvf_utils_env/bin/python
   # every repo Python script; has SQLite 3.50 (login-node python 3.7 breaks state.py)
ROOT=/scratch/users/kaybrand/IGVF_Data_Portal_download_08-2026
   # pseudobulks root; "08-2026" is historical; never rename (state.db stores absolute paths)
STATE=$REPO/resources/igvf_metadata_state.db      # WAL SQLite on Lustre
LOCK=$STATE.warmlock                              # every writer holds: flock $LOCK
CONDA_PREFIX_SHARED=/oak/stanford/groups/engreitz/Users/kaybrand/scE2G_preprint/scE2G/.snakemake/conda
```

- Never use `/oak/.../IGVF-E2GPillarProject/QC_pseudobulks` (retiring) for anything.
- Read `state.db` read-only with `sqlite3.connect("file:...?mode=ro", uri=True)` under `$PY`.
  Only one writer at a time, across hosts (WAL on Lustre).
- The pipeline driver is `workflow/scripts/run_pipeline.py`, always via
  `resources/run_pipeline.sbatch <dataset> <default|local_only> [driver args]`. It bakes in
  `--use-conda --conda-prefix $CONDA_PREFIX_SHARED --rerun-triggers mtime`. If you ever
  call `snakemake` directly you MUST pass both of those flags (omitting them caused a
  multi-hour env-rebuild race), dry-run first, use absolute target paths, and never put a
  positional target after `--omit-from` (it swallows it — a 5.5 h incident).

## The state.db tables you will meet

| table | holds | written by |
|---|---|---|
| `portal_files` | the download ledger: every Portal file, md5, `local_path`, `download_state`; also the source of IGVFDS **file aliases** | downloader |
| `cell_metadata_primary_pseudobulks` / `_principal_pseudobulks` | raw Portal cache of PseudobulkSets (Cell Annotation, term, qualifier, aliases); 24 h TTL | `cell_metadata.fetch_if_stale`, downloader discovery |
| `cell_annotations` | per (dataset, cluster) derived annotation | `derive_scopes`, every default run |
| `cell_annotation_overrides` | hand corrections that survive every refresh | only `set_cell_annotation_override.py` |
| `uploads` | per manifest row: pending / uploaded / failed, payload hash, Portal id | `manage_igvf_metadata.py` |

The 24 h cache exists so a run reuses one fetch instead of hammering the Portal; never
remove the save step from `fetch_if_stale`. Cell Annotation truth is the Portal cache —
never `lab_annotations_with_cl.tsv` or `cell_annotations_by_dataset_cluster.tsv`. A cluster
missing CellAnnotation, SampleTermName or SampleTermID is an error to report.

## Step 1 — Download from the Portal (job)

```bash
flock $LOCK $PY workflow/scripts/download_portal_pseudobulks.py \
    --download-root $ROOT --state-db $STATE --datasets <IGVFDS…> --dry-run   # discovery only
# then the same without --dry-run
```

- Defaults are the minimum multiome set: `fragments.tsv.gz`, `rna_counts_mtx.h5ad`,
  `per_cell_qc.tsv.gz`, only for sets having both fragments and a cell-by-gene matrix
  (`--exclude-file-formats bigBed`, `--require-content-types "fragments,cell by gene matrix"`).
  The `fragments.bb` bigBed duplicates the tsv; skip it.
- Discovery also refreshes the Cell Annotation cache. It only looks at Kundaje-lab sets
  (`DEFAULT_LAB="/labs/anshul-kundaje/"` in `portal_files.py`).
- "Ready" = fragments + matrix, but QC plots also need the per-cell QC report: a dataset
  without it builds 0 QC tables (IGVFDS9947FPYX was "ready" with no reports at first).
- The downloader exits 1 if any file fails, so chain dependents with `afterany` plus a
  ledger check (`download_state`), not `afterok`.
- Layout: `$ROOT/<IGVFDS>/pseudobulks/{annotation}-{IGVFSM…}/`. New-style conventions:
  set alias `anshul-kundaje:pseudobulk-{IGVFDS}-{ann}-{SM}`, per-cell-QC alias suffix
  `_tsv_gz`, per-cell-QC columns `analysis_set_accession` / `barcode_sample`. Barcodes are
  `{16-mer}_{IGVFSM}`. The code handles all of this; do not hand-build aliases — they come
  from the ledger.

## Step 2 — Per-cell QC tables and exploratory plots (job)

```bash
$PY workflow/scripts/build_qc_datatables.py --pseudobulks-root $ROOT --out-dir results \
    --datasets <IGVFDS> --force
# -> results/QC_datatables/<IGVFDS>_data/<cluster>_per_cell_qc.tsv
```

Then in conda env `qc_per_cell`:
`Rscript local_scripts/plot_per_cell_qc.R <table> results/QC_plots/<IGVFDS>/<cluster> [thresholds]`
with placeholder thresholds for exploring:
`--rna-min 1000 --gene-min 1000 --pct-mt-max 30 --atac-min 1000 --nuc-max 1.5 --tss-min 3`
(always pass `--gene-min`: the script's default 0 is not comparable with older datasets;
its default `--frip-min 0` drops FRiP = 0 cells). Output per cluster: `RNA_QC_plots.pdf`,
`ATAC_QC_plots.pdf`, and the guide files.

- **Never** run `results/QC_datatables/run_per_cell_qc_plots.sh` (writes `results/plots/`,
  clobbering igvf9/igvf13 human guides).
- **Humans QC in parallel.** Never re-render a cluster with placeholders unless its current
  `qc_thresholds.tsv` already holds the placeholders (a human guide was overwritten this way
  on 2026-10-21).
- Missing QC values are empty fields, not "NA". ~0.04% of cells have RNA but zero fragments;
  the plot script excludes them.
- **Before making guides, read `references/qc-guides.md`:** threshold-to-flag mapping, guide
  format, merged-cluster tables, and the checks to run.

## Step 3 — Human thresholds → real QC guides

Details: `references/qc-guides.md`. Give the user a TSV sheet (one row per cluster, threshold columns blank) — e.g.
`results/QC_plots/IGVFDS_qc_thresholds.tsv`. When filled, re-run `plot_per_cell_qc.R` per
cluster with the human's values into `results/QC_plots/<IGVFDS>/<cluster>/` and verify every
`qc_thresholds.tsv` equals the sheet and guide rows = cells retained.

- Bash sheet readers: tab-IFS collapses empty fields (`tr '\t' '\037'` first); the last
  line may lack a newline; bash 4.2 + `set -u` rejects empty arrays (`${a[@]+"${a[@]}"}`).
- `plot_per_cell_qc.R` uses strict `>`/`<`; `explore_qc_thresholds.R` uses `>=`/`<=` —
  counts can differ.
- The user decides exclusions. The pipeline's quality gate (`exclusion.auto_thresholds`:
  ≥100 cells, ≥2,000,000 fragments, ≥1,000,000 UMIs; `resolve_exclusions.py`, at parse time)
  drops failing clusters automatically — a failing cluster gets no rules at all, not a
  failed job. Only Multiome files are submitted (`enabled_families`); scATAC outputs are
  made but never uploaded.
- Per cluster the guide dir holds `filtered_barcodes_with_subsamples.tsv.gz` (the guide),
  `filtered_cell_subsample_metrics.tsv`, `qc_thresholds.tsv`. Everything downstream finds
  them via `dirname(qc_guide)`.

## Step 4 — Pipeline config

```bash
R=$REPO/results
$PY local_scripts/generate_pipeline_configs.py --datasets <IGVFDS> --qc-guide-dir $R/QC_plots \
    --data-dir $R --output-dir $R --qc-datatables-dir $R/QC_datatables --force
grep -c QC_pseudobulks/ workflow/config/<IGVFDS>_pipeline_config.yaml   # must be 0
```

Shape (per cluster): `pseudobulk_annotation`, `qc_guide` (absolute), `models:
[multiome_powerlaw_v3, scATAC_powerlaw_v3]`. Absolute `data_dir`/`output_dir` = `$REPO/results`.
Merged cluster: `pseudobulk_annotation: a,b` plus a `submitter_comment` ("Merged a and b
pseudobulks"); add it to `MERGED_CLUSTERS` in the generator (`keep_constituents` to keep
the parts) and build its QC table with `build_one(root, ds, merged, [a, b], out)`.

## Step 5 — Cell Annotation overrides (only if needed)

```bash
$PY workflow/scripts/set_cell_annotation_override.py --state-db $STATE --dataset <IGVFDS> \
    --cluster <c> --cell-annotation "<text>" --cell-qualifier "<text>" --reason "<why, who, date>"
# --list shows overrides; --clear removes one. Hold $LOCK.
```

Used for Portal typos (IGVF7 "DEME" → "DMEM") and merged clusters (a qualifier with no day).
The raw cache keeps the Portal's text.

## Step 6 — Build conda envs ONCE, before any parallel run

Dry runs never build envs. If any `workflow/envs/*` file changed (or the env is new),
concurrent runs will all build the same env, fail, and can **corrupt conda's shared package
cache** (2026-10-22: a later "successful" rebuild had a 0-byte `datetime.py`). So:

```bash
P=$REPO/results/uniformly_processed
E=$(sbatch --parsable --time=03:00:00 --mem=16G --job-name=envbuild resources/run_pipeline.sbatch \
    <IGVFDS> local_only --snakemake-arg=--conda-create-envs-only \
    --snakemake-arg=$P/synapse_manifests/filtered_data_manifest.tsv)   # this target pulls in post_process.yaml
# also run it once with no target (rule all) to cover the other envs
```

Then chain the real runs with `--dependency=afterok:$E`. Verify an env by importing in it,
e.g. `<env>/bin/python -c "import numpy, anndata, h5py"`.

## Step 7 — Dry run, then the real run

```bash
sbatch --time=01:00:00 resources/run_pipeline.sbatch <IGVFDS> default -n
sbatch resources/run_pipeline.sbatch <IGVFDS> default --snakemake-arg=--retries=1
```

- Optional offline preflight: `sbatch resources/predict_fanout_scopes.sbatch <IGVFDS>` (a
  state.db copy, no network) reports `N/M scope(s) resolved` and why the rest don't.
- Stages: preflight → warm (fetch_if_stale + derive_scopes, under the lock) → snakemake
  (Slurm executor) → manifest (**preview** only) → audit. Exit **0** = `RESULT: complete and
  verified` (never weaken this); **1** = output produced but incomplete — read the itemised
  `PROBLEM:`/`DEGRADED:` lines (often one unresolved Cell Annotation); **2** = config/preflight
  error, nothing ran. A Portal outage still runs snakemake, skips manifest/audit, exits 1.
- Status files: `results/report.tsv` (per cluster, `manifest_status`/`manifest_gap_reason`),
  `results/igvf_metadata/<IGVFDS>_cell_annotation_status.tsv`, and
  `results/igvf_manifests/<IGVFDS>/manifest_coverage.tsv`.
- `default -n` needs a cache <24 h old (else SnapshotError) — refresh with
  `sbatch resources/refresh_cell_metadata.sbatch <IGVFDS…>` (or run a manifest preview, which
  fetches). Every Slurm sub-job re-checks the snapshot against the original fetch time, so a
  run longer than (24 h − cache age at start) gets scattered late SnapshotErrors (igvf9: 30
  of 41 gaps). Start right after a refresh; if it happens, refresh and resubmit — only the
  failed jobs redo.
- The dry run cannot evaluate scE2G checkpoints; the real run grows (IGVF7: 218 → 257 steps,
  6 clusters, 6 h 57 m).
- **Always `--retries=1`:** memory (`determine_mem_mb`, `workflow/rules/common.smk:64`) floors
  at 8 GB and doubles per attempt. Without it, big datasets die OOM in `atac_fragment_file`
  (`sort-bed` SIGKILL). Escalate with `--snakemake-arg=--set-resources=atac_fragment_file:mem_mb=…`.
- `local_only` mode: no Portal contact, no manifests — for predictions only. A single target:
  `--snakemake-arg=$P/<IGVFDS>/qc_plots/predictions_qc_report.html` (absolute; the driver
  appends it last, so it is a safe positional).
- Several datasets may run at once (state writes serialize on the lock) — **after** Step 6.
- Debug a failed rule: `grep -h "Error in rule" <driver .err> | sort | uniq -c`, then the per-job
  log `.snakemake/slurm_logs/rule_<rule>/<ds>_<cluster>/<jobid>.log` (Python tracebacks land
  there; the rule's own log may be empty) and `sacct -j <jobid> -o JobID,State,MaxRSS,ReqMem`.
- Completed outputs are kept; re-submitting resumes.
- **Re-runs of an existing dataset** (not the first run):
  - `--rerun-triggers mtime` ignores changed `shell:`/`params:`. Cell metadata reaches the
    reformat rules as params, so stage 1 deletes reformatted files whose
    `# CellAnnotation:` / `# SampleTermName:` / `# SampleTermID:` headers disagree with the
    cache (check alone: `$PY workflow/scripts/stale_reformats.py <IGVFDS>`).
  - Reformat only, no scE2G: `... default --sce2g-modules false --snakemake-arg=portal_reformat`.
    Manifest only: `--skip-snakemake`.
  - If the dry run plans `generate_atac_matrix → compute_kendall → arc_e2g` with
    `reason: Updated input files: .../uniformly_processed/<IGVFDS>`, the directory's mtime
    moved. Reset it (`touch -h -d` to the oldest `atac_matrix.rds` minus 60 s) and dry-run
    again; never `snakemake --touch`, never move `qc_plots/`.
  - Avoid `--forcerun` (rebuilds everything; its `nargs='+'` swallows a following target).

## Step 8 — Manifests: review, validate (with OK), upload (user)

Output: `results/igvf_manifests/<IGVFDS>/round<N>_<table>[_<variant>]_<post|patch>.tsv`
(filename order = upload order) and `manifest_coverage.tsv`. 14 objects per cluster (QC
document, principal pseudobulk set, barcode list, prediction set, filtered fragments + index,
RNA matrix, signal bigWig, elements bed + index, full / thresholded / bedpe predictions +
index).

Regenerate a preview on its own (no snakemake):
`sbatch --time=02:00:00 resources/upload_dataset.sbatch --mode preview <IGVFDS>`
(back up `results/igvf_manifests/<IGVFDS>/` to `$SCRATCH/temp` first).

The same thing by hand (inside a job, holding `$LOCK`):
`$PY workflow/scripts/manage_igvf_metadata.py --pipeline-config workflow/config/<IGVFDS>_pipeline_config.yaml --cluster-keys <IGVFDS> --state-db $STATE --manifest-dir results/igvf_manifests --mode preview`
(a bare dataset in `--cluster-keys` expands to all its clusters; `--tables <name>` limits it).

Reading the files:
- `_post` = not on the Portal at that alias; `_patch` = live with a changed payload hash.
  A table whose rows are all `unchanged` writes no file — check `manifest_coverage.tsv`
  before calling a missing file a gap.
- `unchanged` (ledger `uploaded` + same hash) skips the Portal entirely, so Portal-side
  deletions are invisible. The post/patch lookup uses the search index and can lag a fresh
  write. Confirm with a read-only `get_by_alias` before resetting any ledger row (resetting a
  live row emits a spurious PATCH).
- Preview never marks rows `uploaded`; only a real upload's read-back does.
- A PATCH never removes a field that is absent from the payload.

Check before handing over:
- every cluster present, all rows `planned-post` (or expected patches), 0 invalid;
- cell_type / cell_qualifier as expected (overrides applied); merged sets list every
  constituent's sets and files;
- compare against a backup by row content, not bytes (row order is not deterministic).

Then tell the user. **On their OK only:** `--mode validate` (= `iu_register.py --dry-run`;
md5s every file, ~10 GB/dataset). Validate cannot catch server-side enum or permission
rules — only a real POST does. Upload is the user's:
`resources/upload_dataset.sbatch --mode upload --until-done <IGVFDS>`.

How upload behaves (so you can read its logs): passes over tables in registry order; each
table plans → writes TSV → `iu_register.py` → reads every object back from the Portal
database → marks `uploaded`. Rows whose dependency isn't `uploaded` in the ledger wait for
the next pass (~6 passes per dataset); `--until-done` loops and stops on any failure or on
a pass with no progress. Preview/validate deliberately build and look up the WHOLE chain,
because they approve an `--until-done` upload. A `failed` row is re-planned on the next run.
An `iu_register.py` exit 0 does not prove every row landed; the per-row read-back does.
Synapse is a separate destination with its own ledger (`synapse.dry_run: true` by default);
never touch it unless asked.

## Portal behaviour to expect

- **Schema enums change.** 2026-08: prediction set `file_set_type` → `element-gene links`;
  2026-10: pseudobulk set `pseudobulk analysis` removed → `gene expression pseudobulk`
  (`workflow/scripts/igvf_metadata/tables/principal_pseudobulk_set.py:88`). Symptom: a 422 on
  upload round N, everything downstream stuck. Check `https://api.data.igvf.org/profiles/<type>.json`
  and its changelog. Changing a constant changes the payload hash of already-uploaded rows,
  so other datasets would plan PATCHes — warn the user; never run live datasets to "fix" it.
- **Load.** `/multireport/` fails above a few hundred rows; the cache uses paged, filtered
  `/search/` (200 rows, 10 s apart, 5xx retried). A cold single-object GET can return HTTP
  500 after ~40 s; it is transient (retry later). How to handle 5xx on alias lookups is an
  open question the user has pinned — do not change it unasked.
- The `uploads` ledger lags the Portal by design; confirm with a read-only `get_by_alias`
  before concluding it is stale.

## Checklist

- [ ] Read CLAUDE.md and the newest baton; confirm the dataset is multiome and not live.
- [ ] Download (dry run, then real) under the lock; ledger shows all files `done`.
- [ ] QC tables + exploratory plots; give the user the blank threshold sheet.
- [ ] Real guides from the user's thresholds; every `qc_thresholds.tsv` verified.
- [ ] Config generated; `grep -c QC_pseudobulks/` = 0; overrides set if needed.
- [ ] Envs built once (Step 6) if any env file changed.
- [ ] Cache fresh → `default -n` passes → real run with `--retries=1` exits 0.
- [ ] Manifests reviewed against expectations; user informed. Validate only on OK. Never upload.
- [ ] Update the baton (`.baton/<date>-<slug>.md`) with what happened.
