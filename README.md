# QC-and-Predictions

End-to-end Snakemake pipeline for the IGVF E2G Pillar Project: takes one or
more QC-filtered clusters from pseudobulk to shareable scE2G data products.
On this branch (`igvf-portal-submission`) the primary destination is the IGVF
Data Portal; the Synapse path runs alongside it, not instead of it.

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

One command per dataset, through the driver:

```bash
sbatch resources/run_pipeline.sbatch {dataset}              # default: everything
sbatch resources/run_pipeline.sbatch {dataset} local_only   # no IGVF Portal contact at all
sbatch resources/run_pipeline.sbatch {dataset} default -n   # dry run first
```

Several datasets at once is safe and needs no ordering between them:

```bash
for ds in igvf0 igvf1 igvf2; do sbatch resources/run_pipeline.sbatch "$ds"; done
```

What makes that safe is an `flock` on `{state_db}.warmlock`, held by every
stage that touches `state.db` — the warm stage *and* stage 3, because
`orchestrator.run()` calls `refresh_if_stale` and so writes too. It is
deliberately **not** held across stage 2, which is hours of scE2G compute and
would serialise every dataset. Verified 2026-08-25: eight drivers launched at
once across four hosts, three of them requesting the lock in the same second
and being granted it in strict sequence (3 s and 4 s waits, each acquire
landing exactly on the previous release), `PRAGMA integrity_check` clean
afterwards. Cross-host locking works because Oak is mounted with `flock`, not
`localflock` — worth re-checking with `mount | grep oak` if `state.db` ever
moves to another filesystem, since `localflock` would make locks node-local
and let two hosts write WAL concurrently while both logged success.

`workflow/scripts/run_pipeline.py` runs five stages in order — preflight, cache
warm, Snakemake, manifest preview, audit — so producing a dataset's IGVF Portal
manifests no longer requires three commands in a specific undocumented order.
It bakes in `--conda-prefix` and `--rerun-triggers mtime` (omitting either has
caused real multi-hour incidents) and never passes `--omit-from`.

**Exit 0 means complete and verified, and nothing else does.** A portal outage
degrades the run rather than aborting it: all scE2G work still completes, the
manifest stages are skipped, and the exit code is non-zero with an itemised
reason. The driver can never upload — `--manifest-mode` accepts only `preview`
and `validate`.

| Exit | Meaning |
|------|---------|
| `0` | every manifest-eligible cluster resolved and every expected row is present |
| `1` | output was produced but something is incomplete — portal unreachable, a cluster's annotation unresolved, a manifest gap. Always itemised as `DEGRADED:` / `PROBLEM:` lines |
| `2` | preflight/config error — nothing ran |

Exit 1 is routine and not necessarily a failure to chase: a cluster that is
manifest-eligible but has no resolvable CellAnnotation produces one `PROBLEM:`
line and exit 1 while every other cluster in the dataset completes normally.
Read the itemised list before concluding the run was wasted.

### Environments

Two conda envs, on purpose, and you need both:

- **`run_snakemake9`** — has `snakemake`, `conda`, `mamba`. No `igvf_utils`.
- **`igvf_utils_env`** — has `igvf_utils` (required for the Portal GET) and
  SQLite 3.50. No `snakemake`.

`resources/run_pipeline.sbatch` activates `run_snakemake9` (which puts
`snakemake` on `PATH` for `--use-conda`) and then runs the driver under
`igvf_utils_env`'s interpreter; the driver shells out to `snakemake` from
`PATH`, so both halves get the interpreter they need. Running any of this
Python by hand means picking the right interpreter yourself — in particular
`state.py`'s UPSERTs need SQLite ≥ 3.24 and the login node's system Python
ships 3.7.17, which fails in a way that does not obviously point at SQLite.

### Two modes

- **`default`** — everything, including the live IGVF Portal GET that populates
  the CellAnnotation cache, portal-format reformatting, and manifest generation.
- **`local_only`** — everything that needs no Portal contact whatsoever
  (filtered fragments/RNA matrices, all scE2G predictions/candidates/features,
  QC report, IGV tracks, EnhancerList and RNA matrix packaging). It never opens
  `state.db`, never reads the CellAnnotation snapshot, and requests no
  reformat target — so Portal availability cannot affect it. Set per-run with
  `pipeline_mode` in the config or `--config pipeline_mode=local_only`.

Running `snakemake` directly still works, but in `default` mode it aborts before
the first job if the driver hasn't warmed the cache, rather than silently
producing zero reformat targets the way it used to.

Before committing real compute to a dataset you haven't run, you can answer
"how many of its clusters will even resolve?" offline, with zero network
contact and against a *copy* of `state.db`, so a diagnostic never mutates the
shared ledger:

```bash
sbatch resources/predict_fanout_scopes.sbatch igvf1 igvf2 igvf3   # wraps verify_derive_scopes.py
```

It reports `N/M scope(s) resolved` per dataset plus the exact reason for each
one that didn't (`no_matching_primary_alias: ...`, `unreadable_qc_guide: ...`).
The numbers are a lower bound — a fresh Portal fetch can only add rows — and
in practice they have matched the real run's derive counts exactly.

Always dry-run first (`-n`) and inspect the plan, especially before the
first run for a new dataset. Synapse uploads and stale-entity deletions
default to dry-run/disabled (`synapse.dry_run: true`,
`synapse.confirm_delete: false`, `synapse.confirm_overwrite: false` in the
config) — flip these deliberately, one at a time, after reviewing what the
pipeline reports it would do.

### Knowing what happened

- `{output_dir}/report.tsv` — one row per cluster: quality gate, exclusion
  reason, CellAnnotation, predictions, reformatting, and the manifest roll-up
  (`manifest_status`, `manifest_gap_reason`). Answers "which clusters are in,
  and why isn't this one" without opening a manifest.
- `{output_dir}/igvf_manifests/{dataset}/manifest_coverage.tsv` — the per-row
  detail behind that: one row per `(cluster, table, variant)` with its outcome
  (`planned-post`, `skipped-missing-file:<path>`, `invalid:<msg>`, …).
- `{output_dir}/igvf_metadata/{dataset}_cell_annotation_status.tsv` — per
  cluster, whether its CellAnnotation resolved and the exact reason if not.

See `workflow/config/example_pipeline_config.yaml` for the full config
schema and inline documentation of every field.

### Troubleshooting: a run wants to recompute Kendall for everything

If a dry run plans `generate_atac_matrix` → `compute_kendall` → `arc_e2g` for
clusters whose outputs already exist and are current, check the *directory*
mtime before letting it run. Measured 2026-08-25: **97 of 115 clusters** across
nine datasets were stale this way, for byte-identical output.

The cause is upstream and is not a staleness signal at all. scE2G's
`generate_atac_matrix` declares `cell_barcodes_path = RESULTS_DIR` — the whole
`{output_dir}/uniformly_processed/{dataset}` directory — as an input whenever
`RNA_matrix_filtered: True` (scE2G's default). `generate_atac_matrix.R` then
guards it with `if (file_test("-f", cell_bc_path))`, which is false for a
directory, so **the input is never read and its mtime carries no correctness
information.** But a directory's mtime bumps whenever an entry is added or
removed in it — and `rules/qc_stats.smk` creates
`{dataset}/qc_plots/` there. Under `--rerun-triggers mtime` that invalidates
every cluster's `generate_atac_matrix` and cascades.

It presents as intermittent because it is pure ordering: the bump only matters
if it lands *after* the outputs were written. Writing files *inside* an
existing subdirectory never propagates to the parent.

Diagnose and fix:

```bash
cd {output_dir}/uniformly_processed
stat -c %y igvf4                                    # directory mtime
find igvf4 -mindepth 3 -maxdepth 3 -name atac_matrix.rds -printf '%TY-%Tm-%Td %TH:%TM\n' | sort | head -1
# stale if the directory is NEWER than the oldest atac_matrix.rds:
oldest=$(find igvf4 -mindepth 3 -maxdepth 3 -name atac_matrix.rds -printf '%T@\n' | sort -n | head -1 | cut -d. -f1)
touch -h -d "@$((oldest-60))" igvf4
```

Reset the *directory* mtime, not the outputs — ten metadata writes rather than
rewriting mtimes on ~97 large files, which is what `snakemake --touch` would do
and which can mask something genuinely stale. Confirm with a dry run: the job
count should drop to just the genuinely-missing work. Note the reset is not
permanent — it recurs whenever a new top-level entry appears under a dataset
directory, and it has been observed to revert. **Verify with a dry run
immediately before each real run rather than assuming a previous fix still
holds.** `qc_plots/` cannot simply be moved elsewhere: scE2G's still-imported
`hover_plots` rule reads that exact path (see `rules/qc_stats.smk`'s docstring).

The durable fix is one word upstream — `return ancient(RESULTS_DIR)` — and
belongs in a PR to `EngreitzLab/scE2G`, not here.

## IGVF Data Portal uploads

`igvf-portal-submission` branch only (not present on the Synapse-only
`synapse-submission` branch). Manifest generation is now sequenced by
`workflow/scripts/run_pipeline.py` (see Running above) rather than being a
standalone step run by hand; the CLI below is still the way to do a one-off
pass, and is the ONLY way to upload for real. Eleven metadata tables are
registered (Prediction Set, Prediction Tabular Files, Principal Pseudobulk
Set, Filtered Barcode Lists, Filtered ATAC Fragment Files, Filtered Matrix
Files, Signal Files, ATAC Index File, Elements BED Index File, BEDPE Index
File, Documents).

Eleven tables, but **fifteen** manifest files per fully-covered dataset: a
table can register several *variants*, each of which gets its own file.
Prediction Tabular Files alone has five (`elements_bed`, `genes`, `full`,
`thresholded`, `bedpe`). Don't read "11 tables" and "15 files" as a
discrepancy — see [Reading the manifest directory](#reading-the-manifest-directory).

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

### Reading the manifest directory

`{output_dir}/igvf_manifests/{dataset}/` holds two different kinds of file, and
mistaking one for the other is the easiest way to misread a run.

**`round{N}_{table}_{variant}_{post,patch}.tsv` — the submission files.**
Ephemeral, and rewritten every run. The round number encodes dependency order:
**sorting by filename *is* the upload order**, because a later round's payloads
reference aliases earlier rounds create. A fully-covered dataset produces
fifteen of these.

- `_post` — no object exists at this alias; create it.
- `_patch` — it exists and the payload hash *changed*; update it.
- A table can legitimately have **neither**. If a row's alias is already live
  and its payload hash is unchanged, the outcome is `unchanged` and nothing is
  written, because a patch would ask the Portal to overwrite fields with the
  values they already hold. **Fourteen round files instead of fifteen is not
  automatically a failure** — check `manifest_coverage.tsv` for `unchanged`
  before treating it as a gap.

Writing a PATCH retires a stale sibling POST; writing a POST does **not** retire
a PATCH — it warns and keeps both. That asymmetry is deliberate. `plan_table`'s
alias lookup is the Elasticsearch-backed read, which can lag a write and
falsely report "not found", so the POST-while-PATCH-exists direction can be
fabricated spuriously and must never silently delete real work.

**`{object_type}.tsv` — the durable accumulator.** Named for the *Portal's*
object type, not our table name, which is why nine `QC_documents` rows land in
`document.tsv`. Several of our tables share one object type (`filtered_barcode_list`,
`filtered_atac_fragment_file` and every `prediction_tabular_files` variant are all
`tabular_file`), and grouping this way is what makes a future bulk field edit
across a whole object type one edit-and-resubmit against one file, matching
`iu_register.py`'s one-profile-per-invocation constraint. It carries `record_id`s,
never shrinks, and is **never handed to the registrar** — it is a reference, not
a submission.

One consequence worth knowing before you go hunting: the `unchanged` branch
short-circuits on the local ledger **before** the live alias lookup. So a row
can be skipped with zero Portal contact on the strength of `state.db` alone. If
objects were deleted Portal-side, the ledger will not know. Confirm against the
Portal with a read-only `get_by_alias` before concluding the ledger is stale —
and equally, before "fixing" it, since resetting a ledger row whose object is
genuinely still live just generates a spurious PATCH of identical values.

### One wholesale GET, then per-dataset derivation

`cell_metadata.fetch_if_stale` issues a single multireport GET covering EVERY
primary pseudobulk on the portal (24h TTL, tracked by one global row) and saves
all of it to `state.db`'s raw `cell_metadata_primary_pseudobulks` /
`cell_metadata_principal_pseudobulks` tables. `cell_metadata.derive_scopes` then
turns that raw cache into the per-`(dataset, cluster)` `cell_annotations` view,
using each cluster's own `qc_guide` to pick a winner where contributing
primaries disagree on Cell Annotation — which is why derivation is per-config
and can't be prepopulated wholesale.

Those two used to be one function gated on the same TTL, which was a real bug:
warming for one dataset made the TTL fresh, so every later dataset returned
early before deriving anything and stayed permanently uncached. Split apart,
each dataset derives its own rows offline from whatever the raw cache already
holds, and no dataset depends on another running first.

### Reformat eligibility is gated on the CellAnnotation cache, not a preview TSV

`reformat_predictions`/`reformat_predictions_thresholded`/`reformat_element_list`/
`reformat_gene_list` (the rules that turn scE2G's own output into portal-format
Prediction Tabular Files; the single `{meta}`-wildcarded `reformat_lists` rule was
split into the element/gene pair on 2026-08-19, see `rules/reformat.smk`)
only run for a `(dataset, cluster)` with a real, cached
CellAnnotation row already in `state.db` — computed at Snakemake parse time
as `common.smk`'s `REFORMAT_ELIGIBLE_CLUSTERS` (read from the driver-written
snapshot described below), keyed through
`cell_annotations.py`'s `annotation_lookup_key()` (which resolves an
ATAC-only variant cluster's real CellAnnotation under its base, non-suffixed
name). **A quality-passing cluster with no cached CellAnnotation yet still
gets full scE2G predictions/candidates/features** — it's simply not
reformatted into portal format until `manage_igvf_metadata.py` has warmed
the cache for it. Core file generation is never blocked by Portal metadata
availability, and this is intentional: `state.db`'s cache is populated only
by the driver's warm stage, never by a Snakemake rule.

**No rule under `workflow/rules/` opens `state.db`.** `common.smk` used to, at
parse time — and every Slurm worker re-parses the workflow, so several
concurrent dataset drivers meant worker-node reads on many hosts overlapping
another driver's writes. `state.db` is on Lustre in WAL mode, whose
shared-memory index is only supported on a single host. Snakemake now reads
`{output_dir}/igvf_metadata/{dataset}_cell_annotations.tsv`, a temp snapshot
the driver writes immediately beforehand and deletes afterwards, carrying the
portal fetch timestamp and a digest of the cluster set — both enforced on read,
so a stale or foreign file raises instead of quietly reporting nothing
annotated. See `workflow/scripts/cell_annotation_snapshot.py`.

### Portal-facing packaging generated directly by this pipeline

A few files exist purely so `manage_igvf_metadata.py`'s uploader tables have
something to find on disk — none of these require Portal contact to
produce:

- `Neighborhoods/EnhancerList.bed.gz` + `.tbi` (`rules/reformat.smk`, wired up in
  `rules/common.smk`) — bgzip/tabix of scE2G's own (already coordinate-sorted)
  `EnhancerList.bed`.
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
