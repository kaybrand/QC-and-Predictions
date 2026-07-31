# IGVF Cell Annotation Report -- User Guide

`get_igvf_cell_annotation_report.py` queries the IGVF Data Portal directly
and writes one TSV row per `{dataset}-{cluster}`, covering every primary
pseudobulk (`PseudobulkSet`) currently on the portal: CellAnnotation,
SampleTermID, SampleTermName, CellQualifier, contributing subsamples, and
release status.

It needs only an IGVF API key/secret pair -- no access to this lab's
pipeline, config files, or anything else, in its default mode. An optional
second mode additionally takes this lab's own per-cluster QC guide
directory, if you have it, to resolve a single value per cluster instead of
reporting every value that disagrees. **Pick a specific commit before you
run either mode** -- see below.

## Which mode should I use?

The script has exactly one required input (an IGVF key pair) and one
optional input (`--qc-guide-dir`). Both modes write the same 9 columns; the
only difference is how much they can resolve.

| | No `--qc-guide-dir` | `--qc-guide-dir <path>` |
|---|---|---|
| Works for | Everyone with an IGVF key pair | Only clusters this lab has shared a QC guide for |
| `CellAnnotation`/`SampleTermID`/`SampleTermName`/`CellQualifier` | Every distinct value seen, `" | "`-joined -- a disagreement report | A single resolved value, wherever a matching guide is found |
| `Subsamples` order | Alphabetical (no contribution info available) | Descending by cell count, from the guide |
| `QCGuideFile` | Always blank | The guide's filename, when one was used to resolve that row -- **blank on any row that mode couldn't resolve** (no matching guide found, or the row's own contributing primaries disagreed on SampleTermID/SampleTermName, which should never happen and isn't trusted if it does) |

You always get a full, valid table either way -- passing `--qc-guide-dir`
never causes the script to fail or skip a dataset-cluster; it only adds
resolved values where it can. **`QCGuideFile` tells you, per row, which
behavior actually applied** -- blank means "disagreement reported, not
resolved," non-blank means "resolved from that guide file."

## Pick a specific commit

This repo (and this script) changes over time. Before you run it, pick a
commit and stick with it, rather than tracking a moving branch -- that way
your results are reproducible and you know exactly what logic produced
them:

```bash
git clone https://github.com/kaybrand/QC-and-Predictions.git
cd QC-and-Predictions
git log --oneline -- igvf_cell_annotation_report/    # find a commit you're happy with
git checkout <commit-sha>
```

Record `<commit-sha>` alongside any output you share from this script.

## Requirements

- Python 3.8+
- The `igvf-utils` package, which also pulls in `requests`:

  ```bash
  pip install igvf-utils
  ```

## Credentials

**This repository is public.** The script authenticates via `igvf_utils`,
which reads three environment variables -- it never takes credentials as
command-line arguments, never reads them from a file, and never writes them
anywhere. Never commit a key/secret into this repo (or any fork/copy of
it), a script, a ticket, or a chat message.

```bash
export IGVF_API_KEY=your_key_id_here
export IGVF_SECRET_KEY=your_secret_key_here
export IGVF_MODE=prod   # or: staging, sandbox
```

- Get a key/secret pair from your IGVF DACC data wrangler if you don't
  already have one.
- `IGVF_MODE=prod` points at the real production portal
  (`https://api.data.igvf.org/`); use `staging`/`sandbox` only if you were
  specifically given credentials for one of those instances.
- Setting these in your shell (as above), or in a local `.env`-style file
  you `source` yourself and never `git add`, are both fine. If you keep such
  a file inside this folder, double check `git status` before committing
  anything -- don't rely solely on `.gitignore`.

## Usage

Mode 1 -- works for anyone with a key pair, reports disagreement rather than
resolving it:

```bash
python get_igvf_cell_annotation_report.py
```

Writes `cell_annotations_by_dataset_cluster.tsv` in the current directory.
Use `-o/--output` to write somewhere else.

Mode 2 -- additionally resolves a single value per cluster, wherever a
matching QC guide is found:

```bash
python get_igvf_cell_annotation_report.py --qc-guide-dir /path/to/qc_guides
```

`--qc-guide-dir` must point at a directory structured
`{qc-guide-dir}/{dataset}/{cluster}/filtered_barcodes_with_subsamples.tsv.gz`
(this lab's standard per-cluster filtered-barcode QC guide filename and
layout -- e.g. `.../igvf1/hap1/filtered_barcodes_with_subsamples.tsv.gz`).
If that exact filename isn't present for a given `{dataset}/{cluster}`, but
exactly one other `.tsv`/`.tsv.gz` file is, that one is used instead; if
neither condition holds (no matching directory, no file, or more than one
ambiguous candidate), that row falls back to mode 1's behavior -- the run
never fails because one cluster's guide is missing or oddly named.

Note: a cluster's directory name under `--qc-guide-dir` must match its name
on the portal exactly (e.g. `hl60`, not some renamed variant like
`hl60_concordant_subsamples`) to be found -- a mismatch there falls back to
mode 1 for that cluster, same as a genuinely missing guide.

The script prints a short progress/summary log to stderr, including which
mode is active, how many rows were resolved via a guide, and how many still
have a disagreement flag.

## Output columns

| Column | Meaning |
|---|---|
| `Dataset_Cluster` | The `{dataset}-{cluster}` this row summarizes |
| `Subsamples` | `" | "`-joined list of contributing subsample accessions -- alphabetical in mode 1, descending-by-cell-count (only subsamples contributing >=1 cell) in mode 2 |
| `N_Subsamples` | Count of subsamples in the `Subsamples` column |
| `CellAnnotation` | Resolved value (mode 2, when `QCGuideFile` is set) or every distinct value seen, `" | "`-joined (otherwise); >1 value means disagreement |
| `SampleTermID` | Resolved or distinct portal `cell_type.term_id` value(s) (CURIE form, e.g. `CL:0000235`) |
| `SampleTermName` | Resolved or distinct portal `cell_type.term_name` value(s) (e.g. `macrophage`) |
| `CellQualifier` | Resolved or distinct `cell_qualifier` value(s); often blank |
| `Status` | `released` if every contributing primary pseudobulk is `released`, else `in progress` |
| `QCGuideFile` | The QC guide filename used to resolve this row, or blank if this row reports disagreement instead (always blank in mode 1) |

Multi-value columns are joined with `" | "`, not a comma: several of these
are free-text portal fields that routinely contain a literal comma
themselves, e.g. a `SampleTermName` of `"B cell, CD19-positive"` -- joining
with a comma would make a single consistent value indistinguishable from
two disagreeing ones.

## Troubleshooting

- **`KeyError: 'IGVF_MODE'`** -- you haven't set the `IGVF_MODE` environment
  variable (and didn't pass `--igvf-mode`). Set it as shown above.
- **401/403 response** -- your API key/secret is missing, wrong, or expired.
  Double-check the exported values and confirm the key pair is active with
  your DACC data wrangler.
- **A certificate/TLS verification error** -- this is a real error (a
  network/proxy/cert issue), not something to work around by disabling
  verification. Investigate the underlying cause instead.
- **A cluster I expected to resolve has a blank `QCGuideFile`** -- check the
  stderr log for that dataset/cluster: either no `{dataset}/{cluster}`
  directory was found under `--qc-guide-dir`, its guide file wasn't the
  expected name and wasn't unambiguous, or its contributing primaries
  disagreed on `SampleTermID`/`SampleTermName` (logged as a `WARNING` --
  this should never happen and isn't trusted if it does).
