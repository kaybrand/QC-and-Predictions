# check_dataset_accession_mapping.py

Standalone script, no dependency on the rest of this repo. Requires an IGVF
API key/secret pair.

## Credentials

Set as environment variables -- never as a command-line argument, never in a
file, never committed (this repo is public):

```bash
export IGVF_API_KEY=your_key_id_here
export IGVF_SECRET_KEY=your_secret_key_here
export IGVF_MODE=prod   # or staging/sandbox
```

Ask your IGVF DACC data wrangler for a key/secret pair if you don't have
one; READ access is sufficient.

## What it gives you

Queries every `PseudobulkSet` on the IGVF Data Portal and, for each
`dataset` (parsed from a primary pseudobulk's
`{lab}:{dataset}-{cluster}-{subsample}` alias), reports which
`input_file_sets` accession (i.e. principal analysis set accession) it maps
to. Warns by name on stderr if any dataset maps to more than one accession.

## Run it

```bash
pip install igvf-utils
python check_dataset_accession_mapping.py --json-out dataset_to_principal_analysis_set_accession.json
```

Writes `{"dataset": ["ACCESSION", ...], ...}` -- a one-element list per
dataset in the normal case.
