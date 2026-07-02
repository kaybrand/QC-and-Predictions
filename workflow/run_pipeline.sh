#!/bin/bash
# Wrapper for the QC-and-Predictions pipeline.
#
# Required because scE2G's own workflow/Snakefile computes
# SCRIPTS_DIR = os.path.join(os.getcwd(), "workflow", "scripts") -- it assumes
# it is always run from its own repo root. Importing it as a Snakemake `module`
# from this pipeline's Snakefile does NOT change that: os.getcwd() is a real
# process-level property, unaffected by which Snakefile Snakemake was told to
# run. If invoked from anywhere else, scE2G's model-application scripts
# silently resolve to the wrong repo and fail deep inside its own rules. This
# wrapper cd's into scE2G_dir (read from the given config) before invoking
# snakemake, so you never have to remember this yourself. (A proper fix
# belongs upstream in scE2G -- see the companion task doc for that PR.)
#
# Usage:
#   mamba activate run_snakemake9
#   workflow/run_pipeline.sh workflow/config/{name}_pipeline_config.yaml [extra snakemake args...]
#
# Example (dry run):
#   workflow/run_pipeline.sh workflow/config/igvf10_pipeline_config.yaml -n -p
# Example (real run):
#   workflow/run_pipeline.sh workflow/config/igvf10_pipeline_config.yaml --executor slurm --profile slurm.smk9 --use-conda -p

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <pipeline_config.yaml> [extra snakemake args...]" >&2
    exit 1
fi

CONFIG_FILE=$(readlink -f "$1")
shift

SCE2G_DIR=$(python -c "import sys, yaml; print(yaml.safe_load(open(sys.argv[1]))['scE2G_dir'])" "$CONFIG_FILE")
THIS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

cd "$SCE2G_DIR"
exec snakemake -s "$THIS_DIR/Snakefile" --configfile "$CONFIG_FILE" "$@"
