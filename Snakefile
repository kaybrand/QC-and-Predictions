"""
QC filter IGVF multiome pseudobulks (ATAC + RNA)

Filters ATAC fragment files and RNA count matrices to QC-passing barcodes
for each dataset/cell type combination, then writes a per-dataset config
table for downstream analysis.

Author: Kayla Brand
Date:   Feb 26, 2026

Usage:
    snakemake --configfile config/config_qc_pseudobulks.yaml --cores <N>

Config keys (see config/config_qc_pseudobulks.yaml):
    QC_plots_dir    : directory containing per-cell-type QC barcode lists
    pseudobulk_dir  : root directory of the pseudobulking pipeline output
    out_dir         : root output directory for filtered files
    chrom_sizes     : path to chromosome sizes file (defines sort order)
    datasets        : mapping of dataset -> list of cell types
                      (can be generated with generate_dataset_yaml.py)
"""

import csv
import os
import yaml

configfile: "config/config_QC_pseudobulks.yaml"

# ---------------------------------------------------------------------------
# Directories
# ---------------------------------------------------------------------------
QC_DIR   = config["QC_plots_dir"]   # per-cell-type filtered barcode lists
DATA_DIR = config["pseudobulk_dir"] # pseudobulking pipeline output root
OUT_DIR  = config["out_dir"]        # filtered output root
TRANSCRIPTOME = config["transcriptome"] # GENCODE 43 transcriptome reference

# ---------------------------------------------------------------------------
# Datasets and cell types
# ---------------------------------------------------------------------------
# config["datasets"] is a dict: { dataset: [cell_type, ...], ... }
# Tip: generate this block with generate_dataset_yaml.py and paste into
# config/config_qc_pseudobulks.yaml, or point to a separate yaml and load it:
#     with open(config["datasets_yaml"]) as f:
#         DATASETS = yaml.safe_load(f)["datasets"]
DATASETS = config["datasets"]

# Flat list of all (dataset, cell_type) pairs — used in expand() calls
PAIRS        = [(ds, ct) for ds, cts in DATASETS.items() for ct in cts]
ALL_DATASETS = [p[0] for p in PAIRS]
ALL_CTYPES   = [p[1] for p in PAIRS]


# ---------------------------------------------------------------------------
# Helper: cell types for a given dataset (used in rule make_config_table)
# ---------------------------------------------------------------------------
def cell_types_for(dataset):
    return DATASETS[dataset]

MAX_MEM_MB = 250 * 1000  # 250GB
def determine_mem_mb(wildcards, input, attempt, min_gb=8):
	# Memory resource calculator for snakemake rules
	input_size_mb = input.size_mb
	if ".gz" in str(input):
		input_size_mb *= 8  # assume gz compressed the file <= 8x
	attempt_multiplier = 2 ** (attempt - 1)  # Double memory for each retry
	mem_to_use_mb = attempt_multiplier *  max(4 * input_size_mb, min_gb * 1000)
	return min(mem_to_use_mb, MAX_MEM_MB)


# ---------------------------------------------------------------------------
# rule all — declare final targets to trigger the full workflow
# ---------------------------------------------------------------------------
rule all:
    input:
        # One config table per dataset
        expand(
            os.path.join(OUT_DIR, "config", "tables", "{dataset}_config.tsv"),
            dataset=list(DATASETS.keys())
        ),
        # ATAC index files (presence implies fragment file + index are done)
        expand(
            os.path.join(OUT_DIR, "{dataset}", "{cell_type}",
                         "atac_fragments_{dataset}_{cell_type}.tsv.gz.tbi"),
            zip, dataset=ALL_DATASETS, cell_type=ALL_CTYPES
        ),
        # RNA matrix files
        expand(
            os.path.join(OUT_DIR, "{dataset}", "{cell_type}",
                         "rna_count_matrix_{dataset}_{cell_type}", "matrix.mtx.gz"),
            zip, dataset=ALL_DATASETS, cell_type=ALL_CTYPES
        ),


# ---------------------------------------------------------------------------
# rule fragment_file — filter ATAC fragments for one dataset/cell_type
# ---------------------------------------------------------------------------
rule atac_fragment_file:
    input:
        barcodes   = os.path.join(QC_DIR, "{dataset}", "{cell_type}", "filtered_barcodes_with_subsamples.tsv.gz"),
        pseudobulks = os.path.join(DATA_DIR, "{dataset}", "pseudobulks"),
    output:
        filtered_ATAC = os.path.join(
            OUT_DIR, "{dataset}", "{cell_type}",
            "atac_fragments_{dataset}_{cell_type}.tsv.gz"
        ),
        index = os.path.join(
            OUT_DIR, "{dataset}", "{cell_type}",
            "atac_fragments_{dataset}_{cell_type}.tsv.gz.tbi"
        ),
    params:
        chrom_sizes = config["chrom_sizes"],
    conda:
        "workflow/envs/filter_multiome_env.yaml"
    resources:
        mem_mb = determine_mem_mb
    shell:
        """
        python workflow/scripts/filter_atac_fragments.py \
            --qc-guide    {input.barcodes} \
            --pseudobulks {input.pseudobulks} \
            --cell-type   {wildcards.cell_type} \
            --chrom-sizes {params.chrom_sizes} \
            --out         {output.filtered_ATAC} 
        """


# ---------------------------------------------------------------------------
# rule rna_count_matrix — filter RNA counts for one dataset/cell_type
# ---------------------------------------------------------------------------
rule rna_count_matrix:
    input:
        barcodes    = os.path.join(QC_DIR, "{dataset}", "{cell_type}", "filtered_barcodes_with_subsamples.tsv.gz"),
        pseudobulks = os.path.join(DATA_DIR, "{dataset}", "pseudobulks"),
    output:
        filtered_RNA = os.path.join(
            OUT_DIR, "{dataset}", "{cell_type}",
            "rna_count_matrix_{dataset}_{cell_type}", "matrix.mtx.gz"
        ),
        out_barcodes = os.path.join(
            OUT_DIR, "{dataset}", "{cell_type}",
            "rna_count_matrix_{dataset}_{cell_type}", "barcodes.tsv.gz"
        ),
        out_genes    = os.path.join(
            OUT_DIR, "{dataset}", "{cell_type}",
            "rna_count_matrix_{dataset}_{cell_type}", "features.tsv.gz"
        ),
    params:
        # filter_rna_counts.py derives the output directory from this path
        rna_out = os.path.join(
            OUT_DIR, "{dataset}", "{cell_type}",
            "rna_count_matrix_{dataset}_{cell_type}.mtx"
        ),
        # gtf file to map Ensembl IDs to gene symbols
        transcriptome = TRANSCRIPTOME
    conda:
        "workflow/envs/filter_multiome_env.yaml"
    resources:
        mem_mb = determine_mem_mb
    shell:
        """
        python workflow/scripts/filter_rna_counts.py \
            --qc-guide    {input.barcodes} \
            --pseudobulks {input.pseudobulks} \
            --cell-type   {wildcards.cell_type} \
            --out         {params.rna_out} \
            --gtf {params.transcriptome} \
            --standard-chromosomes-only
        """


# ---------------------------------------------------------------------------
# rule make_config_table — write per-dataset TSV for downstream tools
# ---------------------------------------------------------------------------
rule make_config_table:
    input:
        atacs = lambda wc: expand(
            os.path.join(OUT_DIR, wc.dataset, "{cell_type}",
                         "atac_fragments_" + wc.dataset + "_{cell_type}.tsv.gz"),
            cell_type=cell_types_for(wc.dataset)
        ),
        rnas  = lambda wc: expand(
            os.path.join(OUT_DIR, wc.dataset, "{cell_type}",
                         "rna_count_matrix_" + wc.dataset + "_{cell_type}", "matrix.mtx.gz"),
            cell_type=cell_types_for(wc.dataset)
        ),
    output:
        table = os.path.join(OUT_DIR, "config", "tables", "{dataset}_config.tsv"),
    run:
        os.makedirs(os.path.dirname(output.table), exist_ok=True)
        cell_types = cell_types_for(wildcards.dataset)

        # Map cell_type -> paths for convenient lookup
        atac_map = {
            ct: os.path.join(OUT_DIR, wildcards.dataset, ct,
                             f"atac_fragments_{wildcards.dataset}_{ct}.tsv.gz")
            for ct in cell_types
        }
        rna_map  = {
            ct: os.path.join(OUT_DIR, wildcards.dataset, ct,
                             f"rna_count_matrix_{wildcards.dataset}_{ct}")
            for ct in cell_types
        }

        headers = [
            "cluster", "rna_matrix_file", "atac_frag_file",
            "HiC_file", "HiC_type", "HiC_resolution",
            "alt_TSS", "alt_genes", "model_dir"
        ]
        with open(output.table, "w", newline="") as tsvfile:
            writer = csv.DictWriter(tsvfile, fieldnames=headers, delimiter="\t")
            writer.writeheader()
            for ct in cell_types:
                writer.writerow({
                    "cluster":         ct,
                    "rna_matrix_file": rna_map[ct],
                    "atac_frag_file":  atac_map[ct],
                    "model_dir":       "models/multiome_powerlaw_v3,models/scATAC_powerlaw_v3",
                })
            
 