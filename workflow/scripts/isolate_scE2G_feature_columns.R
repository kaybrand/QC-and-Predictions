## Extract scE2G feature columns for sharing with E2G Pillar Project
## Kayla Brand
## July 23, 2025
##
## Forked from IGVF/workflow/scripts/isolate_scE2G_feature_columns.R to accept
## -m/--model instead of hardcoding "multiome_powerlaw_v3" in the header, so
## ATAC-only clusters (scATAC_powerlaw_v3) get a correct header too.

library(optparse)

# Process input arguments --------------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to transcripts input file", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Path to output file", metavar = "character"),
  make_option(c("-c", "--cell_type"), type = "character", default = NULL,
              help = "Cell type", metavar = "character"),
  make_option(c("-d", "--term_id"), type = "character", default = NULL,
              help = "IGVF Sample Term ID", metavar = "character"),
  make_option(c("-s", "--summary"), type = "character", default = NULL,
              help = "Short description of the sample including treatments", metavar = "character"),
  make_option(c("-m", "--model"), type = "character", default = NULL,
              help = "scE2G model that produced these predictions (e.g. multiome_powerlaw_v3, scATAC_powerlaw_v3)", metavar = "character")
)

# parse arguments
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# function to check for required arguments
check_required_args <- function(arg, opt, opt_parser) {
  if (is.null(opt[[arg]])) {
    print_help(opt_parser)
    stop(arg, " argument is required!", call. = FALSE)
  }
}

# check that all required parameters are provided
required_args <- c("input_file", "output_file", "model")
for (i in required_args) {
  check_required_args(i, opt = opt, opt_parser = opt_parser)
}

# Process file -------------------------------------------------------------------------------------

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# load input file
pred <- fread(opt$input_file)

# get all score columns (all columns except EG-pair defining columns)
message("Reformatting predictions...")

# create header lines
header <- c(
  paste("# Source: scE2G", opt$model),
  "# GenomeReference: IGVFDS0280IQAI",
  "# URL: https://github.com/EngreitzLab/scE2G/tree/main",
  "# Assays: 10x Multiome",
  "# SampleAgnostic: False",
  paste("# SampleTermName:", opt$cell_type),
  paste("# SampleTermID:", opt$term_id),
  paste("# SampleSummaryShort:", opt$summary)
)

# extract raw candidate columns
pred <- pred %>%
mutate(ElementChr = chr,
        name = paste0(ElementChr, ":", start, "-", end)) %>%
select(ElementChr,
        ElementStart = start,
        ElementEnd = end,
        ElementName = name,
        ElementClass = class,
        GeneSymbol = TargetGene,
        GeneTSS = TargetGeneTSS,
        GeneEnsemblID = TargetGeneEnsembl_ID,
        isSelfPromoter = isSelfPromoter,
        CellType,
        E2G_Distance = distance,
        everything())

# save to output file
message("Writing to output file...")

output_dir <- dirname(opt$output_file)
# Create the directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (tools::file_ext(opt$output_file) == "gz") {
# Define the output file path and ensure directories exist
tmp_file <- tools::file_path_sans_ext(opt$output_file)

  # save to gzip compressed file
  writeLines(header, con = tmp_file)
  fwrite(pred, file = tmp_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  system2("gzip", args = c("-f", tmp_file))

} else {
  # save to uncompressed file
  writeLines(header, con = opt$output_file)
  fwrite(pred, file = opt$output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
}

message("Done!")
