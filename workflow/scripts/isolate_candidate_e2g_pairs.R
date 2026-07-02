## Extract candidate E2G links for sharing
## June 17, 2025

library(optparse)

# Process input arguments --------------------------------------------------------------------------

# create arguments list
option_list = list(
  make_option(c("-i", "--input_file"), type = "character", default = NULL,
              help = "Path to transcripts input file", metavar = "character"),
  make_option(c("-o", "--output_file"), type = "character", default = NULL,
              help = "Path to output file", metavar = "character"),
  # make_option(c("-c", "--cell_type"), type = "character", default = NULL,
  #             help = "Cell type", metavar = "character"),
  # make_option(c("-d", "--term_id"), type = "character", default = NULL,
  #             help = "IGVF Sample Term ID", metavar = "character"),
  make_option(c("-s", "--summary"), type = "character", default = NULL,
              help = "Short description of the sample including treatments", metavar = "character")
  # make_option(c("-v", "--version"), type = "character", default = NULL,
  #             help = "E2G method version", metavar = "character"),
  # make_option(c("-l", "--portal_link"), type = "character", default = NULL,
  #             help = "Link to metadata of this file on IGVF data portal", metavar = "character")
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
required_args <- c("input_file", "output_file", "summary")
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

# set summarized sample description if specified
if (!is.null(opt$summary)) {
  pred$SampleSummaryShort <- opt$summary
}

# create header lines
# header <- c(
#   paste("# Source: ABC"),
#   # paste("# Version:", opt$version),
#   "# GenomeReference: IGVFDS0280IQAI",
#   "# URL: https://abc-enhancer-gene-prediction.readthedocs.io/en/latest/usage/methods.html#defining-candidate-elements",
#   "# Assays: 10x Multiome",
#   "# SampleAgnostic: False",
#   paste("# SampleTermName:", opt$cell_type),
#   paste("# SampleTermID:", opt$term_id),
#   paste("# SampleSummaryShort:", opt$summary)
# )

# # add link to where the metadata of the file is stored on the IGVF data portal if available
# if (!is.null(opt$portal_link)) {  
#     header <- c(header, "# Metadata:", opt$portal_link)
# }

# extract raw candidate columns
pred <- pred %>%
mutate(ElementChr = chr,
        name = paste0(ElementChr, ":", start, "-", end)) %>%
select(ElementChr,
        ElementStart = start, ElementEnd = end, 
        ElementName = name, ElementClass = class,  
        GeneTSS = TargetGeneTSS, GeneSymbol = TargetGene,
        GeneEnsemblID = TargetGeneEnsembl_ID,
        SampleSummaryShort, E2G_Distance = distance, 
        isSelfPromoter = isSelfPromoter)

# save to output file
message("Writing to output file...")
if (tools::file_ext(opt$output_file) == "gz") {
  
  # save to gzip compressed file
  tmp_file <- tools::file_path_sans_ext(opt$output_file)
  # writeLines(header, con = tmp_file)
  fwrite(pred, file = tmp_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  system2("gzip", args = c("-f", tmp_file))
  
} else {
  
  # save to uncompressed file
  # writeLines(header, con = opt$output_file)
  fwrite(pred, file = opt$output_file, sep = "\t", quote = FALSE, na = "NA", append = TRUE,
         col.names = TRUE)
  
}

message("Done!")
