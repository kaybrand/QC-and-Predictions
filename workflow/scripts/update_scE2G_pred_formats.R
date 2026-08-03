## Reformat pillar project predictors
## April 24, 2025

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
  make_option(c("-m", "--method"), type = "character", default = "scE2G",
              help = "E2G method that produced the predictions", metavar = "character"),
  make_option(c("-v", "--version"), type = "character", default = NULL,
              help = "E2G method version", metavar = "character"),
  make_option(c("--threshold"), type = "character", default = NULL,
              help = "Used score threshold if applicable", metavar = "character"), 
  make_option(c("-l", "--portal_link"), type = "character", default = NULL,
              help = "Link to metadata of this file on IGVF data portal", metavar = "character"), 
  make_option(c("-a", "--all_columns"), action = "store_true", default = FALSE,
              help = "Include all columns in element/gene lists")
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
required_args <- c("input_file", "output_file", "cell_type", "method", "version")
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

# set cell annotation if specified
if (!is.null(opt$summary)) {
  pred$CellAnnotation <- opt$summary
}

# create header lines
header <- c(
  paste("# Source:", opt$method),
  paste("# Version:", opt$version),
  "# GenomeReference: IGVFDS0280IQAI",
  "# URL: https://github.com/EngreitzLab/scE2G/tree/main",
  "# Assays: 10x Multiome",
  "# SampleAgnostic: False",
  paste("# SampleTermName:", opt$cell_type),
  paste("# SampleTermID:", opt$term_id),
  paste("# CellAnnotation:", opt$summary)
)

# add threshold if applicable
if (!is.null(opt$threshold)) {
  header <- c(header, paste("# ScoreThreshold:", opt$threshold))
}

# add ScoreType if these are predictions (not metadata List)
if (!grepl("_list", opt$input_file)) {
  header <- c(header, "# ScoreType: positive_score")
}

# add link to where the metadata of the file is stored on the IGVF data portal if available
# (opt$portal_link is the file's own IGVF alias; spelled out into a full portal URL here)
if (!is.null(opt$portal_link)) {
    header <- c(header, paste0("# Metadata: https://data.igvf.org/tabular-files/", opt$portal_link))
}

# add additional columns and extract output columns
if (grepl("element_list", opt$input_file) & !(opt$all_columns)) { # if file is an element list
  pred <- pred %>%
    mutate(ElementChr = chr,
           name = paste0(ElementChr, ":", start, "-", end)) %>%
    select(ElementChr,
           ElementStart = start,
           ElementEnd = end,
           ElementName = name,
           ElementClass = class)
} else if (grepl("element_list", opt$input_file) & (opt$all_columns)) { # include all element list columns
  pred <- pred %>%
    rename(
      ElementChr = chr,
      ElementStart = start,
      ElementEnd = end,
      ElementClass = class
    ) %>%
    mutate(
      ElementName = paste0(ElementChr, ":", ElementStart, "-", ElementEnd)
    ) %>%
    select(ElementChr, ElementStart, ElementEnd, ElementName, ElementClass, everything())
} else if (grepl("gene_list", opt$input_file) & !(opt$all_columns)) { # if file is a gene list
  pred <- pred %>%
    mutate(
        TSSStart = tss-250,
        TSSEnd = tss+250
    ) %>%
    select(
        TSSChr = chr,
        TSSStart,
        TSSEnd,
        GeneSymbol = name,
        GeneEnsemblID = Ensembl_ID,
        GeneStrand = strand)
} else if (grepl("gene_list", opt$input_file) & (opt$all_columns)) { # include all gene list columns
  pred <- pred %>%
    mutate(
        TSSStart = tss-250,
        TSSEnd = tss+250
    ) %>%
    rename(TSSChr = chr,
           TSS = tss,
           GeneSymbol = name,
           GeneEnsemblID = Ensembl_ID,
           GeneStrand = strand) %>% 
    select(TSSChr, TSSStart, TSSEnd, GeneSymbol, GeneEnsemblID, GeneStrand, everything())
} else if (grepl("ATAC", opt$method) & is.null(opt$threshold)) {
  # scATAC unthresholded ("full") predictions -- wider column set (2026-07-23):
  # every model-feature column gets a "feature." prefix, ABC.Score (scATAC's
  # own final integrative score, in place of Multiome's ARC.E2G.Score -- ARC
  # additionally integrates Kendall, which has no scATAC equivalent) is the
  # last column in that prefixed block, and ElementName keeps the raw "name"
  # column's own "class|chrN:start-end" value as-is rather than the narrower
  # thresholded format's lossy chrN:start-end-only reconstruction below. Score
  # sits right after CellAnnotation per the consortium-standard column order
  # (2026-07-24) -- every remaining column is "[Optional columns]", kept in
  # the same relative order as before, just moved after Score. RNA_pseudobulkTPM/
  # RNA_meanLogNorm/RNA_percentCellsDetected are dropped (2026-07-24): those come
  # from an RNA modality that scATAC-only runs don't have, and scE2G doesn't
  # populate them here unless Multiome was also run for this cluster.
  pred <- pred %>%
    select(ElementChr = chr,
           ElementStart = start,
           ElementEnd = end,
           ElementName = name,
           ElementClass = class,
           GeneSymbol = TargetGene,
           GeneEnsemblID = TargetGeneEnsembl_ID,
           GeneTSS = TargetGeneTSS,
           CellAnnotation,
           Score = E2G.Score.qnorm,
           distance,
           feature.normalizedATAC_prom = normalizedATAC_prom,
           feature.numTSSEnhGene = numTSSEnhGene,
           feature.numNearbyEnhancers = numNearbyEnhancers,
           feature.ubiqExpressed = ubiqExpressed,
           feature.numCandidateEnhGene = numCandidateEnhGene,
           feature.ABC.Score = ABC.Score,
           isSelfPromoter,
           normalizedATAC_enh)
} else if (grepl("ATAC", opt$method)) { # scATAC thresholded predictions -- unchanged narrow column set
    pred <- pred %>%
    mutate(ElementChr = chr,
          name = paste0(ElementChr, ":", start, "-", end)) %>%
    select(ElementChr,
          ElementStart = start,
          ElementEnd = end,
          ElementName = name,
          ElementClass = class,
          GeneSymbol = TargetGene,
          GeneEnsemblID = TargetGeneEnsembl_ID,
          GeneTSS = TargetGeneTSS,
          CellAnnotation,
          Score = E2G.Score.qnorm,
          isSelfPromoter = isSelfPromoter)
} else if (is.null(opt$threshold)) {
  # scE2G_multiome unthresholded ("full") predictions -- wider column set
  # (2026-07-23): every model-feature column gets a "feature." prefix, ending
  # in feature.ARC.E2G.Score (Multiome-only -- integrates ABC.Score and
  # Kendall). ABC.Score is NOT a feature for Multiome (unlike scATAC, where
  # it IS feature.-prefixed and sits last in the feature block in place of
  # ARC.E2G.Score) -- kept unprefixed. Kendall (also Multiome-only) likewise
  # stays unprefixed. ElementName keeps the raw "name" column's own
  # "class|chrN:start-end" value as-is. Score sits right after CellAnnotation
  # per the consortium-standard column order (2026-07-24) -- every remaining
  # column is "[Optional columns]", kept in the same relative order as
  # before, just moved after Score.
  pred <- pred %>%
    select(ElementChr = chr,
          ElementStart = start,
          ElementEnd = end,
          ElementName = name,
          ElementClass = class,
          GeneSymbol = TargetGene,
          GeneEnsemblID = TargetGeneEnsembl_ID,
          GeneTSS = TargetGeneTSS,
          CellAnnotation,
          Score = E2G.Score.qnorm,
          distance,
          feature.normalizedATAC_prom = normalizedATAC_prom,
          feature.numTSSEnhGene = numTSSEnhGene,
          feature.numNearbyEnhancers = numNearbyEnhancers,
          feature.ubiqExpressed = ubiqExpressed,
          feature.numCandidateEnhGene = numCandidateEnhGene,
          feature.ARC.E2G.Score = ARC.E2G.Score,
          isSelfPromoter,
          Kendall,
          normalizedATAC_enh,
          RNA_pseudobulkTPM,
          RNA_meanLogNorm,
          RNA_percentCellsDetected,
          ABC.Score)
} else { # scE2G_multiome thresholded predictions -- unchanged narrow column set
  pred <- pred %>%
    mutate(ElementChr = chr,
          name = paste0(ElementChr, ":", start, "-", end)) %>%
    select(ElementChr,
          ElementStart = start, ElementEnd = end, ElementName = name,
          ElementClass = class, GeneSymbol = TargetGene,
          GeneEnsemblID = TargetGeneEnsembl_ID, GeneTSS = TargetGeneTSS,
          CellAnnotation, Score = E2G.Score.qnorm,
          isSelfPromoter = isSelfPromoter)
}

# save to output file
message("Writing to output file...")
if (tools::file_ext(opt$output_file) == "gz") {
  
  # save to gzip compressed file
  tmp_file <- tools::file_path_sans_ext(opt$output_file)
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
