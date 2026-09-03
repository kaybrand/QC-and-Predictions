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
              help = "Include all columns in element/gene lists"),
  # Every other output mode is chosen by sniffing opt$input_file for
  # "gene_list"/"element_list". That cannot work here: this mode's input IS the
  # predictions file, whose name contains neither, so it would fall through to
  # the predictions branches and silently emit the wrong shape. Hence an
  # explicit flag rather than another filename pattern.
  make_option(c("--genes_from_predictions"), action = "store_true", default = FALSE,
              help = paste("Build the gene list from a scE2G_predictions file instead of a",
                           "scE2G_gene_list file. Yields only genes actually scored (those with",
                           "considered_for_predictions=TRUE), not the full annotation. Requires",
                           "--strand_annotation, because no scE2G_predictions file carries strand")),
  make_option(c("--strand_annotation"), type = "character", default = NULL,
              help = paste("processed_genes_file.bed, supplying GeneStrand for",
                           "--genes_from_predictions. Columns: chr/start/end/symbol/score/strand/",
                           "EnsemblID/gene_type"), metavar = "character")
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
} else if (opt$genes_from_predictions) {
  # Gene list built from scE2G_predictions rather than scE2G_gene_list.
  #
  # WHY: the predictions file is the authoritative statement of which genes this
  # prediction set actually scored. Verified on catlas Bgl_51 (scATAC): the
  # predictions carry 14942 distinct genes, exactly equal to the
  # considered_for_predictions=TRUE subset of scE2G_gene_list's 20531 -- so this
  # mode drops the 5589 genes that were annotated but never scored.
  #
  # WHY A SECOND INPUT: no scE2G_predictions file carries strand, in EITHER model
  # family -- multiome's 27 columns are a strict superset of scATAC's 20 and
  # neither includes it (checked 2026-09-02). Strand is a property of the gene
  # annotation, not of a prediction, so it is joined from processed_genes_file.bed.
  # That file is a pure genome annotation and was byte-identical (md5
  # e00a3a154a2e...) across all 157 cluster directories checked, so it contributes
  # nothing cluster-specific.
  if (is.null(opt$strand_annotation)) {
    stop("--genes_from_predictions requires --strand_annotation ",
         "(no scE2G_predictions file carries strand)", call. = FALSE)
  }
  genes <- pred %>%
    select(chr, GeneSymbol = TargetGene, GeneEnsemblID = TargetGeneEnsembl_ID,
           TSS = TargetGeneTSS) %>%
    distinct()
  # One TSS per gene is an assumption worth enforcing rather than silently
  # collapsing: 0 genes had >1 distinct (chr, TargetGeneTSS) when verified, so a
  # violation means scE2G changed something upstream.
  dup <- genes %>% count(GeneEnsemblID, GeneSymbol) %>% filter(n > 1)
  if (nrow(dup) > 0) {
    stop(nrow(dup), " gene(s) have more than one distinct (chr, TargetGeneTSS) in ",
         opt$input_file, "; first: ", dup$GeneSymbol[1], call. = FALSE)
  }
  ann <- fread(opt$strand_annotation,
               col.names = c("annChr", "annStart", "annEnd", "GeneSymbol",
                             "annScore", "GeneStrand", "GeneEnsemblID", "annType")) %>%
    select(GeneSymbol, GeneEnsemblID, GeneStrand, annChr, annStart, annEnd) %>%
    distinct() %>%
    mutate(annRow = row_number())
  joined <- genes %>% inner_join(ann, by = c("GeneSymbol", "GeneEnsemblID"))
  # An inner join that loses rows would quietly ship an incomplete gene list.
  # The join resolved 14942/14942 with 0 missing and 0 ambiguous strand when
  # verified, so any shortfall is a broken assumption, not an expected case.
  if (nrow(joined) != nrow(genes)) {
    missing <- anti_join(genes, ann, by = c("GeneSymbol", "GeneEnsemblID"))
    stop(nrow(genes) - nrow(joined), " of ", nrow(genes),
         " predicted gene(s) have no strand in ", opt$strand_annotation,
         "; first: ", missing$GeneSymbol[1], " (", missing$GeneEnsemblID[1], ")", call. = FALSE)
  }
  # The annotation's own TSS (strand-aware gene bound) must agree with the
  # predictions' TargetGeneTSS -- they did for all 14942. Disagreement means the
  # two inputs describe different gene models and the output would be incoherent.
  bad <- joined %>%
    mutate(annTSS = ifelse(GeneStrand == "+", annStart, annEnd)) %>%
    filter(chr != annChr | annTSS != TSS)
  if (nrow(bad) > 0) {
    stop(nrow(bad), " gene(s) disagree on TSS/chr between ", opt$input_file,
         " and ", opt$strand_annotation, "; first: ", bad$GeneSymbol[1], call. = FALSE)
  }
  # Ordered by the annotation's own row order, so the output is positionally
  # sorted and byte-stable across runs -- predictions row order is not.
  pred <- joined %>%
    arrange(annRow) %>%
    mutate(TSSStart = TSS - 250, TSSEnd = TSS + 250) %>%
    select(TSSChr = chr, TSSStart, TSSEnd, GeneSymbol, GeneEnsemblID, GeneStrand)
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
  # normalizedATAC_enh is NOT selected here (unlike Multiome's equivalent
  # branch below, which does have it): scATAC_powerlaw_v3's own predictions
  # table only ever has normalizedATAC_prom -- confirmed no
  # normalizedATAC_enh column exists in real scATAC_powerlaw_v3 output
  # (2026-08-13) -- enhancer-level ATAC normalization isn't a feature this
  # model computes, unlike Multiome's DHS-based features.
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
           isSelfPromoter)
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
