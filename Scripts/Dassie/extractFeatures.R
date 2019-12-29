#'---
#' title: Extract regions of specified length upstream and downstream of TSS and PAS
#' author: Michaela Müller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - tssWindow: '`sm config["aberrantExpression"]["dassie"]["tssWindow"]`'
#'   - pasWindow: '`sm config["aberrantExpression"]["dassie"]["pasWindow"]`'
#'  input:
#'   - txdb: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/txdb.db"`'
#'  output:
#'   - featureRegions: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/dassie/count_ranges.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "dassie_features.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/dassie_features.snakemake")

suppressPackageStartupMessages({
  library(tidyr)
  library(GenomicRanges)
  library(GenomicFeatures)
  #library(DASSIE)
  devtools::load_all("DASSIE")
})

txdb <- loadDb(snakemake@output$txdb)

#' TODO: exon mapping
#' 1. get exons by transcript
#' 2. map exons to tx to gene
#' 3. remove exons that are not filtered
exon_by_tx_dt <- exon_to_tx(txdb)
setnames(exon_by_tx_dt, "transcript_name", "transcript_id")
gencode_exons_dt <- merge(exon_by_tx_dt, gene_mapping, by = "transcript_id")
gencode_exons_dt[, transcriptID := transcript_id]
exon_annotation

# TODO: select most representative transcript
# Take longest transcript per gene
exon_annotation[, exons_per_tx := .N, by = transcript_id]
exon_annotation[, tx_length := sum(width), by = transcript_id]
setorder(exon_annotation, -exons_per_tx)
setorder(exon_annotation, -tx_length)
longest_isoform <- exon_annotation[, .SD[1], by = gene_id]

count_ranges <- readRDS(snakemake@input$count_ranges)
exon_annotation <- readRDS(snakemake@input$count_ranges)

# Extract all regions around TSS and PAS
tss <- tssRegions(exon_annotation, window = snakemake@params$TSS_WINDOW)
pas <- pasRegions(exon_annotation, window = snakemake@params$PAS_WINDOW)
feature_dt <- rbind(tss, pas)
feature_dt[, feature_id := 1:.N]
saveRDS(feature_dt, snakemake@output$allRegions)

feature_dt <- feature_dt[transcript_id %in% longest_isoform$transcript_id]
# Remove redundant regions
feature_dt <- feature_dt[, .SD[1], by = c("seqnames", "start", "end", "strand", "feature", "location")]

# Convert to GenomicRanges and save
feature_regions <- makeGRangesFromDataFrame(feature_dt, keep.extra.columns = T)
names(feature_regions) <- feature_regions$feature_id
saveRDS(feature_regions, snakemake@output$featureRegions)
