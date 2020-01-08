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
#'   - dassieRegions: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/dassie/count_ranges.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "dassie_features.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/dassie_features.snakemake")

source("Scripts/_helpers/load_packages.R")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
})

txdb <- loadDb(snakemake@input$txdb)

exon_dt <- getExons(txdb, filter_length = "tx_length", filter_exon = TRUE)

# Extract all regions around TSS and PAS
dassie_features <- getDassieRegions(exon_dt, tx_id = "tx_id",
                                    tssWindow = snakemake@params$tssWindow,
                                    pasWindow = snakemake@params$pasWindow)
saveRDS(dassie_features, snakemake@output$dassieRegions)

# TODO: features for quality control
