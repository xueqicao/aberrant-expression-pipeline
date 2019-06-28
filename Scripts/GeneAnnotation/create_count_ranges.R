#'---
#' title: Create count ranges
#' author: mumichae
#' wb:
#'  input:
#'   - txdb: '`sm parser.getProcDataDir() + "/{annotation}/txdb.db"`'
#'   - gene_name_mapping: '`sm parser.getProcDataDir() + "/{annotation}/gene_name_mapping.Rds"`'
#'  output:
#'   - count_object: '`sm parser.getProcDataDir() + "/{annotation}/count_ranges.Rds"`'
#'  type: script
#'---


saveRDS(snakemake,  "tmp/count_ranges.snakemake")
# snakemake <- readRDS("tmp/count_ranges.snakemake")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicFeatures)
  library(GenomicRanges)
})

txdb <- loadDb(snakemake@input$txdb)
txdb <- keepStandardChromosomes(txdb)
count_object <- exonsBy(txdb, by = "gene")

saveRDS(count_object, snakemake@output$count_object)
