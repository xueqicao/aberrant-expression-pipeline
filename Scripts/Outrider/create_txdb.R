#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - gtf: '`sm lambda wildcards: parser.getGeneAnnotationFile(wildcards.annotation) `'
#'  output:
#'   - txdb: '`sm parser.getProcDataDir() + "/{annotation}/txdb.db"`'
#'  type: script
#'---

saveRDS(snakemake,  "tmp/txdb.snakemake")
# snakemake <- readRDS("tmp/txdb.snakemake")

suppressPackageStartupMessages({
  library(GenomicFeatures)
})

gencode_txdb <- makeTxDbFromGFF(snakemake@input$gtf)
seqlevelsStyle(gencode_txdb) <- snakemake@config$CHROMOSOME_FORMAT
txdb <- keepStandardChromosomes(txdb)

saveDb(txdb, snakemake@output$txdb)