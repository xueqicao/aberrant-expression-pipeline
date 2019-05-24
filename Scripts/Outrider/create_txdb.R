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

txdb <- makeTxDbFromGFF(snakemake@input$gtf)
seqlevelsStyle(txdb) <- snakemake@config$CHROMOSOME_FORMAT
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, snakemake@output$txdb, overwrite=TRUE)