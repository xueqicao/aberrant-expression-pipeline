#'---
#' title: Create txdb from GTF
#' author: mumichae
#' wb:
#'  input:
#'   - gtf: '`sm lambda wildcards: parser.getGeneAnnotationFile(wildcards.annotation) `'
#'  output:
#'   - txdb: '`sm parser.getProcDataDir() + "/{annotation}/mytxdb.db"`'
#'  type: script
#'---

saveRDS(snakemake,  "tmp/txdb.snakemake")
# snakemake <- readRDS("tmp/txdb.snakemake")

suppressPackageStartupMessages({
  library(GenomicFeatures)
})

txdb <- makeTxDbFromGFF(snakemake@input$gtf)
seqlevelsStyle(txdb) <- snakemake@config$CHROMOSOME_FORMAT_rna
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, snakemake@output$txdb, overwrite=TRUE)