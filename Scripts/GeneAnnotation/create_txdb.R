#'---
#' title: Create txdb from GTF
#' author: mumichae
#' wb:
#'  input:
#'   - gtf: '`sm lambda wildcards: parser.getGeneAnnotationFile(wildcards.annotation) `'
#'  output:
#'   - txdb: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/txdb.db"`'
#'  type: script
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/txdb.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantExpression/txdb.snakemake") )

suppressPackageStartupMessages({
  library(GenomicFeatures)
})

txdb <- makeTxDbFromGFF(snakemake@input$gtf)
seqlevelsStyle(txdb) <- snakemake@config$CHROMOSOME_FORMAT_rna
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, snakemake@output$txdb, overwrite=TRUE)