#'---
#' title: Preprocess Gene Annotations
#' author: mumichae
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - gtf: '`sm lambda wildcards: parser.getGeneAnnotationFile(wildcards.annotation) `'
#'  output:
#'   - txdb: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/txdb.db"`'
#'   - gene_ranges: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/expression/count_ranges.Rds" `'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "annotation.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/annotation.snakemake")

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(data.table)
})

## Create txdb
txdb <- makeTxDbFromGFF(snakemake@input$gtf)
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, snakemake@output$txdb, overwrite=TRUE)

# save count ranges
gene_ranges <- exonsBy(txdb, by = "gene")
saveRDS(gene_ranges, snakemake@output$gene_ranges)
