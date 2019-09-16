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

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/annotation.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantExpression/annotation.snakemake") )

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(rtracklayer)
  library(data.table)
  library(magrittr)
  library(tidyr)
})


## Create txdb
txdb <- makeTxDbFromGFF(snakemake@input$gtf)
txdb <- keepStandardChromosomes(txdb)

tmpFile <- tempfile()
saveDb(txdb, tmpFile)
R.utils::copyFile(tmpFile, snakemake@output$txdb, overwrite=TRUE)

# save count ranges
count_object <- exonsBy(txdb, by = "gene")
saveRDS(count_object, snakemake@output$count_object)


## Create Gene Name mapping
gtf_dt <- import(snakemake@input$gtf) %>% as.data.table
if (!"gene_name" %in% colnames(gtf_dt)) {
  gtf_dt[gene_name := gene_id]
}
gtf_dt <- gtf_dt[type == "gene", .(seqnames, start, end, strand, gene_id, gene_name, gene_type, gene_status)]

# make gene_names unique
gtf_dt[, N := 1:.N, by = gene_name] # warning message
gtf_dt[, gene_name_orig := gene_name]
gtf_dt[N > 1, gene_name := paste(gene_name, N, sep = '_')]
gtf_dt[, N := NULL]

fwrite(gtf_dt, snakemake@output$gene_name_mapping, na = NA)
