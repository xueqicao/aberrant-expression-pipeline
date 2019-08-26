#'---
#' title: OUTRIDER Results
#' author: mumichae
#' wb:
#'  input:
#'   - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'  output:
#'   - results: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_all: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results_all.Rds"`'
#'   - results_public: '`sm config["webDir"] + "/aberrant_expression/results/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'  type: script
#'---

#+ echo=F
saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/outrider_results.snakemake"))
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantExpression/outrider_results.snakemake"))

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

ods <- readRDS(snakemake@input$ods)
res <- OUTRIDER::results(ods, all = TRUE)
res[, FC := round(2^l2fc, 2)]
res[, geneID := toupper(geneID)]
#res <- add_gene_type(res, gene_name_col = 'geneID')
#saveRDS(res[,.(geneID, sampleID, pValue, padjust, zScore, l2fc, rawcounts, normcounts, meanCorrected, theta, aberrant, AberrantBySample, AberrantByGene, padj_rank, FC, gene_type)], snakemake@output$results_all)
saveRDS(res[,.(geneID, sampleID, pValue, padjust, zScore, l2fc, rawcounts, normcounts, meanCorrected, theta, aberrant, AberrantBySample, AberrantByGene, padj_rank, FC)], snakemake@output$results_all)

# Subset to significant results
res <- res[padjust <= .05]
#res <- add_all_gene_info(res, gene_name_col = 'geneID', dis_genes = F, gene_type = F)  # gene_type already added before

# Save results 
fwrite(res, snakemake@output$results, sep = "\t", quote = F)
fwrite(res, snakemake@output$results_public, sep = "\t", quote = F)
