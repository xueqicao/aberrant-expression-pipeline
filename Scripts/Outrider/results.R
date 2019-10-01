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

# Save all the results and significant ones
saveRDS(res, snakemake@output$results_all)

# Subset to significant results
res <- res[padjust <= snakemake@config$outrider_padjCutoff & 
               abs(zScore) > snakemake@config$outrider_zscoreCutoff]

# Save results 
fwrite(res, snakemake@output$results, sep = "\t", quote = F)
fwrite(res, snakemake@output$results_public, sep = "\t", quote = F)
