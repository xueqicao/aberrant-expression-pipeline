#'---
#' title: OUTRIDER Results
#' author: mumichae
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'  output:
#'   - results: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_all: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results_all.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "outrider_results.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/outrider_results.snakemake")

suppressPackageStartupMessages({
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(SummarizedExperiment)
    library(OUTRIDER)
})

ods <- readRDS(snakemake@input$ods)
res <- results(ods, all = TRUE)

# Add fold change
res[, foldChange := round(2^l2fc, 2)]

# Save all the results and significant ones
saveRDS(res, snakemake@output$results_all)

# Subset to significant results
res <- res[padjust <= snakemake@config$aberrantExpression$padjCutoff & 
               abs(zScore) > snakemake@config$aberrantExpression$zScoreCutoff]

# Save results 
fwrite(res, snakemake@output$results, sep = "\t", quote = F)

web_dir <- snakemake@config$webDir
if (!is.null(web_dir)) {
    pub_res <- paste0(web_dir, 
                      "/aberrant_expression/results/",{snakemake@wildcards$annotation},"/outrider/",
                      {snakemake@wildcards$dataset},"/OUTRIDER_results.tsv")
    fwrite(res, pub_res, sep = "\t", quote = F)
}
