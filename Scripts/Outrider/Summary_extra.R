#'---
#' title: OUTRIDER Summary
#' author: mumichae, vyepez, Xueqi
#' wb:
#'  input:
#'   - ods: '`sm parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - results: '`sm parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_public: '`sm config["webDir"] + "/results/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/Outrider/{annotation}/Summary_{dataset}.html"`'
#'  type: noindex
#'---

#+ echo=F
saveRDS(snakemake, "tmp/outrider_summary.snakemake")
# snakemake <- readRDS("tmp/outrider_summary.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(cowplot)
    library(data.table)
    library(dplyr)
    library(ggbeeswarm)
    library(ggthemes)
})




#' ## Read ods object
ods <- readRDS(snakemake@input$ods)

# Number of samples and genes
dim(ods)




#' ## Visualize
#' ### Parameters
barplot(sort(sizeFactors(ods)), main = paste('Size Factors (', snakemake@wildcards$dataset, ')'), xaxt = 'n', xlab = 'rank', ylab = 'Size Factors')
plotEncDimSearch(ods)

#' ### Aberrant samples
plotAberrantPerSample(ods, main = snakemake@wildcards$dataset)

#' ### Batch correction
#+ heatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, normalized = FALSE, main = paste('Raw Counts (', snakemake@wildcards$dataset, ')'))
plotCountCorHeatmap(ods, normalized = TRUE, main = paste('Normalized Counts (', snakemake@wildcards$dataset, ')'))

#+ heatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, normalized = FALSE, main = paste('Raw Counts (', snakemake@wildcards$dataset, ')'))
plotCountCorHeatmap(ods, normalized = TRUE, main = paste('Normalized Counts (', snakemake@wildcards$dataset, ')'))

#' ### BCV - Biological Cofficient of Variation
#' ods2 <- OutriderDataSet(countData = counts(ods))
# BCV before autoencoder
estimateThetaWithoutAutoCorrect = function(ods){

  ods1 <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
  # use rowMeans as expected means
  normalizationFactors(ods1) <- matrix(rowMeans(counts(ods1)), ncol=ncol(ods1), nrow=nrow(ods1))
  ods1 <- fit(ods1)
  theta(ods1)
  
  
  return(theta(ods1))
}
# boxplot
boxplot( 1/sqrt( estimateThetaWithoutAutoCorrect( ods )) ,1/sqrt( theta( ods )), 
         names = c("Before","After"),
         main = paste("BCV - Before and After Autoencoder (",snakemake@wildcards$dataset, ')'),
         ylab = "BCV")




#' ## Results
res <- fread(snakemake@input$results)

#' ### How many samples with at least one gene
res[, uniqueN(sampleID)]

#' ### Positive/Negative outlier
res[, zScore > 0]
res[, zScore < 0]

#' ### Aberrant samples
if (nrow(res) > 0) {
    ab_table <- res[AberrantBySample > nrow(ods)/1000, .N, by = .(sampleID)] %>% unique
    if (nrow(ab_table) > 0) {
      setorder(ab_table, N)
      DT::datatable(ab_table)
    } else {
      print("no aberrant samples")
    }
} else {
    print('no results')
}

#' ### Download Aberrant Samples Table
#' TODO
#results_link <- paste0('https://i12g-gagneurweb.informatik.tu-muenchen.de/project/genetic_diagnosis/results/', snakemake@wildcards$annotation,'/OUTRIDER_results_', snakemake@wildcards$dataset, '.tsv')
# #' [Download OUTRIDER results table](`r results_link`)
#DT::datatable(res, caption = "OUTRIDER results", style = 'bootstrap', filter = 'top')
