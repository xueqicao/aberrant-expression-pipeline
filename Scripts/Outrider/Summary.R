#'---
#' title: "OUTRIDER Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - results: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'   - results_public: '`sm config["webDir"] + "/aberrant_expression/results/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html"`'
#'  type: noindex
#'---

#+ echo=F
saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/outrider_summary.snakemake") )
# snakemake <- readRDS(snakemake@config$tmpdir, "/AberrantExpression/outrider_summary.snakemake"))

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


#' ## Read the ods object
ods <- readRDS(snakemake@input$ods)
# Number of samples and genes
dim(ods)


#' ## Visualize
#' ### Parameters
plotEncDimSearch(ods) +
  theme_cowplot() +
  background_grid() +
  scale_color_brewer(palette = "Set1")


#' ### Aberrant samples
plotAberrantPerSample(ods, main = snakemake@wildcards$dataset)


#' ### Batch correction
#+ heatmap1, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, normalized = FALSE, main = paste('Raw Counts (', snakemake@wildcards$dataset, ')'))
plotCountCorHeatmap(ods, normalized = TRUE, main = paste('Normalized Counts (', snakemake@wildcards$dataset, ')'))


#' ### gene sample heatmap
# subset the genes with top BCV for gene sample heatmap 
bcv = 1/sqrt( theta( ods ))
bcv_sub = bcv > quantile( bcv, probs = c(0.95))
ods_sub = ods[ bcv_sub, ]

# #+ heatmap2, fig.height=15, fig.width=6
# plotCountGeneSampleHeatmap(ods_sub, normalized = FALSE, nGenes = nrow(ods_sub),
#                           main = paste('Raw Counts (', snakemake@wildcards$dataset, ')'))
#plotCountGeneSampleHeatmap(ods_sub, normalized = TRUE, nGenes = nrow(ods_sub),
#                           main = paste('Normalized Counts (', snakemake@wildcards$dataset, ')'))


#' ### BCV - Biological Cofficient of Variation
# function to calculate BCV before autoencoder
estimateThetaWithoutAutoCorrect = function(ods){
  
  ods1 <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
  # use rowMeans as expected means
  normalizationFactors(ods1) <- matrix(rowMeans(counts(ods1)), ncol=ncol(ods1), nrow=nrow(ods1))
  ods1 <- fit(ods1)
  theta(ods1)
  
  return(theta(ods1))
}

# boxplot of BCV Before and After Autoencoder
boxplot( 1/sqrt( estimateThetaWithoutAutoCorrect( ods )) ,1/sqrt( theta( ods )), 
         names = c("Before","After"),
         main = paste("BCV - Before and After Autoencoder (",snakemake@wildcards$dataset, ')'),
         ylab = "BCV")


#' ## Results
res <- fread(snakemake@input$results)

#' ### How many samples with at least one gene
res[, uniqueN(sampleID)]

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


#' ## Results table
#+echo=F
res[, pValue := format(pValue, scientific = T, digits = 2)]
res[, padjust := format(padjust, scientific = T, digits = 2)]
DT::datatable(res, caption = "OUTRIDER results", style = 'bootstrap', filter = 'top')

#' ### Download Aberrant Samples Table
results_link <- snakemake@input$results_public
#' [Download OUTRIDER results table](`r results_link`)
