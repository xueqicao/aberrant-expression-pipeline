#'---
#' title: "OUTRIDER Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: mumichae, vyepez
#' wb:
#'  input:
#'   - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'   - results: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv"`'
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

# used for most plots
dataset_title <- paste("Dataset:", snakemake@wildcards$dataset)

#' ## Read the ods object
ods <- readRDS(snakemake@input$ods)
#' Number of samples: `r ncol(ods)`  
#' Number of genes: `r nrow(ods)`  

#'
#' ## Visualize
#' ### Parameters
plotEncDimSearch(ods) +
    labs(title = dataset_title) +
    theme_cowplot() +
    background_grid() +
    scale_color_brewer(palette = "Set1")


#' ### Aberrant samples
plotAberrantPerSample(ods, main = dataset_title)


#' ### Batch correction
#+ countCorHeatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, normalized = FALSE, 
                    main = paste0('Raw Counts (', dataset_title, ')'))
plotCountCorHeatmap(ods, normalized = TRUE, 
                    main = paste0('Normalized Counts (', dataset_title, ')'))


#' ### Expression by gene per sample
#+ geneSampleHeatmap, fig.height=15, fig.width=6
plotCountGeneSampleHeatmap(ods, normalized = FALSE, nGenes = 50,
                           main = paste0('Raw Counts (', dataset_title, ')'),
                           bcvQuantile = .95, show_names = 'row')
plotCountGeneSampleHeatmap(ods, normalized = TRUE, nGenes = 50,
                           main = paste0('Normalized Counts (',dataset_title,')'),
                           bcvQuantile = .95, show_names = 'row')


#' ### BCV - Biological Cofficient of Variation
# function to calculate BCV before autoencoder
estimateThetaWithoutAutoCorrect <- function(ods){
  
  ods1 <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
  # use rowMeans as expected means
  normalizationFactors(ods1) <- matrix(rowMeans(counts(ods1)), ncol=ncol(ods1), nrow=nrow(ods1))
  ods1 <- fit(ods1)
  theta(ods1)
  
  return(theta(ods1))
}

before <- data.table(when = "Before",
                     BCV = 1/sqrt(estimateThetaWithoutAutoCorrect(ods)))
after <- data.table(when = "After", BCV = 1/sqrt( theta( ods )))
bcv_dt <- rbind(before, after)

# boxplot of BCV Before and After Autoencoder
#+ BCV, fig.height=5, fig.width=6
ggplot(bcv_dt, aes(when, BCV)) +
    geom_boxplot() +
    theme_clean(base_size = 16) +
    labs(x = "",
         title = paste0("BCV - Before and After Autoencoder (", 
                        dataset_title, ")"))


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
web_dir <- snakemake@config$webDir

if (!is.null(results_link)) {
    results_link <- paste0(web_dir, "/aberrant_expression/results/{annotation}/outrider/{dataset}/OUTRIDER_results.tsv")
else {
    results_link <- snakemake@input$results
}
#' [Download OUTRIDER results table](`r results_link`)

