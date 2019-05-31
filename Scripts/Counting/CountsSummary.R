#'---
#' title: Counts Summary
#' author: Daniela Andrade, Michaela Muller
#' wb:
#'  input: 
#'    - counts: '`sm parser.getProcDataDir() + "/{annotation}/counts/{dataset}/total_counts.Rds"`'
#'    - ods: '`sm parser.getProcResultsDir() + "/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  output:
#'   - wBhtml: "Output/html/Counting/{annotation}/CountingSummary_{dataset}.html"
#'  type: noindex
#'---


saveRDS(snakemake, "tmp/counting_summary.snakemake")

snakemake <- readRDS("../../tmp/counting_summary.snakemake")

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(SummarizedExperiment)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  #library(ggbeeswarm)
  library(ggthemes)
  library(cowplot)
  library(GenomicAlignments)
  library(dplyr)
  library(tidyr)
  #library(DASSIE)
  devtools::load_all("../genetic-diagnosis-tools")
})


#' ## Read counts object
counts_obj <- readRDS(snakemake@input$counts)
class(counts_obj)
dim(counts_obj)
colnames(counts_obj)

#' ## Read filtered counts object
ods <- readRDS(snakemake@input$ods)
class(ods)
dim(ods)
colnames(ods)

#' # Counting Results
#+ plotheatmap, fig.height=8, fig.width=8
#plotCountCorHeatmap(ods, main = "Raw Counts Correlation", normalized = F)
#plotCountCorHeatmap(ods, main = "Normalised Counts Correlation")

#+ globalqq
#plotQQ(ods, global = T)

#feature_dt <- data.table(gene_id = rownames(counts_obj), meanCounts = rowMeans(counts(counts_obj)))

#' # Count Statistics
#' ## Mean Count Distribution
#' The following plots are based on mean counts over all samples.  

# dist_hist <- ggplot(feature_dt, aes(x = meanCounts, y = stat(count))) +
#   geom_density() +
#   scale_x_log10() + 
#   labs(x = "Mean Counts", y = "Frequency") +
#   guides(col = guide_legend(title = NULL)) +
#   theme(legend.position = "bottom")
# 
# dist_box <- ggplot(feature_dt, aes(x = "", y = meanCounts)) + 
#   geom_boxplot() +
#   scale_y_log10() + 
#   geom_text(data = feature_dt[, .(median = round(median(meanCounts), 2))],
#             aes(x = "", y = median, label = median), 
#             size = 3, vjust = -1.5) +
#   labs(x = "", y = "Mean Counts") +
#   theme(legend.position = "none")

#+ meanCounts, fig.height=7, fig.width=14
#plot_grid(dist_hist, dist_box)
















