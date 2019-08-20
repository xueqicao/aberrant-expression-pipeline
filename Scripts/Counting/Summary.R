#'---
#' title: Counts Summary
#' author: Daniela Andrade, Michaela Muller
#' wb:
#'  input: 
#'    - counts: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/counts/{dataset}/total_counts.Rds"`'
#'    - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html"`'
#'  type: noindex
#'---


saveRDS(snakemake, "tmp/counting_summary.snakemake")
#snakemake <- readRDS("tmp/counting_summary.snakemake")

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
  devtools::load_all("../genetic-diagnosis-tools")
})


counts_obj <- readRDS(snakemake@input$counts)
ods <- readRDS(snakemake@input$ods)

#' # Counting Results
#+ plotheatmap, fig.height=8, fig.width=8
plotCountCorHeatmap(ods, main = "Raw Counts Correlation", normalized = F)

#' ## Filtering
# TODO: number reads counted, histogram(colSums(ods)) (before and after filtering)

# # of covered genes per sample before and after filtering. 
# Define covered genes using 2 cutoffs: FPKM>1 and raw_counts > 10. 
# Then make a histogram (before and after filtering) overlapping the distributions.
hist(colSums(fpkm(ods)>1))
hist(colSums(counts(ods, normalized = F)>1))
hist(colSums(counts(ods, normalized = F)>10))

#' # Count Statistics
feature_dt <- data.table(gene_id = rownames(counts_obj), meanCounts = rowMeans(assay(counts_obj)))
#' ## Mean Count Distribution
#' The following plots are based on mean counts over all samples.  

dist_hist <- ggplot(feature_dt, aes(x = meanCounts, y = stat(count))) +
  geom_histogram() + #geom_density() +
  scale_x_log10() +
  labs(x = "Mean Counts", y = "Frequency") +
  guides(col = guide_legend(title = NULL)) +
  theme(legend.position = "bottom")

dist_box <- ggplot(feature_dt, aes(x = "", y = meanCounts)) +
  geom_boxplot() +
  scale_y_log10() +
  geom_text(data = feature_dt[, .(median = round(median(meanCounts), 2))],
            aes(x = "", y = median, label = median),
            size = 3, vjust = -1.5) +
  labs(x = "", y = "Mean Counts") +
  theme(legend.position = "none")

#+ meanCounts, fig.height=7, fig.width=14
plot_grid(dist_hist, dist_box)
















