#'---
#' title: "Counts Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: Daniela Andrade, Michaela Mueller
#' wb:
#'  input: 
#'    - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html"`'
#'  type: noindex
#'---


saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/counting_summary.snakemake") )
#snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantExpression/counting_summary.snakemake") )

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(SummarizedExperiment)
  library(GenomicAlignments)
  library(ggplot2)
  #library(ggbeeswarm)
  library(ggthemes)
  library(cowplot)
  library(data.table)
  library(dplyr)
  library(tidyr)
})

ods <- readRDS(snakemake@input$ods)

#' # Count Quality Control
#' 
#' Compare number of records vs. read counts
#' 

#' # Filtering
quant <- .95
cnts_mtx <- counts(ods, normalized = F)
filter_mtx <- list(
  all = cnts_mtx,
  passed_FPKM = cnts_mtx[rowData(ods)$passedFilter,],
  min_1 = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 1, ],
  min_10 = cnts_mtx[rowQuantiles(cnts_mtx, probs = quant) > 10, ]
)
summary_dt <- lapply(names(filter_mtx), function(filter_name) {
  mtx <- filter_mtx[[filter_name]]
  data.table(gene_ID = rownames(mtx), median_counts = rowMeans(mtx), filter = filter_name)
}) %>% rbindlist
summary_dt[, filter := factor(filter, levels = c('all', 'passed_FPKM', 'min_1', 'min_10'))]

binwidth <- .2
p_hist <- ggplot(summary_dt, aes(x = median_counts, fill = filter)) +
  geom_histogram(binwidth = binwidth) +
  scale_x_log10() +
  facet_wrap(.~filter) +
  labs(x = "Mean counts per gene", y = "Frequency", title = 'Mean Count Distribution') +
  guides(col = guide_legend(title = NULL)) +
  scale_fill_brewer(palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "none")

p_dens <- ggplot(summary_dt, aes(x = median_counts, col = filter)) +
  geom_density(aes(y=binwidth * ..count..)) +
  scale_x_log10() +
  labs(x = "Mean counts per gene", y = "Frequency") +
  guides(col = guide_legend(title = NULL)) +
  scale_color_brewer(palette = "Paired") +
  theme_cowplot() +
  theme(legend.position = "top")

#+ meanCounts, fig.height=7, fig.width=14
plot_grid(p_hist, p_dens)

#' Expressed genes per sample
#+ expressedGenes, fig.height=7, fig.width=9
plotExpressedGenes(ods)
