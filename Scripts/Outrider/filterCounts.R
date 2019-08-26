#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - counts: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/counts/{dataset}/total_counts.Rds"`'
#'   - txdb: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/txdb.db"`'
#'  output:
#'   - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'   - plot: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/filtered_hist.png"`'
#'   - filtered_counts: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/counts/{dataset}/filtered_counts.Rds"`'
#'  type: script
#'---

saveRDS(snakemake,  paste0(snakemake@config$tmpdir, "/AberrantExpression/filter_counts.snakemake") )
# snakemake <- readRDS( paste0(snakemake@config$tmpdir, "/AberrantExpression/filter_counts.snakemake") )

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
})

counts <- readRDS(snakemake@input$counts)
ods <- OutriderDataSet(counts)
txdb <- loadDb(snakemake@input$txdb)

# filter not expressed genes ### Here we filter almost all samples --> zero counts
ods <- filterExpression(ods, gtfFile=txdb, filter=FALSE, fpkmCutoff=snakemake@config$fpkmCutoff)
g <- plotFPKM(ods) + theme_bw(base_size = 14)
ggsave(snakemake@output$plot, g)

# change row names from gene ID to gene name
rownames(ods) <- rowRanges(ods)$gene_name

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# Save the ods object before filtering, so as to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)

# Save the filtered count matrix (as a matrix)
ods <- ods[mcols(ods)$passedFilter,]
saveRDS(counts(ods), snakemake@output$filtered_counts)



