#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - counts: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/total_counts.Rds"`'
#'   - txdb: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/txdb.db"`'
#'  output:
#'   - ods: '`sm parser.getProcResultsDir() + "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  type: script
#'---

saveRDS(snakemake,  file.path(snakemake@params$tmpdir, "filter_counts.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/filter_counts.snakemake")

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(GenomicFeatures)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
})

counts <- readRDS(snakemake@input$counts)
ods <- OutriderDataSet(counts)
txdb <- loadDb(snakemake@input$txdb)

# filter not expressed genes
fpkmCutoff <- snakemake@config$aberrantExpression$fpkmCutoff
ods <- filterExpression(ods, gtfFile=txdb, filter=F, fpkmCutoff=fpkmCutoff,
                        addExpressedGenes=T)

# change row names from gene ID to gene name
if (snakemake@config$aberrantExpression$useGeneNames) {
  rownames(ods) <- rowData(ods)$gene_name
}

# add column for genes with at least 1 gene
rowData(ods)$counted1sample = rowSums(assay(ods)) > 0

# Save the ods before filtering to preserve the original number of genes
saveRDS(ods, snakemake@output$ods)
