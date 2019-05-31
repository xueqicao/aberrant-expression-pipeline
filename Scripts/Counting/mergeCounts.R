#'---
#' title: Merge the counts for all samples
#' author: Michaela Muller
#' wb:
#'  input: 
#'    - counts: '`sm lambda wildcards: parser.getCountFileByOutriderGroup(wildcards.annotation, wildcards.dataset)`'
#'    - gene_annot_dt: '`sm lambda wildcards: parser.getGeneInfoFile(wildcards.annotation) `'
#'  output:
#'    - counts: '`sm parser.getProcDataDir() + "/{annotation}/counts/{dataset}/total_counts.Rds"`'
#'  threads: 30
#'  type: script
#'---

saveRDS(snakemake, "tmp/count_all.snakemake")
# snakemake <- readRDS("tmp/count_all.snakemake")

suppressPackageStartupMessages({
    library(BiocParallel)
    library(SummarizedExperiment)
    library(data.table)
    library(dplyr)
})

register(MulticoreParam(snakemake@threads))

# Read counts
counts_list <- bplapply(snakemake@input$counts, readRDS)
names(counts_list) <- snakemake@config$outrider[[snakemake@wildcards$dataset]]
message(paste("read", length(counts_list), 'files'))

# merge counts
## more efficient merging matrices instead of complete SE
merged_assays <- do.call(cbind, bplapply(counts_list, assay, withDimnames = F))
total_counts <- SummarizedExperiment(assays = list(counts = merged_assays))
# total_counts <- do.call(cbind, counts_list)
colnames(total_counts) <- names(counts_list)
rownames(total_counts) <- rownames(counts_list[[1]])
rowRanges(total_counts) <- rowRanges(counts_list[[1]])

# Add gene annotation data (rowData)
gene_annot_dt <- fread(snakemake@input$gene_annot_dt)
row_data <- data.table(gene_id_unique = rownames(total_counts))
row_data <- left_join(row_data, gene_annot_dt[,.(gene_id_unique, gene_name_unique, gene_type, gene_status)], by = "gene_id_unique")
rowData(total_counts) <- row_data

# Add sample annotation data (colData)

col_data <- data.table(RNA_ID = colnames(total_counts))
# Removed adding sample_anno info. Do this later in the Analysis script when needed
# sample_anno <- fread(snakemake@config$SAMPLE_ANNOTATION)
# col_data <- left_join(col_data, unique(sample_anno[, .(RNA_ID, GENDER, BATCH, TISSUE, GROWTH_MEDIUM)]), by = snakemake@config$rna_assay)

setnames(col_data, 'RNA_ID', 'sampleID')
colData(total_counts) <- DataFrame(col_data, row.names = col_data$sampleID)

saveRDS(total_counts, snakemake@output$counts)
