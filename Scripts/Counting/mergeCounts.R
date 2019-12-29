#'---
#' title: Merge the counts for all samples
#' author: Michaela Müller
#' wb:
#'  py:
#'   - |
#'    def getCountFiles(annotation, feature_type, dataset):
#'        ids = parser.outrider_ids[dataset]
#'        file_stump = parser.getProcDataDir() + f"/aberrant_expression/{annotation}/{feature_type}/counts/"
#'        return expand(file_stump + "{sampleID}.Rds", sampleID=ids)
#'  params:
#'    - ids: '`sm parser.outrider_ids`'
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input: 
#'    - counts: '`sm lambda wildcards: getCountFiles(wildcards.annotation, wildcards.feature_type, wildcards.dataset)`'
#'  output:
#'    - counts: '`sm parser.getProcDataDir() +
#'               "/aberrant_expression/{annotation}/{feature_type}/{dataset}/total_counts.Rds"`'
#'  threads: 30
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "merge_counts.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/merge_counts.snakemake")

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(BiocParallel)
    library(SummarizedExperiment)
})

register(MulticoreParam(snakemake@threads))

# Read counts
counts_list <- bplapply(snakemake@input$counts, readRDS)
names <- snakemake@params$ids[[snakemake@wildcards$dataset]]
names(counts_list) <- names

message(paste("read", length(counts_list), 'files'))

# merge counts
## more efficient merging matrices instead of complete SE
merged_assays <- do.call(cbind, lapply(counts_list, assay, withDimnames=FALSE))
total_counts <- SummarizedExperiment(assays=list(counts=merged_assays))

# total_counts <- do.call(cbind, counts_list)
colnames(total_counts) <- names(counts_list)
rownames(total_counts) <- rownames(counts_list[[1]])
rowRanges(total_counts) <- rowRanges(counts_list[[1]])

# Add sample annotation data (colData)
sample_anno <- fread(snakemake@config$sampleAnnotation)
col_data <- data.table(RNA_ID = colnames(total_counts))
col_data <- left_join(col_data, sample_anno, by = "RNA_ID")
rownames(col_data) <- col_data$RNA_ID
colData(total_counts) <- as(col_data, "DataFrame")

saveRDS(total_counts, snakemake@output$counts)
