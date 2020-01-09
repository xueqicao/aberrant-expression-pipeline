#'---
#' title: Preprocess Counts for running OUTRIDER
#' author: Michaela Müller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - filtered_genes: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/expression/{dataset}/ods_unfitted.Rds"`'
#'   - counts: '`sm parser.getProcDataDir() +
#'              "/aberrant_expression/{annotation}/dassie/{dataset}/total_counts.Rds"`'
#'  output:
#'   - ods: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/dassie/{dataset}/ods_unfitted.Rds"`'
#'   - filter_stats: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/dassie/{dataset}/filter_stats.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "dassie_filter.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/dassie_filter.snakemake")

source("Scripts/_helpers/load_packages.R")
suppressPackageStartupMessages({
  library(OUTRIDER)
})

ods_gene <- readRDS(snakemake@input$filtered_genes)
feature_counts <- readRDS(snakemake@input$counts)
feature_info <- featureInfo(feature_counts, summary = TRUE)

# by zero counts
count_matrix <- assay(feature_counts)
feature_info[, zero_count_filtered := rowSums2(count_matrix == 0) == 0]
# by_counts <- rowSums2(count_matrix == 0) > 0.75 * ncol(count_matrix)

# by gene
filtered_genes <- rownames(ods_gene)[rowData(ods_gene)$passedFilter == TRUE]
feature_info[, gene_filtered := FALSE]
feature_info[gene_id %in% filtered_genes, gene_filtered := TRUE]

# create Outrider Dataset
ods <- OutriderDataSet(feature_counts)

# add filters
filtered <- feature_info[, zero_count_filtered == TRUE & gene_filtered == TRUE]
mcols(ods)$passedFilter <- filtered

saveRDS(ods, snakemake@output$ods)

# filtering info
parameters <- c("All genes",
                "All features",
                "Features with all zero counts",
                "Features of FPKM-filtered genes",
                "Remaining genes",
                "Features passing all filters")

values <- c(nrow(ods_gene), # all genes
            nrow(feature_info), # all features
            uniqueN(feature_info[zero_count_filtered == FALSE, feature_id]), 
            uniqueN(feature_info[gene_filtered == TRUE, feature_id]), 
            uniqueN(feature_info[gene_filtered == TRUE, gene_id]), 
            sum(filtered))

filter_stats <- data.table(description = parameters, value = values)

fwrite(filter_stats, snakemake@output$filter_stats)
