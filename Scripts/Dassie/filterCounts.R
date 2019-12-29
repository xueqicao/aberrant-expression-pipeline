#'---
#' title: Preprocess Counts for running OUTRIDER
#' author: Michaela Müller
#' wb:
#'  input:
#'   - filtered_genes: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/expression/{dataset}/ods_unfitted.Rds"`'
#'   - counts: '`sm parser.getProcDataDir() +
#'              "/aberrant_expression/{annotation}/dassie/{dataset}/total_counts.Rds"`'
#'  output:
#'   - ods: '`sm parser.getProcDataDir() +
#'           "/aberrant_expression/{annotation}/dassie/{dataset}/ods_unfitted.Rds"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "dassie_filter.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/dassie_filter.snakemake")

suppressPackageStartupMessages({
  library(OUTRIDER)
  library(DASSIE)
})

ods_gene <- readRDS(snakemake@input$ods_gene)
feature_annotation <- readRDS(snakemake@input$feature_annotation)
feature_counts <- readRDS(snakemake@input$feature_counts)

# Filter
# by zero counts
count_matrix <- assay(feature_counts)
# by_counts <- rowSums2(count_matrix == 0) > 0.75 * ncol(count_matrix)
by_counts <- rowSums2(count_matrix == 0) > 1

# by gene
feature_info <- featureInfo(feature_counts, summary = T)
feature_info <- merge(feature_info, unique(feature_annotation[, .(feature_id, gene_id)]), by = "feature_id")
genewise_filtered <- feature_info[gene_name %in% rownames(ods_gene) 
                                  | gene_id %in% rownames(ods_gene)]

# apply filters
ods <- OutriderDataSet(feature_counts)
filter <- rownames(ods) %in% genewise_filtered$feature_id
filter[by_counts] <- F
mcols(ods)$passedFilter <- filter
ods <- ods[mcols(ods)$passedFilter == T]

saveRDS(ods, snakemake@output$ods)

# TODO: put into summary HTML
# filtering info
parameters <- c("All extracted genes",
                "All extracted features",
                "Features with zero counts",
                "FPKM selected genes (incl. biomart & singletons)",
                "Genes after genewise FPKM selection",
                "Features after genewise FPKM selection",
                "Final number of genes",
                "Final number of features")
values <- c(uniqueN(feature_info$gene_name), # all genes
            nrow(feature_info), # all features
            sum(by_counts), # all zero counts
            nrow(ods_gene), # all fpkm selected
            uniqueN(genewise_filtered$gene_name), # genes fpkm w/o singletons (only gencode 19)
            nrow(genewise_filtered), # features after fpkm w/o singletons (only gencode 19)
            uniqueN(mcols(ods)$gene_name), # final genes
            nrow(ods)) # final features
