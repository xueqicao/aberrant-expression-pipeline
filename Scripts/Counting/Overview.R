#'---
#' title: Counts Overview
#' author: Michaela Muller
#' wb:
#'  input: 
#'  - summaries: '`sm expand(config["htmlOutputPath"] + "/Counting/{annotation}/Summary_{dataset}.html", annotation=config["GENE_ANNOTATION"].keys(), dataset=parser.outrider_filtered)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, "tmp/counting_overview.snakemake")
# snakemake <- readRDS("tmp/counting_overview.snakemake")

groups <- names(snakemake@config$outrider_filtered)
gene_annotation_names <- names(gene_annotation_names, snakemake@config$GENE_ANNOTATION)
summaries_titles <- paste( groups)
summaries <- paste('[', summaries_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
