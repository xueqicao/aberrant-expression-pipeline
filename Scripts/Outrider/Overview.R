#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  input:
#'  - summaries: '`sm expand(config["htmlOutputPath"] + "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html", annotation=list(config["GENE_ANNOTATION"].keys()) , dataset=parser.outrider_filtered)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, "tmp/outrider_overview.snakemake")
# snakemake <- readRDS("tmp/outrider_overview.snakemake")

groups <- names(snakemake@config$outrider_filtered)
gene_annotation_names <- names(snakemake@config$GENE_ANNOTATION)
summaries_titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', summaries_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
