#'---
#' title: Counts Overview
#' author:  mumichae, salazar
#' wb:
#'  input: 
#'  - summaries: '`sm expand(config["htmlOutputPath"] + "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html", annotation=list(config["GENE_ANNOTATION"].keys()), dataset=parser.outrider_filtered)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/counting_overview.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantExpression/counting_overview.snakemake")

groups <- names(snakemake@config$outrider_all)
gene_annotation_names <- names(snakemake@config$GENE_ANNOTATION)
summaries_titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', summaries_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries:  `r summaries`
