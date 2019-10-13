#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  params:
#'    - ids: '`sm parser.outrider_ids`'
#'  input:
#'    - summaries: '`sm expand(config["htmlOutputPath"] + "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html",
#'    annotation=config["geneAnnotation"].keys() , dataset=parser.outrider_ids)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/outrider_overview.snakemake")  )
# snakemake <- readRDS( paste0(snakemake@config$tmpdir, "/AberrantExpression/outrider_overview.snakemake") )


groups <- names(snakemake@params$ids)
gene_annotation_names <- names(snakemake@config$geneAnnotation)
summaries_titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', summaries_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
