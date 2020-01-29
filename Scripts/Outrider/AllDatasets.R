#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  params:
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'    - summaries: '`sm expand(config["htmlOutputPath"] + "/AberrantExpression/Outrider/{annotation}/Summary_{dataset}.html",
#'    annotation=config["geneAnnotation"].keys() , dataset=config["aberrantExpression"]["groups"])`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "outrider_overview.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/outrider_overview.snakemake")


groups <- snakemake@config$aberrantExpression$groups 
gene_annotation_names <- names(snakemake@config$geneAnnotation)
summaries_titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', summaries_titles ,'](', 
                   gsub(snakemake@config$htmlOutputPath, ".", 
                        snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries: `r summaries`
