#'---
#' title: Counts Overview
#' author:  mumichae, salazar
#' wb:
#'  params:
#'   - ids: '`sm parser.outrider_ids`'
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input: 
#'   - summaries: '`sm expand(config["htmlOutputPath"] + "/AberrantExpression/Counting/{annotation}/Summary_{dataset}.html",
#'    annotation=list(config["geneAnnotation"].keys()), dataset=parser.outrider_ids)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "counting_overview.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/counting_overview.snakemake")

groups <- names(snakemake@params$ids)
gene_annotation_names <- names(snakemake@config$geneAnnotation)
summaries_titles <- paste(gene_annotation_names, groups)
summaries <- paste('[', summaries_titles ,'](', gsub(snakemake@config$htmlOutputPath, ".", snakemake@input$summaries), ')', sep = '')
summaries <- paste(summaries, sep = '\n')
#' Summaries:  `r summaries`
