#'---
#' title: Results Overview
#' author: mumichae
#' wb:
#'  params:
#'    - ids: '`sm parser.outrider_ids`'
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'    - expression: '`sm expand(config["htmlOutputPath"] + 
#'                  "/AberrantExpression/Outrider/{annotation}/Expression_Summary_{dataset}.html",
#'                  annotation=config["geneAnnotation"].keys() , dataset=parser.outrider_ids)`'
#'    - dassie: '`sm expand(config["htmlOutputPath"] + 
#'               "/AberrantExpression/Outrider/{annotation}/DASSIE_Summary_{dataset}.html",
#'               annotation=config["geneAnnotation"].keys() , dataset=parser.outrider_ids)`'
#' output:
#'   html_document:
#'    code_folding: hide
#'    code_download: TRUE
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "outrider_overview.snakemake")  )
# snakemake <- readRDS(".drop/tmp/AE/outrider_overview.snakemake")

groups <- names(snakemake@params$ids)
gene_annotation_names <- names(snakemake@config$geneAnnotation)
summaries_titles <- paste(gene_annotation_names, groups)

markdown_links <- function(html_files, summaries_titles) {
  rel_path <- gsub(snakemake@config$htmlOutputPath, ".", html_files)
  summaries <- paste('[', summaries_titles ,'](', rel_path, ')', sep = '')
  summaries <- paste(summaries, sep = '\n')
}

#' ## OUTRIDER Summaries 
#' 
#' `r markdown_links(snakemake@input$expression)`
#'
#' ## DASSIE Summaries 
#' 
#' `r markdown_links(snakemake@input$dassie)`
#' 
