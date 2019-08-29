#'---
#' title: Count quality control values for all samples
#' author: Michaela Mueller
#' wb:
#'  input:
#'    - bam_files: '`sm lambda wildcards: parser.getFilePaths(wildcards.dataset, isRNA=True) `'
#'    - count_ranges: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/count_ranges.Rds" `'
#'  output:
#'    - QC: '`sm parser.getProcDataDir() + "/aberrant_expression/{annotation}/counts/{dataset}/qc.Rds"`'
#'  threads: 30
#'  type: script
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantExpression/qc.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantExpression/qc.snakemake"))

suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicAlignments)
})

count_ranges <- readRDS(snakemake@input$count_ranges)
for (file in snakemake@input$bam_files) {
  bam_file <- BamFile(file, yieldSize = 2e6)
  countBam(bam_file, 
           param=ScanBamParam(which = unlist(count_ranges)))
}
