#'---
#' title: Filter Counts for OUTRIDER
#' author: Michaela Mueller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - ods: '`sm parser.getProcResultsDir() + 
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods_unfitted.Rds"`'
#'  output:
#'   - ods: '`sm parser.getProcResultsDir() + 
#'           "/aberrant_expression/{annotation}/outrider/{dataset}/ods.Rds"`'
#'  type: script
#'  threads: 30
#'---

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(SummarizedExperiment)
    library(ggplot2)
    library(data.table)
    library(dplyr)
    library(magrittr)
})

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "outrider.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/outrider.snakemake")

ods <- readRDS(snakemake@input$ods)
implementation <- snakemake@config$aberrantExpression$implementation
register(MulticoreParam(snakemake@threads))

## subset filtered and estimate
ods <- ods[mcols(ods)$passedFilter,] 
ods <- estimateSizeFactors(ods)

## find optimal encoding dimension
a <- 5 
b <- min(ncol(ods), nrow(ods)) / 3   # N/3
Nsteps <- min(20, ncol(ods)/3, nrow(ods)/3)   # Do at most 20 steps or N/3
# Do unique in case 2 were repeated
pars_q <- round(exp(seq(log(a),log(b),length.out = Nsteps))) %>% unique
ods <- findEncodingDim(ods, params = pars_q, lnorm = TRUE, implementation = implementation)

## fit OUTRIDER
ods <- OUTRIDER(ods, implementation = implementation)
message("outrider fitting finished")

# Save the new ods with a date stamp
op <- snakemake@output$ods
op_date <- paste0(strsplit(op, "\\.")[[1]][1], "-", format(Sys.time(), "%Y%m%d") , ".Rds")
saveRDS(ods, op_date)

# Create a link to the previous file
if(file.exists(op)) file.remove(op)
file.symlink(op_date, op)
