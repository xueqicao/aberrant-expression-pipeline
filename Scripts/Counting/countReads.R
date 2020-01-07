#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - strand: '`sm lambda wildcards: hp.getCountSettings(parser,
#'               wildcards.sampleID, wildcards.feature_type, "strand")`'
#'   - count_mode: '`sm lambda wildcards: hp.getCountSettings(parser,
#'                   wildcards.sampleID, wildcards.feature_type, "count_mode")`'
#'   - paired_end: '`sm lambda wildcards: hp.getCountSettings(parser,
#'                    wildcards.sampleID, wildcards.feature_type, "paired_end")`'
#'   - overlap: '`sm lambda wildcards: hp.getCountSettings(parser,
#'                wildcards.sampleID, wildcards.feature_type, "overlap")`'
#'  input:
#'   - sample_bam: '`sm lambda wildcards: parser.getFilePath(wildcards.sampleID, 
#'                                                  file_type="RNA_BAM_FILE") `'
#'   - count_ranges: '`sm parser.getProcDataDir() + 
#'        "/aberrant_expression/{annotation}/{feature_type}/count_ranges.Rds" `'
#'  output:
#'   - counts: '`sm parser.getProcDataDir() + 
#'    "/aberrant_expression/{annotation}/{feature_type}/counts/{sampleID}.Rds"`'
#'  type: script
#'  threads: 1
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "counts.snakemake"))
# snakemake <- readRDS(".drop/tmp/AE/counts.snakemake")

suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(BiocParallel)
  library(GenomicAlignments)
})

# get count settings
strand     <- snakemake@params$strand
count_mode <- snakemake@params$count_mode
paired_end <- as.logical(snakemake@params$paired_end)
overlap    <- as.logical(snakemake@params$overlap)

# infer preprocessing and strand info
preprocess_reads <- NULL
if (strand == "yes") {
  strand_spec <- T
} else if (strand == "no"){
  strand_spec <- F
} else if (strand == "reverse") {
  # set preprocess function for later
  preprocess_reads <- invertStrand
  strand_spec <- T
} else {
  stop(paste("invalid strand information", strand))
}

# read files
bam_file <- BamFile(snakemake@input$sample_bam, yieldSize = 2e6)
count_ranges <- readRDS(snakemake@input$count_ranges)
# set chromosome style
seqlevelsStyle(count_ranges) <- seqlevelsStyle(bam_file)

# show info
message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))
message(paste('\tcount mode:', count_mode, sep = "\t"))
message(paste('\tpaired end:', paired_end, sep = "\t"))
message(paste('\toverlap:', overlap, sep = "\t"))
message(paste('\tstrand:', strand, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))
message(paste(seqlevels(count_ranges), collapse = ' '))

# start counting
message("\ncounting")
starttime <- Sys.time()
se <- summarizeOverlaps(
  count_ranges
    , bam_file
    , mode = count_mode
    , singleEnd = !paired_end
    , ignore.strand = !strand_spec  # FALSE if done strand specifically
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = !overlap # TRUE: reads mapping to multiple features are dropped
    , preprocess.reads = preprocess_reads
    , BPPARAM = MulticoreParam(snakemake@threads)
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # print time taken
print(sum(assay(se))) # total counts

