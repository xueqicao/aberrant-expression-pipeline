#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - sample_bam: '`sm lambda wildcards: parser.getFilePath(wildcards.sampleID, "rna_assay") `'
#'   - count_ranges: '`sm parser.getProcDataDir() + "/{annotation}/count_ranges.Rds" `'
#'  output:
#'   - counts: '`sm parser.getProcDataDir() + "/{annotation}/counts/{sampleID,[^/]+}.Rds"`'
#'  type: script
#'  threads: 5
#'---

saveRDS(snakemake, "tmp/counts.snakemake")
# snakemake <- readRDS("tmp/counts.snakemake")
suppressPackageStartupMessages({
  library(data.table)
  library(Rsamtools)
  library(BiocParallel)
  library(GenomicAlignments)
})

# Get strand specific information from sample annotation
sampleID <- snakemake@wildcards$sampleID
sample_anno <- fread(snakemake@config$SAMPLE_ANNOTATION)
sample_anno <- sample_anno[get(snakemake@config$rna_assay) == sampleID]

count_mode <- sample_anno[, get(snakemake@config$count_mode_column)]
paired_end <- sample_anno[, get(snakemake@config$paired_end_column)]
inter_feature <- sample_anno[, get(snakemake@config$inter_feature_column)]
strand <- sample_anno[, get(snakemake@config$strand_column)]

# infer preprocessing and strand info
preprocess_reads <- NULL
if (strand == "yes") {
  strand_spec <- T
} else if (strand == "no"){
  strand_spec <- F
} else if (strand == "reverse") {
  preprocess_reads <- invertStrand
  strand_spec <- T
} else {
  stop(paste("invalid strand information", strand))
}

# read files
bam_file <- BamFile(snakemake@input$sample_bam, yieldSize = 2e6)
count_ranges <- readRDS(snakemake@input$count_ranges)

# show info
message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))
message(paste('\tcount mode:', count_mode, sep = "\t"))
message(paste('\tpaired end:', paired_end, sep = "\t"))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
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
    , inter.feature = inter_feature # TRUE, reads mapping to multiple features are dropped
    , preprocess.reads = preprocess_reads
    , BPPARAM = MulticoreParam(snakemake@threads)
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(assay(se))) # total counts


