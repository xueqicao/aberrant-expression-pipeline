#'---
#' title: Count reads
#' author: Michaela Mueller
#' wb:
#'  input:
#'   - sample_bam: '`sm lambda wildcards: parser.getFilePath(wildcards.sampleID, "rna_assay") `'
#'   - features: '`sm lambda wildcards: parser.getCountRangesFile(wildcards.annotation) `'
#'  output:
#'   - counts: '`sm parser.getProcResultsDir() + "/{annotation}/counts/{sampleID,[^/]+}.Rds"`'
#'  type: script
#'---

# #'   - sample_bam: '`sm standardFileNames("Data/helmholtz/{sampleID}/RNAout/paired-endout/stdFilenames/{sampleID}.bam")`'
# #'   - sample_bam: '`sm lambda wildcards: parser.getFilePath(wildcards.sampleID,config["rna_assay"]) `'



#source(".wBuild/wBuildParser.R")
#parseWBHeader("Scripts/counting/countReads.R")
saveRDS(snakemake, "tmp/counts.snakemake")
# snakemake <- readRDS("tmp/counts.snakemake")
suppressPackageStartupMessages({
    library(data.table)
})

# import count settings from config
anno <- snakemake@wildcards$annotation

# import sample annotation
sampleID <- snakemake@wildcards$sampleID
sample_anno <- fread(snakemake@config$SAMPLE_ANNOTATION)
# Get strand specific information from sample annotation

rna_assay <- snakemake@config$rna_assay
strand_spec <- sample_anno[get(rna_assay) == sampleID, get(snakemake@config$strand_column)]
inter_feature <- sample_anno[get(rna_assay) == sampleID, get(snakemake@config$inter_feature_column)]

# show info
message(paste("input:", snakemake@input$features))
message(paste("output:", snakemake@output$counts))
message(paste('\tinter.feature:', inter_feature, sep = "\t"))
message(paste('\tstrand specific:', strand_spec, sep = "\t"))


# read files
bam_file <- Rsamtools::BamFile(snakemake@input$sample_bam, yieldSize = 2e6)
feature_regions <- readRDS(snakemake@input$features)

# start counting
message("counting")
starttime <- Sys.time()
se <- GenomicAlignments::summarizeOverlaps(
    feature_regions
    , bam_file
    , mode = 'IntersectionStrict'
    , singleEnd = F
    , ignore.strand = !strand_spec  # FALSE if done strand specifically
    , fragments = F
    , count.mapped.reads = T
    , inter.feature = inter_feature # TRUE, reads mapping to multiple features are dropped
)
saveRDS(se, snakemake@output$counts)
message("done")

print(format(Sys.time()- starttime)) # time taken
print(sum(SummarizedExperiment::assay(se))) # total counts

