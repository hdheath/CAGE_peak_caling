#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
bam_list     <- strsplit(args[1], ",")[[1]]
out_prefix   <- args[2]
minThreshold <- as.numeric(args[3])
maxDist      <- as.numeric(args[4])

# load CAGEr
if (! requireNamespace("CAGEr", quietly=TRUE)) {
  stop("CAGEr package not found—please install it in your cager-env.")
}
suppressPackageStartupMessages(library(CAGEr))

# genome package
genomePkg <- "BSgenome.Scerevisiae.UCSC.sacCer3"
if (! requireNamespace(genomePkg, quietly=TRUE)) {
  stop(paste0("Genome package '", genomePkg, "' not installed—please install with BiocManager."))
}

# construct
if ("CAGEset" %in% ls("package:CAGEr")) {
  cs <- CAGEset(inputFiles      = bam_list,
                inputFilesType  = "bam",
                sampleLabels    = out_prefix,
                genomeName      = genomePkg)
} else {
  cs <- CAGEexp(inputFiles      = bam_list,
                inputFilesType  = "bam",
                sampleLabels    = out_prefix,
                genomeName      = genomePkg)
}

# 1) CTSS extraction + raw-count normalization
cs <- getCTSS(cs)
cs <- normalizeTagCount(cs, method = "none")

# 2) clustering on raw counts
cs <- clusterCTSS(
  object           = cs,
  threshold        = minThreshold,
  thresholdIsTpm   = FALSE,
  method           = "distclu",
  maxDist          = maxDist,
  removeSingletons = TRUE
)

# 3) extract clusters (works on CAGEset or CAGEexp)
if ("tagClustersGR" %in% ls("package:CAGEr")) {
  tc <- tagClustersGR(cs, sample = out_prefix)
} else {
  tc <- tagClusters(cs)
}

# 4) write BED
out_bed <- file.path(getwd(), paste0(out_prefix, "_CAGEr_clusters.bed"))
bed <- data.frame(
  seqnames = as.character(seqnames(tc)),
  start    = start(tc),
  end      = end(tc),
  name     = names(tc),
  score    = mcols(tc)$score,      # total tag count per cluster
  strand   = as.character(strand(tc))
)
write.table(bed, file = out_bed, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
cat("[CAGEr] Wrote tag clusters to", out_bed, "\n")
