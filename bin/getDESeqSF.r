#!/usr/bin/env Rscript
.libPaths(setdiff(.libPaths(), normalizePath(Sys.getenv("R_LIBS_USER"))))

args <- commandArgs(TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

## Calculate scaling factor based on DESeq2 method
stopifnot(require(DESeq2))
count.table <- args[1]
d <- read.table(count.table, comment.char="", header=TRUE, sep="\t", check.names=FALSE)
d <- d[,-c(1:3), drop=FALSE]
colnames(d) = sub("\\_.*", "", colnames(d))
dds <- DESeqDataSetFromMatrix(d, colData=DataFrame(type=rep("chip", ncol(d))), design=~ 1)
dds <- estimateSizeFactors(dds)

write.table(1/sizeFactors(dds), file=gsub(".tab", ".sf", count.table), quote=FALSE, col.names=FALSE, sep=',')
