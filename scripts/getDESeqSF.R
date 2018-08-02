rm(list=ls())

args <- commandArgs(TRUE)
la <- length(args)
if (la > 0){
  for (i in 1:la)
    eval(parse(text=args[[i]]))
}

## Calculate scaling factor based on DESeq2 method
stopifnot(require(DESeq2))
d <- read.table(count.table, comment.char="", header=TRUE, sep="\t", check.names=FALSE)
d <- d[,-c(1:3)]
dds <- DESeqDataSetFromMatrix(d, colData=DataFrame(type=rep("chip", ncol(d))), design=~ 1)
dds <- estimateSizeFactors(dds)

write.table(sizeFactors(dds), file=gsub(".tab", ".sf", count.table), quote=FALSE, col.names=FALSE)
