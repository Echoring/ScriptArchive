#!/usr/bin/env Rscript
# $1=count.matrix $2=condition $3=minGroupSize $4=minCount $5=outFile
# Don't use - in sample name.
suppressPackageStartupMessages(suppressWarnings({
    library("DESeq2")
}))
args <- (commandArgs(TRUE))

# Import data
cts <- read.table(args[1], header=TRUE, sep="\t", row.name=1)
cts <- round(cts)
coldata <- read.table(args[2], header=TRUE, sep="\t", row.name=1)
if (!all(rownames(coldata) %in% colnames(cts))) {
    print("Sample name not match in matrix and condition.")
    q()
}
cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ condition)

# Filter
smallestGroupSize <- args[3]
keep <- rowSums(counts(dds) >= args[4]) >= smallestGroupSize
dds <- dds[keep,]

# Run
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$pvalue),]
write.csv(as.data.frame(resOrdered), file=args[5])
