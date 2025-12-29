#!/usr/bin/env Rscript
# Usage: Rscript enricher.R <genelist> <geneclass> <pvalue> <prefix>
suppressPackageStartupMessages(suppressWarnings({
    library(dplyr)
    library(clusterProfiler)
    library(ggplot2)
    library(stringr)
}))

args <- (commandArgs(TRUE))
genelist <- args[1]
geneclass <- args[2]
pvalue <- as.numeric(args[3])
prefix <- args[4]

gene <- read.table(genelist, header=FALSE, stringsAsFactors=FALSE) %>% pull(V1)
term2gene <- read.table(geneclass, header=FALSE, stringsAsFactors=FALSE) %>% select(V2, V1) %>% rename(term = V2, gene = V1)
enrichment_result <- enricher(gene, TERM2GENE=term2gene, pvalueCutoff=pvalue)
write.table(as.data.frame(enrichment_result), file=paste0(prefix,".enrichment_result.tsv"), sep="\t", row.names=FALSE, quote=FALSE)

test=as.data.frame(enrichment_result)
test=test[1:min(10, nrow(enrichment_result)),]
gr1 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,1])
gr2 <- as.numeric(str_split(test$GeneRatio,"/",simplify = T)[,2])
test$fold <- (gr1/gr2)
test <- arrange(test,fold)
test$Description = factor(test$Description,levels = test$Description,ordered = T)
pdf(file=paste0(prefix,".enrichment_dotplot.pdf"))
ggplot(test,aes(x = fold, y = Description)) +
    geom_point(aes(color = pvalue,size = Count)) +
    scale_color_gradient(low = "#EB746A", high = "#43BBBC") + # change color here
    theme_bw() +
    guides(color = guide_colorbar(reverse = TRUE)) +
    xlab("GeneRatio")
dev.off()

pdf(file=paste0(prefix,".enrichment_barplot.pdf"))
barplot(enrichment_result,showCategory=10,color="pvalue")
dev.off()
pdf(file=paste0(prefix,"enrichment_dotplot.pdf"))
dotplot(enrichment_result,showCategory=10,color="pvalue")
dev.off()
