#' * Libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(here)
  library(pheatmap)
  library(RColorBrewer)
  library(VennDiagram)
})

#' * Helpers
source(here("volcanoPlot.R"))

#' * Graphics
pal <- brewer.pal(8,"Dark2")  

#' * Load the data
load(here("dds.6h.rda"))

#' # Differential expression analysis
dds <- DESeq2::DESeq(dds)
DESeq2::plotDispEsts(dds)

design(dds)
colData(dds)[,c("condition","replicate")]

#possible contrasts
resultsNames(dds)

#concave (repeat for convex)
res <- results(dds,name="condition_concave_vs_control")

levels(dds$condition)
levels(dds$replicate)

head(res)
mcols(res, use.names = TRUE)
summary(res)
hist(res$pvalue,breaks=seq(0,1,.01))

#' Following Schurch _et al._, (RNA, 2016) recommandations, we would have for 3 replicates per conditions:
resSchurch <- results(dds, name="condition_concave_vs_control", lfcThreshold = 0.5, alpha = 0.01)
summary(resSchurch)

#save tables_cv (repeat for convex)
resSig <- subset(resSchurch, resSchurch$padj < 0.01)
head(resSig)
write.table(resSig, "cv.6h.sig.txt", sep = ';', col.names = NA)

resSig_up<-subset(resSig, resSig$log2FoldChange>0.5)
resSig_down<-subset(resSig, resSig$log2FoldChange<0.5)

write.table(resSig_up, "cv.6h.up.txt", sep = ';',col.names = NA)
write.table(resSig_down, "cv.6h.down.txt", sep = ';',col.names = NA)

#' ### Visualisation
volcanoPlot(resSchurch)

#' #### Individual genes
plotCounts(dds, gene = "XM_002325714.4", intgroup = "condition", 
           normalized = TRUE, transform = FALSE)

#' ## MA plot with DESeq2
DESeq2::plotMA(res, ylim = c(-5, 5))

#' ## Heatmap of the most significant genes
#correcting for library size (variance-stabilized values)
vsd <- DESeq2::vst(dds)
mat <- assay(vsd)[head(order(res$padj), 100), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("condition", "replicate")])
pheatmap(mat, annotation_col = df)
