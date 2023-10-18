#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(emoji)
  library(gplots)
  library(gtools)
  library(here)
  library(hyperSpec)
  library(limma)
  library(magrittr)
  library(parallel)
  library(patchwork)
  library(PCAtools)
  library(pheatmap)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(vsn)
})

#' * Helper functions
source(here("featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("cyan","magenta","yellow"))(100)

#' # Data
#' * Sample information
samples <- read_tsv(here("Sample_metadata.txt"),
                    col_types=cols(.default=col_factor()))

#' * Raw data
filelist <- list.files(here("/Users/mohamedkouhen/Library/CloudStorage/OneDrive-UniversitàdegliStudidelMolise/UNIMOL/RNA sequencing/Data Analysis/pop_6hbs/QA_DESeq2"), 
                          pattern = "*.tabular",
                          full.names = TRUE)

#' * Sanity check to ensure that the data is sorted according to the sample info
stopifnot(all(sub("CONTROL","ctrl",
                  toupper(sub("\\.","",
                              gsub(".*_|\\.tabular","",basename(filelist))))) == 
                samples$sample_name))

#' * add filelist to samples as a new column
samples %<>% mutate(Filenames = filelist) %>% rename(SampleID="sample_name")

#' Read the expression at the gene level
txi <- suppressMessages(tximport(files = samples$Filenames,
                                 type = "salmon",
                                 txOut=TRUE))
counts <- txi$counts
# replace NA by zero
counts[is.na(txi)] <- 0
head(counts)

colnames(counts) <- samples$SampleID

#' # Quality Control
#' * "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Sequencing depth
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

# Check for rows with missing values
rows_with_missing <- which(rowSums(is.na(dat)) > 0)

ggplot(dat,aes(x,y,fill=condition)) + 
  geom_col() + 
  scale_y_continuous(name="Reads") +
  facet_grid(~ factor(replicate), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' ## per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 

ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' ## Per-sample expression

dat <- as.data.frame(log10(counts)) %>% 
  utils::stack() %>% 
  mutate(condition=samples$condition[match(ind,samples$SampleID)])

ggplot(dat,aes(x=values,group=ind,col=condition)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' * Export raw expression data
dir.create(here("/Users/mohamedkouhen/Library/CloudStorage/OneDrive-UniversitàdegliStudidelMolise/UNIMOL/RNA sequencing/Data Analysis/pop_6hbs/DESeq2"),showWarnings=FALSE,recursive=TRUE)
write.csv(counts,file=here("//Users/mohamedkouhen/Library/CloudStorage/OneDrive-UniversitàdegliStudidelMolise/UNIMOL/RNA sequencing/Data Analysis/pop_6hbs/DESeq2/raw-unormalised-gene-expression_data.csv"))
#' 
#' <hr />
#' &nbsp;
#' 
#' # Data normalisation 
#' ## Preparation
#' 
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
dds <- DESeqDataSetFromTximport(
  txi=txi,
  colData = samples,
  design = ~ condition+replicate)

colnames(dds) <- samples$SampleID

save(dds,file=here("/Users/mohamedkouhen/Library/CloudStorage/OneDrive-UniversitàdegliStudidelMolise/UNIMOL/RNA sequencing/Data Analysis/pop_6hbs/QA_DESeq2/dds.6h.rda"))

#' ## size factors 
#' (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds) %>% 
  suppressMessages()

boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")
abline(h=1, col = "Red", lty = 3)

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' Assess whether there might be a difference in library size linked to a
#' given metadata
boxplot(split(normalizationFactors(dds),dds$condition),las=2,
        main="Sequencing libraries size factor by condition",outline=FALSE)

plot(colMeans(normalizationFactors(dds)),log10(colSums(counts(dds))),
     ylab="log10 raw depth",xlab="average scaling factor",
     col=rainbow(n=nlevels(dds$condition))[as.integer(dds$condition)],pch=19)
legend("bottomright",fill=rainbow(n=nlevels(dds$condition)),
       legend=levels(dds$condition),cex=0.6)

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' ## Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After:
meanSdPlot(vst[rowSums(vst)>0,])

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).
#' 

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Using PCAtools
p <- pca(vst,colData(dds))

#' ### Scree plot
#' 
#' We define the number of variable of the model: ```r 2```
nvar <- 1

#' An the number of possible combinations
nlevel=nlevels(dds$condition)

#' We devise the optimal number of components using two methods
horn <- suppressWarnings(parallelPCA(vst))
elbow <- findElbowPoint(p$variance)

#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent),p=percent),aes(x=x,y=y)) +
  geom_line() + geom_col(aes(x,p)) + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component",breaks=1:length(percent),minor_breaks=NULL) + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",linewidth=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",linewidth=0.5) +
  geom_vline(xintercept=c(horn$n,elbow),colour="black",linetype="dotted",linewidth=0.5) +
  geom_label(aes(x = horn$n + 1, y = cumsum(percent)[horn$n],label = 'Horn', vjust = 1)) +
  geom_label(aes(x = elbow + 1, y = cumsum(percent)[elbow],label = 'Elbow', vjust = 1))


#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    PC3=pc$x[,3],
                    as.data.frame(colData(dds)))

#PC1 vs PC2
p1 <- ggplot(pc.dat,aes(x=PC1,y=PC2,col=condition,shape=replicate,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p1) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
install.packages("ggalt")
biplot(p,
       colby = 'condition',
       colLegendTitle = 'Condition',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)


#PC1 vs PC3
p2 <- ggplot(pc.dat,aes(x=PC1,y=PC3,col=condition,shape=replicate,text=SampleID)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC3 (",percent[3],"%)",sep=""))) %>% suppressWarnings()

#' The same as a biplot
biplot(p,x = 'PC1', y = 'PC3',
       colby = 'condition',
       colLegendTitle = 'condition',
       encircle = TRUE,
       encircleFill = TRUE,
       legendPosition = 'top', 
       legendLabSize = 16, legendIconSize = 8.0)

#' `r emoji("point_right")` ** **

#' ```{r subplot, out.width = '100%'}
#' subplot(style(p1, showlegend = FALSE), p2,
#'         titleX = TRUE, titleY = TRUE, nrows = 1, margin = c(0.05, 0.05, 0, 0))
#' ```

#' ### Pairs plot
#' This allows for looking at more dimensions, five by default
#' 
suppressMessages(pairsplot(p,colby="condition",shape="replicate"))

#' `r emoji("point_right")` **There could be a tentative batch effect in PC2, dots, then squares, then triangles**

#' ### Loadings
#' Loadings, _i.e._ the individual effect of every gene in a component can be studied. Here the most important ones are visualized throughout the different PCs
plotloadings(p,
             rangeRetain = 0.01,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'PC1 to PC5',
             caption = 'Top 1% variables',
             drawConnectors = TRUE)

#' ### Correlation
#' This is a plot showing the correlation between the PC and the model variables. 
#' Note that condition and replicate are categorical variables, not numerical variables. Hence, to achieve a correlation,
#' I used the PCA plots to order them accordingly. This is of course somewhat circular, but it should be nonetheless informative.
#' The replicates are sorted based on component 2. The condition are sorted based on component 1, which happens to be the order of the 
#' categorical variable already
p$metadata$Reffect <- ifelse(as.integer(dds$replicate)==1,3,as.integer(dds$replicate)-1)

#' Plotting the two relevant variables. As expected the Batch effect is associated to PC2, but is also contributing to PC5, and 7.
#' The condition is associated with PC1, and PC7 and 9. 
suppressWarnings(eigencorplot(p,components=getComponents(p,1:9),metavars=c('Reffect','condition')))

#' ### Samples Distance
hpal <- colorRampPalette(c("blue","white","red"))(100)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(sampleDistMatrix) <- dds$SampleID
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=hpal)


#' ## Sequencing depth
#' The figures show the number of genes expressed per condition at different expression cutoffs. The scale on the lower plot is the same as on the upper.
#' The first plot is a heatmap showing the number of genes above a given cutoff. The second plot shows it as a ratio of the number of genes expressed for (a)
#' given variable(s) divided by the average number of genes expressed at that cutoff across all variable(s). The latter plot is of course biased at higher cutoff 
#' as the number of genes becomes smaller and smaller.
#' The point of these two plots is to assert whether the number of genes expressed varies between conditions, as this would break some assumptions for normalisation and differential expression.
conds <- dds$condition
dev.null <- rangeSamplesSummary(counts=vst,
                                conditions=conds,
                                nrep=3)

#' ## Heatmap
#' 
#' Here we want to visualise all the informative genes as a heatmap. We first filter the genes to remove those below the selected noise/signal cutoff. 
#' The method employed here is naive, and relies on observing a sharp decrease in the number of genes within the first few low level of expression. 
#' Using an independent filtering method, such as implemented in DESeq2 would be more accurate, but for the purpose of QA validation, a naive approach is sufficient.
#' Note that a sweet spot for computation is around 20 thousand genes, as building the hierarchical clustering for the heatmap scales non-linearly.
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=3) %>% 
  suppressWarnings()

vst.cutoff <- 1
nn <- nrow(vst[sels[[vst.cutoff+1]],])
tn <- nrow(vst)
pn <- round(nn * 100/ tn, digits=1)

mat <- t(scale(t(vst[sels[[vst.cutoff+1]],])))
hm <- pheatmap(mat,
               color = hpal,
               border_color = NA,
               clustering_distance_rows = "correlation",
               clustering_method = "ward.D2",
               show_rownames = FALSE,
               labels_col = conds,
               angle_col = 90,
               legend = FALSE)

#' ## Clustering of samples
#' ```{r echo=FALSE,eval=FALSE}
#' # Developer: This wouldonly works with the gplots heatmap.2, not the pheatmap
#' plot(as.hclust(hm$colDendrogram),xlab="",sub="")
#' ```
#'
#' 
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 100, parallel = TRUE)

plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)

plot(hm.pvclust, labels = dds$SampleID)
