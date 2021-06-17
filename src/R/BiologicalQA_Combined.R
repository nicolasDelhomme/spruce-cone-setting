#' ---
#' title: "Spruce cone development Biological QA combined sets"
#' author: "Veronika Nordal & Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/onilsson/cone-development")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/onilsson/cone-development")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")

#' Create a palette
pal <- brewer.pal(8,"Dark2")

#' Register the default plot margin
mar <- par("mar")

#' Read the sample information
#  Read both sample files
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/samples.csv",
                    row.names=1,as.is=TRUE)
samples.2 <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/samples_2.csv",
                    row.names=1,as.is=TRUE)

# Combine the sample files
samples <- rbind(samples,samples.2)

#' Read the data
set1 <- read.csv("set1/analysis/HTSeq/raw-unormalised_tech-rep-combined_data.csv",row.names = 1)
# The same for set 2
set2 <- read.csv("set2/analysis/HTSeq/raw-unormalised-data.csv",row.names = 1)

# combine them by column
count.table <- cbind(set1,set2)

#' # Data normalisation 
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' Create the dds object
conditions <- colnames(count.table)
colnames(count.table) <- make.unique(colnames(count.table))
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=conditions),
  design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
#' There is no big variation, a Variance Stabilizing Transformation can
#' be used (over a Relative Log Transformation)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(count.table)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#' Perform the VST
colData(dds)$condition <- factor(colData(dds)$condition,
                                 levels=unique(conditions))
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(count.table)
vst <- vst - min(vst)
dir.create("analysis/HTSeq",showWarnings = FALSE, recursive = TRUE)
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd relationship.
meanSdPlot(vst[rowSums(count.table)>0,])

#' # QC on the normalised data
#' 
#' ## PCA
#' 
#' First perform a Principal Component Analysis (PCA) of the data
#'to do a quick quality assessment; i.e. replicate should cluster
#' and the first 2-3 dimensions shouldbe explainable by biological means.
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Plot the PCA 3 first dimensions, They are colored by sample 
#' (using repeating colors)
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(factor(samples$Sample))],
              pch=17)
par(mar=mar)

#' Then the first two dimensions
#' The first dimension separates... 
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Sampling))],
     pch=samples$Sex,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottom", cex=0.75,
       col=pal[1:nlevels(factor(samples$Sampling))],
       pch=19,
       legend=levels(factor(samples$Sampling)))
#' And the 2nd and 3rd dims...
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$Sampling))],
     pch=samples$Sex,
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("bottomleft", inset=0.1,
       col=pal[1:nlevels(factor(samples$Sampling))],
       pch=19,
       legend=levels(factor(samples$Sampling)))
#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.8, labCol=paste(samples$Sex,samples$Sampling,sep="-") )

#' # Conclusion
#' Both sets are affected by the time of sampling and the tissue type. However,
#' their relative importance is different. The first set is affected first by tissue
#' type and then by sampling time, whereas the second is affected the other way around.
#' Also there seem to be a batch effect, samples have been preped at different times and
#' sequenced in different locations. This affects the clustering and we should keep it in
#' mind for the DE (check if it affects the results).
#' ```{r empty, eval=FALSE, echo=FALSE}
#' ```
#' # Exporting the data 
colData(dds) <- cbind(colData(dds),samples[,c("Sex","Sampling","Sample")])
colData(dds)$Sex <- factor(colData(dds)$Sex)
colData(dds)$Sampling <- factor(colData(dds)$Sampling)
design(dds) <- ~Sampling * Sex
dir.create("analysis/DESeq2",showWarnings = FALSE)
save(dds,file = "analysis/DESeq2/deseq-sampling_time-tissue_type.rda")
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
