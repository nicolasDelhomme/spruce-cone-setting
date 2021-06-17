#' ---
#' title: "Spruce Cone development all samples Biological QA"
#' author: "Nicolas Delhomme and Veronika Nordal"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/u2015029/20170824")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/u2015029/20170824")
#' ```

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create palettes
pal <- brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv2("~/Git/UPSCb/projects/spruce-cone-development/doc/All_samples.csv")

#' ### Original data
orig <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(orig) <- sub("_sortmerna.*","",
                    sapply(strsplit(orig, "/"), .subset, 2))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(orig),samples$Complete_SciLifeID),]

#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = orig, 
                                type = "kallisto", 
                                txOut = TRUE))
kg <- round(tx$counts)

#' ## Raw Data QC analysis
#' ### Original data

#' Check how many genes are never expressed
sel <- rowSums(kg) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kg),digits=1),
        sum(sel),
        nrow(kg))

#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(kg))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by sex. 
plot.multidensity(lapply(1:ncol(kg),function(k){log10(kg)[,k]}),
                  col=c(1,pal)[as.integer(samples$Sex)],
                  legend.x="topright",
                  legend=levels(samples$Sex),
                  legend.col=c(1,pal)[1:nlevels(samples$Sex)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' and by batch. 
plot.multidensity(lapply(1:ncol(kg),function(k){log10(kg)[,k]}),
                  col=pal[as.integer(samples$Batch)],
                  legend.x="topright",
                  legend=levels(samples$Batch),
                  legend.col=pal[1:nlevels(samples$Batch)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("analysis","kallisto"),showWarnings=FALSE)
write.csv(kg,file="analysis/kallisto/raw-unormalised-gene-expression_data.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ### Original data
#' Create the dds object, without giving any prior on the design
dds.kg <- DESeqDataSetFromMatrix(
  countData = kg,
  colData = samples[,c("Complete_SciLifeID","ID","Biol_repl_ID","Sex","Type","Sampling","Tree")],
  design = ~Complete_SciLifeID)

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
names(sizes.kg) <- colnames(kg)
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")

#' ## Variance Stabilising Transformation
#' At the gene level
vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE)
vst.kg <- assay(vsd.kg)
vst.kg <- vst.kg - min(vst.kg)

#' Validate the VST 
meanSdPlot(vst.kg[rowSums(kg)>0,])

#' Export the vst
write.csv(vst.kg,"analysis/kallisto/library-size-normalized_variance-stabilized_gene-expression_data.csv")

#' # QC on the normalised data
#' 
#' ## PCA
".pca" <- function(vst,fact,lgd="topleft",pal=brewer.pal(8,"Dark2")){
  pc <- prcomp(t(vst))
  
  percent <- round(summary(pc)$importance[2,]*100)
  
  #' ### 3 first dimensions
  mar=c(5.1,4.1,4.1,2.1)
  scatterplot3d(pc$x[,1],
                pc$x[,2],
                pc$x[,3],
                xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                color=pal[as.integer(fact)],
                pch=19)
  legend(lgd,pch=19,
         col=pal[1:nlevels(fact)],
         legend=levels(fact))
  par(mar=mar)
}

#' ### Batch
.pca(vst.kg,factor(samples$Batch),lgd="topright")

#' ### Sex
.pca(vst.kg,factor(samples$Sex),lgd="topright",pal=c(1,pal))

#' ### 1st and 2nd dims
pc <- prcomp(t(vst.kg))
percent <- round(summary(pc)$importance[2,]*100)

#' Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Batch))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal,levels(factor(samples$Batch)),pch=19)

#' Sex
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=c(1,pal)[as.integer(factor(samples$Sex))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=c(1,pal),levels(samples$Sex),pch=19)

#' Sampling
samples$Date <- 
  factor(sapply(lapply(strsplit(sub("Oct","10",
               sub("Sep","09",
                   sub("Aug","08",
                       sub("Jun","06",
                           sub("May","05",
                               sub("Apr","04",samples$Sampling)))))),"-"),rev),
         paste,collapse="-"))

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal12[as.integer(factor(samples$Date))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal12,levels(samples$Date),pch=19)

#' ### 2nd and 3rd dims
#' Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$Batch))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal,levels(factor(samples$Batch)),pch=19)

#' Sex
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(1,pal)[as.integer(samples$Sex)],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=c(1,pal),levels(samples$Sex),pch=19)

#' Sampling
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal12[as.integer(factor(samples$Date))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal12,levels(samples$Date),pch=19)

#' ## Heatmap
sel <- featureSelect(vst.kg,samples$Biol_repl_ID,exp = 5,nrep = 2)
heatmap.2(vst.kg[sel,],trace = "none",labRow = FALSE, 
          col=hpal, ColSideColors = c(1,pal)[as.integer(samples$Sex)],
          scale="row")

#' ### Saturate the expression
vst.sat <- vst.kg[sel,]
vst.sat[vst.sat<3] <- 3
vst.sat[vst.sat>8] <- 8

heatmap.2(vst.sat,trace = "none",labRow = FALSE, 
          col=hpal, ColSideColors = c(1,pal)[as.integer(samples$Sex)])

#' ### Sample Hierarchical clustering
hc <- hclust(dist(t(vst.sat)))
plot(hc, labels=samples$ID,
     main = "Hierarchical clustering",cex=0.7)

#' # Merge the technical replicates
counts <- do.call(
  cbind,
  lapply(split.data.frame(t(kg),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),csamples$ID),]
write.csv(csamples,file="~/Git/UPSCb/projects/spruce-cone-development/doc/biological-samples.csv")

#' ## QC 
#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(counts))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by sex. 
plot.multidensity(lapply(1:ncol(counts),function(k){log10(counts)[,k]}),
                  col=c(1,pal)[as.integer(csamples$Sex)],
                  legend.x="topright",
                  legend=levels(csamples$Sex),
                  legend.col=c(1,pal)[1:nlevels(csamples$Sex)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' and by batch. 
plot.multidensity(lapply(1:ncol(counts),function(k){log10(counts)[,k]}),
                  col=pal[as.integer(csamples$Batch)],
                  legend.x="topright",
                  legend=levels(csamples$Batch),
                  legend.col=pal[1:nlevels(csamples$Batch)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Write the count table (tech. rep. combined)
write.csv(counts,file="analysis/kallisto/raw-unormalised-gene-expression-tech-rep-combined_data.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ### Original data
#' Create the dds object, without giving any prior on the design
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = csamples[,c("ID","Biol_repl_ID","Sex","Type","Sampling","Date","Batch","Tree")],
  design = ~Biol_repl_ID)

#' Check the size factors (i.e. the sequencing library size effect)
dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(counts)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")
pander(sort(sizes))

#' ## Variance Stabilising Transformation
#' At the gene level
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Validate the VST 
meanSdPlot(vst[rowSums(vst)>0,])

#' Export the vst
write.csv(vst,"analysis/kallisto/library-size-normalized_variance-stabilized_gene-expression-tech-rep-combined_data.csv")

#' # QC on the normalised data
#' 
#' ## PCA
".pca" <- function(vst,fact,lgd="topleft",pal=brewer.pal(8,"Dark2")){
  pc <- prcomp(t(vst))
  
  percent <- round(summary(pc)$importance[2,]*100)
  
  #' ### 3 first dimensions
  mar=c(5.1,4.1,4.1,2.1)
  scatterplot3d(pc$x[,1],
                pc$x[,2],
                pc$x[,3],
                xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                color=pal[as.integer(fact)],
                pch=19)
  legend(lgd,pch=19,
         col=pal[1:nlevels(fact)],
         legend=levels(fact))
  par(mar=mar)
}

#' ### Batch
.pca(vst,factor(csamples$Batch),lgd="topright")

#' ### Sex
.pca(vst,factor(csamples$Sex),lgd="topright",pal=c(1,pal))

#' ### 1st and 2nd dims
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(csamples$Batch))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal,levels(factor(csamples$Batch)),pch=19)

#' Sex
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=c(1,pal)[as.integer(factor(csamples$Sex))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=c(1,pal),levels(csamples$Sex),pch=19)

#' Sampling
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal12[as.integer(factor(csamples$Date))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal12,levels(csamples$Date),pch=19)

#' ### 2nd and 3rd dims
#' Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(csamples$Batch))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal,levels(factor(csamples$Batch)),pch=19)

#' Sex
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=c(1,pal)[as.integer(csamples$Sex)],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=c(1,pal),levels(csamples$Sex),pch=19)

#' Sampling
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal12[as.integer(factor(csamples$Date))],
     pch=19,
     main="Principal Component Analysis",sub="variance stabilized counts")

legend("topright",bty="n",col=pal12,levels(csamples$Date),pch=19)

#' ## Heatmap
sel <- featureSelect(vst,csamples$Biol_repl_ID,exp = 5,nrep = 2)
heatmap.2(vst[sel,],trace = "none",labRow = FALSE, 
          col=hpal, ColSideColors = c(1,pal)[as.integer(csamples$Sex)],
          scale="row")

#' ### Saturate the expression
vst.sat <- vst[sel,]
vst.sat[vst.sat<3] <- 3
vst.sat[vst.sat>8] <- 8

heatmap.2(vst.sat,trace = "none",labRow = FALSE, 
          col=hpal, ColSideColors = c(1,pal)[as.integer(csamples$Sex)])

#' ### Sample Hierarchical clustering
hc <- hclust(dist(t(vst.sat)))
plot(hc, labels=csamples$ID,
     main = "Hierarchical clustering",cex=0.7)

plot(hc, labels=round(sizes[csamples$ID],digits = 2),
     main = "Hierarchical clustering",cex=0.7)

#' # Conclusion
#' The data quality looks overall good. There is some covariates that will need to be taken into
#' account in the analysis, principaly the sampling date. This is confounded with the batch effect.
#' Some direct comparisons will be counfonded with time and it might be difficult to disentangle that 
#' effect. We could try to block it on the whole dataset, but we will have an unbalanced design.
#' 
#' The technical replicates...
#'  
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
