#' ---
#' title: "Spruce cone development Biological QA"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/onilsson/cone-development/set1")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/onilsson/cone-development/set1")
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
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/samples.csv",
                    row.names=1,as.is=TRUE)

#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

names(res) <- sub(".*_P1387","P1387",sub("_sortmerna.*\\.txt","",dir("htseq",pattern="*.txt")))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(res),samples$SciLifeID),,drop=FALSE]

#' Raw Data QC analysis
addInfo <- c("__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")
sel <- match(addInfo,res[[1]][,1])
count.table <- do.call(cbind,lapply(res,"[",3))[-sel,]
colnames(count.table) <- names(res)
rownames(count.table) <- res[[1]][,1][-sel]
dir.create(file.path("analysis","HTSeq"),recursive = TRUE, showWarnings = FALSE)
write.csv(count.table,"analysis/HTSeq/raw-unormalised-data.csv")

#' Extract the HTSeq stat lines
count.stats <- do.call(cbind,lapply(res,"[",2))[sel,]
colnames(count.stats) <- names(res)
rownames(count.stats) <- sub("__","",addInfo)
count.stats <- rbind(count.stats,aligned=colSums(count.table))
count.stats <- count.stats[rowSums(count.stats) > 0,]

#' Convert them into percentages
pander(apply(count.stats,2,function(co){round(co*100/sum(co))}))

#' Plot the stats
#' 
#' There are no outliers
col <- pal[1:nrow(count.stats)]
par(mar=c(7.1,5.1,4.1,2.1))
barplot(as.matrix(count.stats),col=col,beside=TRUE,las=2,main="read proportion",
        ylim=range(count.stats) + c(0,4e+6),cex.names=.6)
legend("top",fill=col,legend=gsub("_"," ",rownames(count.stats)),bty="n",cex=0.8)
par(mar=mar)

#' The average percentage of aligned reads is 72%
# round(mean(unlist(count.stats["aligned",]/colSums(count.stats))),digits=2)*100
boxplot(unlist(count.stats["aligned",]/colSums(count.stats)),
        main="aligned reads",ylab="percent aligned",ylim=c(0,1))

#' Check how many genes are never expressed
sel <- rowSums(count.table) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(count.table),digits=1),
        sum(sel),
        nrow(count.table))

#' Display the per-gene mean expression
#' 
#' i.e. the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
#' The cumulative coverage is as expected, around 100X
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples
#' 
#' The observed distribution is strikingly similar. 
cols <- rep(sample(pal,length(unique(colnames(count.table))),TRUE),2)
plot.multidensity(log10(count.table),
                  col=cols,
                  legend="",
                  legend.col=NA,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

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
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv")

#' Validate the VST 
#' 
#' Visualize the corrected mean - sd relationship. It is amazingly linear,
#' meaning we can assume homoscedasticity.
#' The slight initial trend / bump is due to genes having few counts in
#' a few subset of the samples and hence having a higher variability. This is
#' expected.
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
#' (using repeating colors), replicates are dots and triangles
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=rep(c(rep(pal,3),1),2),
              pch=rep(c(17,19),each=25))
par(mar=mar)

#' Then the first two dimensions
#' The first dimension separates two set of samples. 
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=rep(c(rep(pal,3),1),2),
     pch=rep(c(17,19),each=25),
     main="Principal Component Analysis",sub="variance stabilized counts")

#' And the 2nd and 3rd dims possibly separate samples on the 3rd dimension also
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=rep(c(rep(pal,3),1),2),
     pch=rep(c(17,19),each=25),
     main="Principal Component Analysis",sub="variance stabilized counts")

#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
#' 
#' At an high level there is clearly 2 separate groups, then at a finer level,
#' there are up to 5 groups. The replicates nicely group together. It is obvious
#' from the scale of the dendrogram, where the technical repliation error appears
#' quite high, that there may not be a large difference overall between all samples
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

#' # First Conclusion
#' The raw quality of the data appears very good. The data normalisation
#' gives satisfying results (as far as the VST fit is concerned). The PCA and 
#' the heatmap identify 2 to 5 groups, which needs to be re-evaluated once we have
#' the sample information. Finally, it is clear from the analysis that the technical 
#' replicates can be combined together, which is what we do next.
count.table <- do.call(cbind,lapply(split.data.frame(t(count.table),samples$SciLifeID),colSums))
write.csv(count.table,"analysis/HTSeq/raw-unormalised_tech-rep-combined_data.csv")

#' # QC using the sample information
#' Now that we have the sample information about the sex, the sampling date and the
#' replicate number we can redo the QC on these.
samples <- samples[match(colnames(count.table),samples$SciLifeID),]
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=colnames(count.table),
                       sex=samples$Sex,
                       date=samples$Sampling,
                       rep=samples$Sample),
  design = ~ condition)

#' Check the size factors (i.e. the sequencing library size effect)
#' Merging the technical replicates had little effect and as
#' there is still no big variation in the size factor, a Variance Stabilizing Transformation can
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
write.csv(vst,"analysis/HTSeq/library-size-normalized_variance-stabilized_tech-rep-combined_data.csv")

#' Redo the heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.6 )

#' labelled by sex
heatmap.2(vst[sel,],labRow = NA,trace = "none", labCol=samples$Sex )

#' labelled by sex and scaled
heatmap.2(vst[sel,],labRow = NA,trace = "none", labCol=samples$Sex,
          scale="row")

#' labelled by ID
heatmap.2(vst[sel,],labRow = NA,trace = "none", labCol=samples$ID)

#' Finally a PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Sex))],
     pch=c(17,19)[as.integer(factor(samples$Sampling))],
     main="Principal Component Analysis",sub="variance stabilized counts")
legend("center",pch=c(17,rep(19,4)),
       col=c(1,1,pal[1:3]),
       legend=c(levels(factor(samples$Sampling)),
                levels(factor(samples$Sex))))

#' PCA looking at the replicate samples, 1st and 2nd dimension
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(factor(samples$Sex))],
     pch=as.character(samples$Sample),
     main="Principal Component Analysis",sub="variance stabilized counts")

#' PCA looking at the replicate samples, 2nd and 3rd dimension
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(factor(samples$Sex))],
     pch=as.character(samples$Sample),
     main="Principal Component Analysis",sub="variance stabilized counts")

#' # Second Conclusion
#' The biological QA of the data looks very promising
#' 
#' The data clusters primarily by Sex, with the male samples being most different
#' from the Vegetative and Female samples. The sampling date seem to have a limited 
#' but visible importance, especially in the PCA. For the 1000 most variable genes - 
#' in the heatmap - the effect is confounded as sample replication over time (e.g. M1 in 10-08 or 09-12)
#' clusters together most frequently than not. Sample replicate 3 and 5 (M or F), 
#' however still show a marked difference over time; e.g. sample M3 and M5 cluster by 
#' sampling date, rather than by replication. Also, interestingly, the vegetative
#' samples cluster with the female sample also according to their replication number;
#' i.e. V3 and V5 clusters with F3 and F5.
#' Finally, the PCA 2nd and 3rd dimension clearly separates the samples 1,2,4 from
#' the samples 3,5. This is certainly not just due to chance, and assuming all
#' samples came from a single tree, might be related to their location on the tree.
#'
#' As an outlook, here is a list of comparision that could be (or not) done:
#' 
#' (i) comparing male vs. female sample is going to reveal the most
#' differentially expressed genes. To increase power, such a design should make
#' sure to block the sampling date effect.
#' 
#' (ii) as above comparing male vs. vegetative
#' 
#' (iii) comparing vegetative with female is most likely to unravel only noise or
#' very little DE genes.
#' 
#' (iv) comparing the sampling date is going to reveal DE genes (there is already 
#' a trend visible in the first component for Male samples - and the 2nd and 3rd
#' component, show also in addition to the 1,2,4 vs. 3,5 sample grouping, an effect
#' of time). Obviously, the effect of sex and sample would have to be blocked, which 
#' may not be optimal in terms of power.
#' 
#' (v) comparing the sample set 1,2,4 and 3,5, as (iv)
#' 

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
