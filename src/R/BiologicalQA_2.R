#' ---
#' title: "Spruce cone development Biological QA"
#' author: "Veronika Nordal and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/onilsson/cone-development/set2")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/onilsson/cone-development/set2")
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
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/samples_2.csv",
                    row.names=1,as.is=TRUE)

#' Read the HTSeq files in a matrix
res <- mclapply(dir("htseq",pattern="*.txt",full.names=TRUE),function(fil){
  read.delim(fil,header=FALSE,stringsAsFactors=FALSE)
},mc.cores=max(mcaffinity()))

names(res) <- sub("_S[0-9]+$","",sub("_sortmerna.*\\.txt","",dir("htseq",pattern="*.txt")))

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

#' The average percentage of aligned reads is...
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
#' The cumulative coverage is...
plot(density(log10(rowMeans(count.table))),col=pal[1],
     main="mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples
#' 
#' The observed distribution is... 
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
#' (using repeating colors).
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
legend("bottomleft",
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
legend("bottomleft",
       col=pal[1:nlevels(factor(samples$Sampling))],
       pch=19,
       legend=levels(factor(samples$Sampling)))
#' ## Heatmap
#' The 1000 most variable genes are selected and plotted as a heatmap
dimnames(vst) <- dimnames(count.table)
sel <- order(apply(vst,1,sd),decreasing=TRUE)[1:1000]
heatmap.2(vst[sel,],labRow = NA,trace = "none",cexCol = 0.8, labCol=paste(samples$Sex,samples$Sampling,sep="-") )
# # First Conclusion
# ...

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
