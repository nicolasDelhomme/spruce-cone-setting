#' ---
#' title: "Tissue and Date heatmaps"
#' author: "Nicolas Delhomme & Veronika Nordal"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' # Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/onilsson/cone-development/20170824")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/onilsson/cone-development/20170824")
#' ```

#' Libs
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(RColorBrewer))

#' Helper
source("~/Git/UPSCb/src/R/expressionSpecificityUtility.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)
pal <- brewer.pal(8,"Dark2")

#' # Data
#' ## Expression
counts <- read.csv("analysis/kallisto/raw-unormalised-gene-expression-tech-rep-combined_data.csv",row.names = 1)


#' ## Samples
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/biological-samples.csv",row.names=1)

#' ## VST
dds.date <- DESeqDataSetFromMatrix(counts,
                                   colData = samples[,c("ID","Sex","Date")],
                                   design=~Date)

dds.tissue <- DESeqDataSetFromMatrix(counts,
                                     colData = samples[,c("ID","Sex","Date")],
                                     design=~Sex)

vsd <- varianceStabilizingTransformation(dds.date,blind = FALSE)
vst.date <- assay(vsd)
vst.date <- vst.date - min(vst.date)
vst.date <- vst.date[featureSelect(vst.date,samples$Biol_repl_ID,exp=4),]

vsd <- varianceStabilizingTransformation(dds.tissue,blind = FALSE)
vst.tissue <- assay(vsd)
vst.tissue <- vst.tissue - min(vst.tissue)
vst.tissue <- vst.tissue[featureSelect(vst.tissue,samples$Biol_repl_ID,exp=5),]

#' # Expression specificity
es.date <- expressionSpecificity(vst.date,as.character(samples[,"Date"]),"local","complete")
es.tissue <- expressionSpecificity(vst.tissue,as.character(samples[,"Sex"]),"local","complete")

#' Sort by expression peak and tissue specificity
ep.date <- colnames(es.date[,-match(c("score","maxn","n"),
                                    colnames(es.date))])[apply(es.date[,-match(c("score","maxn","n"),
                                                                               colnames(es.date))] == es.date[,"maxn"],1,which)]
ep.tissue <- colnames(es.tissue[,-match(c("score","maxn","n"),
                                        colnames(es.tissue))])[apply(es.tissue[,-match(c("score","maxn","n"),
                                                                                       colnames(es.tissue))] == es.tissue[,"maxn"],1,which)]

ord.date <- order(ep.date,(1-es.date[,"score"]))
ord.tissue <- order(ep.tissue,(1-es.tissue[,"score"]))

#' Export
write.csv(cbind(es.date,ep.date)[ord,],file="vst-expression-data_date-aware_expression-specificity-order.csv")
write.csv(cbind(es.tissue,ep.tissue)[ord,],file="vst-expression-data_tissue-aware_expression-specificity-order.csv")

#' # Heatmap
s.ord <-  order(samples$Date)
heatmap.2(vst.date[ord.date,s.ord],Colv=FALSE,Rowv=FALSE,
          dendrogram="none",col=hpal,
          trace="none",scale="row",labRow = FALSE,
          ColSideColors = pal[as.integer(samples$Date)][s.ord],
          RowSideColors = pal[as.integer(factor(ep.date))][ord.date])

s.ord <- order(samples$Sex)
heatmap.2(vst.tissue[ord.tissue,s.ord],Colv=FALSE,Rowv=FALSE,
          dendrogram="none",col=hpal,
          trace="none",scale="row",labRow = FALSE,
          ColSideColors = pal[as.integer(samples$Sex)],
          RowSideColors = pal[as.integer(factor(ep.tissue))][ord.tissue])


## Vst for both tisssue and date together  (codes from Nico)
## Data
#Expression
counts <- read.csv("analysis/kallisto/raw-unormalised-gene-expression-tech-rep-combined_data.csv",row.names = 1)
colnames(counts) <- gsub("\\.","-",colnames(counts))
colnames(counts) <- sub("--","-",colnames(counts))

##Samples
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/biological-samples.csv",row.names=1)

samples.sel <- samples$Date %in% c("08-01","08-19","09-16") &
  samples$Sex %in% c("F","VL") & 
  samples$Type == "WT"

stopifnot(all(colnames(counts) == samples$ID))

#VST
dds <- DESeqDataSetFromMatrix(counts[,samples.sel],
                              colData = samples[samples.sel,
                                                c("ID","Sex","Date")],
                              design=~Date+Sex)

# VSD
vsd <- varianceStabilizingTransformation(dds,blind = FALSE)
vst.date.sex <- assay(vsd)
vst.date.sex <- vst.date.sex - min(vst.date.sex)

#Export
write.csv(vst.date.sex,file="analysis/vst-expression-data_date-sex-aware.csv")

# MODIFIED BY Nico (2018-01-16)
# Read in the IDS
ID <- scan("~/Git/UPSCb/projects/spruce-cone-development/doc/all-diff-genes-from-august-and-Sep.csv",what="character")

## Read shirin's gene ID file
ID <- scan("~/Git/UPSCb/projects/spruce-cone-development/src/R/all diff genes from august and Sep.csv",what="character")

## WARNING Working only on the DE genes of interest
s.vst<-t(scale(t(vst.date.sex[ID,])))
s.ord<-  order(colData(dds)$Date) 

# from this heatmap, we can see that the sex has the strongest effect and that the date also has an effect
# within a sex
heatmap.2(s.vst[,s.ord],Colv=FALSE,
          dendrogram="row",col=hpal,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D2")},
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(colData(dds)$Date)][s.ord])

##new codes by shirin to extract information from heatmap (22.01.2017)
## hclust and cutree to extract the information:

hc<-hclust(dist(s.vst[,s.ord]))
##Plot only dendogram
plot(hc)
## Cut tree to get sub-groups
ct<- cutree(hc, k=15)

## draw red rectangles to mark the subgroups

rect.hclust(hc, k=15)

# Expression specificity
# First by date
es.date <- expressionSpecificity(vst.date.sex[ID,],
                                 as.character(colData(dds)$Date),"local","complete")

#' Sort by expression peak and tissue specificity (new codes)
ep.date <- colnames(es.date[,-match(c("score","maxn","n"),
                                    colnames(es.date))])[apply(es.date[,-match(c("score","maxn","n"),
                                                                               colnames(es.date))] == es.date[,"maxn"],1,which)] 


ord.date<- order(ep.date,(1-es.date[,"score"])) 

heatmap.2(s.vst[ord.date,s.ord],
          Colv=FALSE,
          Rowv=FALSE,
          dendrogram="none",col=hpal,
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(colData(dds)$Date)][s.ord],
          RowSideColors = pal[as.integer(factor(ep.date[ord.date]))])

heatmap.2(s.vst[,s.ord],Colv=FALSE,
          dendrogram="row",col=hpal,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D2")},
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(colData(dds)$Date)][s.ord],
          RowSideColors = pal[as.integer(factor(ep.date))])


# Then by sex
es.sex <- expressionSpecificity(vst.date.sex[ID,],
                                as.character(colData(dds)$Sex),"local","complete")

#' Sort by expression peak and tissue specificity
ep.sex <- colnames(es.sex[,-match(c("score","maxn","n"),
                                  colnames(es.sex))])[apply(es.sex[,-match(c("score","maxn","n"),
                                                                           colnames(es.sex))] == es.sex[,"maxn"],1,which)] 


ord.sex<- order(ep.sex,(1-es.sex[,"score"])) 
s.ord <- order(colData(dds)$Sex) 

#' in the heatmap, the ribbon are sorted in alphabetical order of Sex (F then V)
heatmap.2(s.vst[ord.sex,s.ord],Colv=FALSE,Rowv=FALSE,
          dendrogram="none",col=hpal,
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(colData(dds)$Sex)][s.ord],
          RowSideColors = pal[as.integer(factor(ep.sex[ord.sex]))])

heatmap.2(s.vst[,s.ord],Colv=FALSE,
          dendrogram="row",col=hpal,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D2")},
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(colData(dds)$Sex)][s.ord],
          RowSideColors = pal[as.integer(factor(ep.sex))])


# Combine the biological replicates (take the median)
m.vst <- sapply(split.data.frame(t(vst.date.sex[ID,]),paste0(colData(dds)$Sex,colData(dds)$Date)),
                colMedians)
rownames(m.vst) <- ID
s.m.vst <- t(scale(t(m.vst)))

heatmap.2(s.m.vst,Colv=FALSE,
          dendrogram="row",col=hpal,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D2")},
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(factor(substr(colnames(s.m.vst),1,1)))],
          RowSideColors = pal[as.integer(factor(ep.sex))])

heatmap.2(s.m.vst,Colv=FALSE,
          dendrogram="row",col=hpal,
          distfun = pearson.dist,
          hclustfun = function(X){hclust(X,method="ward.D2")},
          trace="none",labRow = FALSE,
          ColSideColors = pal[as.integer(factor(substr(colnames(s.m.vst),1,nchar(colnames(s.m.vst)))))],
          RowSideColors = pal[as.integer(factor(ep.date))])

##
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

