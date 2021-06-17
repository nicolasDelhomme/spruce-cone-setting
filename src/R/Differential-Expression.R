#' ---
#' title: "Cone Development Differential Expression Analyses"
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
setwd("/mnt/picea/projects/spruce/onilsson/cone-development")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/onilsson/cone-development")
#' ```

#' Libs
#' ```{r drop loading warnings, echo=FALSE}
#' options(warn=-1)
#' ```
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(VennDiagram))
#' ```{r restore warnings, echo=FALSE}
#' options(warn=0)
#' ```

#' Helper files
suppressMessages(source("~/Git/UPSCb/src/R/densityPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))

#' Load saved data
load("analysis/DESeq2/deseq-sampling_time-tissue_type.rda")
vst <- read.csv("analysis/HTSeq/library-size-normalized_variance-stabilized_data.csv",
                row.names=1)

#' Setup graphics
pal=brewer.pal(8,"Dark2")
mar <- par("mar")

# ----
#' # Process

#' Differential Expression
#' This fails because the matrix is not full rank, 
#' as not all interaction exists
#' ```{r not full rank missing interaction, eval=FALSE}
#' dds <- DESeq(dds)
#' ```

#' This also fails because effects are confounded, e.g. P
#' was only sequenced on 05-12
#' ```{r confounded variable, eval=FALSE}
#' design(dds) <- ~Sampling + Sex
#' dds <- DESeq(dds)
#' ```

#' Plotting the information
#' Clearly the O, P and S samples are date confounded, so they 
#' need to be removed for the time being
barplot(sapply(split(colData(dds)$Sex,colData(dds)$Sampling),table),
        beside=TRUE,col=pal[1:nlevels(colData(dds)$Sex)])

legend("topleft",legend=levels(colData(dds)$Sex),
       fill=pal[1:nlevels(colData(dds)$Sex)])

#' Removing the O, P and S samples using their date
dds.s <- dds[,!colData(dds)$Sampling %in% c("05-12","06-11")]
vst <- vst[,!colData(dds)$Sampling %in% c("05-12","06-11")]

#' and adjusting the factor
dds.s$Sampling <- droplevels(dds.s$Sampling)
dds.s$Sex <- droplevels(dds.s$Sex)

#' The absence of V for the 09-12 date makes it impossible to
#' get a full rank interaction matrix. However, looking at the
#' biological QA heatmap and knowing that the meristem of the 
#' cones mature differently between sex, we can add a novel
#' variable
colData(dds.s)$Meristem <- factor(ifelse(colData(dds.s)$Sampling=="08-18",
                                  "1",
                                  ifelse(colData(dds.s)$Sampling 
                                         %in% c("09-12","10-08"),
                                         ifelse(colData(dds.s)$Sex=="M",
                                         "3","2"),
                                    "4")))

#' Furthermore, we are not interested in late meristem samples
dds.m <- dds.s[,colData(dds.s)$Meristem %in% c(1:3)]
vst <- vst[,colData(dds.s)$Meristem %in% c(1:3)]
dds.m$Sampling <- droplevels(dds.m$Sampling)
dds.m$Meristem <- droplevels(dds.m$Meristem)
design(dds.m) <- ~Meristem

#' Do the differential expression
dds.m <- DESeq(dds.m)

#' Dispersion Estimation
#' 
#' The dispersion estimation is adequate
plotDispEsts(dds.m)

#' Look at what contrasts are available
resultsNames(dds.m)

# ----
#' ## Meristem 3 (F|V at 09|10) vs 2 (M at 09|10)
res <- results(dds.m,contrast=c("Meristem","3","2"))

#' assume a 1% FDR
alpha=0.01

#' Plot the Median vs Average
#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
#' There are many genes that appear differentially expressed at a 1% FDR cutoff
#' The log2 fold-change range is relatively broad, with three extreme
#' values
plotMA(res,alpha)

#' Plot the log odds vs. log2 fold change
#' 
#' The volcano plot shows the same results as the MA plot; a
#' large number of genes show significant fold-changes
volcanoPlot(res,alpha=alpha)

#' Plot the adjusted p-value histogram
#' 
#' Which is almost evenly distributed, with an enrichment for 
#' lower p-values (more significant)
hist(res$padj,breaks=seq(0,1,.01))

#' Select genes below alpha
#' 
#' Note the cutoff selection on the log fold change is motivated by the 
#' Schurch et al. RNA, 2016 publication
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
sprintf("There are %s genes differentially expressed at a %s cutoff",
        sum(sel),alpha)

# plot(density(log2(res[sel,"baseMean"])))
# 
# sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5 & 
#   log2(res$baseMean) > 3
# sprintf("There are %s genes differentially expressed at a %s cutoff",
#         sum(sel),alpha)

#' There is more over-expression in Meristem 2 (-1) than in 3 (1)
pander(table(sign(res[sel,"log2FoldChange"])))

# #' Write them out
write(sub("\\.0$","",rownames(res)[sel]),
       file=file.path("analysis/DESeq2",
                      "Meristem3-vs-Meristem2_one-percent-FDR-cutoff_significant-genes.txt"))
 
write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange)==1]),
      file=file.path("analysis/DESeq2",
                     "m3-vs-m2-Meristem3-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))

write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange) == -1]),
      file=file.path("analysis/DESeq2",
                     "m3-vs-m2-Meristem2-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))

write.csv(as.data.frame(res),
          file=file.path("analysis/DESeq2",
                         "Meristem3_vs_meristem2.csv"))

m3.vs.m2 <- rownames(res)[sel]

#' #### Cluster the VST expression of the genotype genes
#' The color code for the column next to the row dendogram (genes)
#' is yellow for significant positive fold change and darkorange for
#' significant negative fold change
heatmap.2(as.matrix(vst[rownames(vst) %in% m3.vs.m2,]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "pink",
                               "blue"),
          labCol=paste(colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))

heatmap.2(as.matrix(vst[rownames(vst) %in% m3.vs.m2,]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "pink",
                               "blue"),
          labCol=paste(colData(dds.m)$Sampling,colData(dds.m)$Sample,
            colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))


# ----
#' ## Meristem 3 (F|V at 09|10) vs 1 (M|F|V at 08)
res <- results(dds.m,contrast=c("Meristem","3","1"))

#' assume a 1% FDR
alpha=0.01

#' There are  many genes that appear differentially expressed. 
#' The log2 fold-change range is relatively broad, with two extreme
#' positive values
plotMA(res,alpha)

#' Plot the log odds vs. log2 fold change
#' 
#' The volcano plot shows the same results as the MA plot; a
#' large number of genes show significant fold-changes
volcanoPlot(res,alpha=alpha)

#' Plot the adjusted p-value histogram
#' 
#' Which is almost evenly distributed, with an enrichment for 
#' lower p-values (most significant)
hist(res$padj,breaks=seq(0,1,.01),xlab="FDR")

#' Select genes below alpha
#' 
#' Note the 0.5 cutoff on the log fold change is motivated by the 
#' Schurch et al. RNA, 2016 publication
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
sprintf("There are %s genes differentially expressed at a %s cutoff",
        sum(sel),alpha)

#' Slightly skewed towards over-expression in meristem3 (1) than in meristem1 (-1)
pander(table(sign(res[sel,"log2FoldChange"])))

#' Write them out
write(sub("\\.0$","",rownames(res)[sel]),
      file=file.path("analysis/DESeq2",
                     "Meristem3-vs-Meristem1_one-percent-FDR-cutoff_significant-genes.txt"))

write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange)==1]),
      file=file.path("analysis/DESeq2",
                     "m3-vs-m1-Meristem1-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))

write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange) == -1]),
      file=file.path("analysis/DESeq2",
                     "m3-vs-m1-Meristem3-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))

write.csv(as.data.frame(res),
          file=file.path("analysis/DESeq2",
                         "Meristem3_vs_meristem1.csv"))

m3.vs.m1 <- rownames(res)[sel]

#' #### Cluster the VST expression of the genotype genes
#' The color code for the column next to the row dendogram (genes)
#' is yellow for significant positive fold change and dark orange for
#' significant negative fold change
heatmap.2(as.matrix(vst[rownames(vst) %in% m3.vs.m1,]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "pink",
                               "blue"),
          labCol=paste(colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))

# ----
#' ## Meristem 2 (M at 09|10) vs 1 (M|F|V at 08)
res <- results(dds.m,contrast=c("Meristem","2","1"))

#' assume a 1% FDR
alpha=0.01

#' There are not many genes that appear differentially expressed. 
#' The log2 fold-change range is relatively broad, with one extreme
#' value.
plotMA(res,alpha)

#' Plot the log odds vs. log2 fold change
#' 
#' The volcano plot shows the same results as the MA plot; a
#' large number of genes show large insignificant fold-changes 
#' in term of the interaction of pcp and final temperature
volcanoPlot(res,alpha=alpha)

#' Plot the adjusted p-value histogram
#' 
#' Which is almost evenly distributed, with an enrichment for 
#' higher p-values (less significant)
hist(res$padj,breaks=seq(0,1,.01),xlab="FDR")

#' Select genes below alpha
#' 
#' Note the 0.5 cutoff on the log fold change is motivated by the 
#' Schurch et al. RNA, 2016 publication
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
        
sprintf("There are %s genes differentially expressed at a %s cutoff",
        sum(sel),alpha)

#' Almost evenly distributed
pander(table(sign(res[sel,"log2FoldChange"])))

#' Write them out
write(sub("\\.0$","",rownames(res)[sel]),
      file=file.path("analysis/DESeq2",
                     "Meristem2-vs-Meristem1_one-percent-FDR-cutoff_significant-genes.txt"))

write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange)==1]),
      file=file.path("analysis/DESeq2",
                     "m2-vs-m1-Meristem2-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))

write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange) == -1]),
      file=file.path("analysis/DESeq2",
                     "m2-vs-m1-Meristem1-down-regulated_one-percent-FDR-cutoff_significant-genes.txt"))

write.csv(as.data.frame(res),
          file=file.path("analysis/DESeq2",
                         "Meristem2_vs_Meristem1.csv"))

m2.vs.m1 <- rownames(res)[sel]

#' #### Cluster the VST expression of the interaction genes
#' The color code for the column next to the row dendogram (genes)
#' is yellow for significant positive fold change and dark orange for
#' significant negative fold change
heatmap.2(as.matrix(vst[rownames(vst) %in% m2.vs.m1,]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "pink",
                               "blue"),
          labCol=paste(colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))

#' ### Check the overlap of the three sets
plot.new()
grid.draw(venn.diagram(list(
  m3.vs.m2=m3.vs.m2,
  m3.vs.m1=m3.vs.m1,
  m2.vs.m1=m2.vs.m1),
  filename=NULL,
  col=pal[1:3],
  category.names=c("m3.vs.m2","m3.vs.m1","m2.vs.m1")))

#' intersect m3.vs.m2 and m2.vs.m1
heatmap.2(as.matrix(vst[rownames(vst) %in% intersect(m3.vs.m2,m2.vs.m1),]),
          scale="row",labRow=NA,trace="none",
          labCol=paste(colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))


heatmap.2(as.matrix(vst[rownames(vst) %in% intersect(m3.vs.m1,m2.vs.m1),]),
          scale="row",labRow=NA,trace="none",
          labCol=paste(colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))

#' intersect m3.vs.m2 and m2.vs.m1 (data is not scaled, so absolute expression is shown)
heatmap.2(as.matrix(vst[rownames(vst) %in% intersect(m3.vs.m1,m2.vs.m1),]),
          labRow=NA,trace="none",
          labCol=paste(colData(dds.m)$Meristem,colData(dds.m)$Sex,sep="-"))

#' Lists of genes from different intersection in the venndiagram
#' number2 <- intersect(m3.vs.m1,intersect(m3.vs.m2,m2.vs.m1))
#' number1 <- setdiff(intersect(m3.vs.m2,m2.vs.m1),number2)
#' number3 <- setdiff(m2.vs.m1,union(m3.vs.m2,m3.vs.m1))
common.M3vsM2.M2vsM1.M3vsM1.gene.list <- intersect(m3.vs.m1,intersect(m3.vs.m2,m2.vs.m1))
common.M3vsM2.M2vsM1.minus.M3vsM1.gene.list <- setdiff(intersect(m3.vs.m2,m2.vs.m1),common.M3vsM2.M2vsM1.M3vsM1.gene.list)
write(common.M3vsM2.M2vsM1.minus.M3vsM1.gene.list, file="analysis/DESeq2/common_M3vsM2_M2vsM1_minus_M3vsM1_gene_list.txt")
common.M3vsM2.M3vsM1.minus.M2vsM1.gene.list <- setdiff(intersect(m3.vs.m2,m3.vs.m1),common.M3vsM2.M2vsM1.M3vsM1.gene.list)
write(common.M3vsM2.M3vsM1.minus.M2vsM1.gene.list, file="analysis/DESeq2/common_M3vsM2_M3vsM1_minus_M2vsM1_gene_list.txt")
common.M3vsM1.M2vsM1.minus.M3vsM2.gene.list <- setdiff(intersect(m3.vs.m1,m2.vs.m1),common.M3vsM2.M2vsM1.M3vsM1.gene.list)
write(common.M3vsM1.M2vsM1.minus.M3vsM2.gene.list, file="analysis/DESeq2/common_M3vsM1_M2vsM1_minus_M3vsM2_gene_list.txt")

#' # Conclusion
#'
#' 
#' ```{r empty,echo=FALSE,eval=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
