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
#' Select only the M and F sample from September
sel <- colData(dds)$Sex %in% c("M","F") & colData(dds)$Sampling == "09-12"
dds.s <- dds[,sel]
vst <- vst[,sel]

#' and adjusting the factor
dds.s$Sampling <- droplevels(dds.s$Sampling)
dds.s$Sex <- droplevels(dds.s$Sex)

#' adapt the design
design(dds.s) <- ~Sex

#' Do the differential expression
dds.s <- DESeq(dds.s)

#' Dispersion Estimation
#' 
#' The dispersion estimation is adequate
plotDispEsts(dds.s)

#' Look at what contrasts are available
resultsNames(dds.s)

# ----
#' ## Sex M vs. F
res <- results(dds.s,contrast=c("Sex","M","F"))

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

#' There is more over-expression in F (-1) than in M (1)
pander(table(sign(res[sel,"log2FoldChange"])))

# # #' Write them out
# write(sub("\\.0$","",rownames(res)[sel]),
#        file=file.path("analysis/DESeq2",
#                       "Meristem3-vs-Meristem2_one-percent-FDR-cutoff_significant-genes.txt"))
#  
# write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange)==1]),
#       file=file.path("analysis/DESeq2",
#                      "m3-vs-m2-Meristem2-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))
# 
# write(sub("\\.0$","",rownames(res)[sel & sign(res$log2FoldChange) == -1]),
#       file=file.path("analysis/DESeq2",
#                      "m3-vs-m2-Meristem3-up-regulated_one-percent-FDR-cutoff_significant-genes.txt"))
# 
# write.csv(as.data.frame(res),
#           file=file.path("analysis/DESeq2",
#                          "Meristem3_vs_meristem2.csv"))

m.vs.f <- rownames(res)[sel]

#' #### Cluster the VST expression of the genotype genes
#' The color code for the column next to the row dendogram (genes)
#' is yellow for significant positive fold change and darkorange for
#' significant negative fold change
heatmap.2(as.matrix(vst[rownames(vst) %in% m.vs.f,]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "pink",
                               "blue"),
          labCol=paste(colData(dds.s)$Sample,
                       colData(dds.s)$Sex,sep="-"))


#' # Conclusion
#'
#' 
#' ```{r empty,echo=FALSE,eval=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
