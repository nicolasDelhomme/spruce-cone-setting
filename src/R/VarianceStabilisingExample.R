#' ---
#' title: "Variance Stabilising Example"
#' author: "Nicolas Delhomme"
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

#' # Data
#' ## Expression
counts <- read.csv("analysis/kallisto/raw-unormalised-gene-expression-tech-rep-combined_data.csv",row.names = 1)
colnames(counts) <- gsub("\\.","-",colnames(counts))
colnames(counts) <- sub("--","-",colnames(counts))

#' ## Samples
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/biological-samples.csv",row.names=1)

samples.sel <- samples$Date %in% c("08-01","08-19","09-16") &
  samples$Sex %in% c("F","VL") & 
  samples$Type == "WT"

stopifnot(all(colnames(counts) == samples$ID))

#' ## VST
dds <- DESeqDataSetFromMatrix(counts[,samples.sel],
                              colData = samples[samples.sel,
                                                c("ID","Sex","Date")],
                              design=~Date+Sex)

vsd <- varianceStabilizingTransformation(dds,blind = FALSE)
vst.date.sex <- assay(vsd)
vst.date.sex <- vst.date.sex - min(vst.date.sex)

#' Export
write.csv(vst.date.sex,file="analysis/vst-expression-data_date-sex-aware.csv")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```

