#' ---
#' title: "Spruce Cone development all samples Differential Expression"
#' author: "Nicolas Delhomme and Veronika Nordal"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' # Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/u2015029/20170824")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/u2015029/20170824")
#' ```

#' Libs
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(Glimma))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))

#' Helper files
suppressMessages(source("~/Git/UPSCb/src/R/densityPlot.R"))
#suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plot.multidensity.R"))

#' Load saved data
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/biological-samples.csv",row.names = 1)

#' Annotation
annot <- read.delim("/mnt/picea/storage/reference/Picea-abies/v1.1/annotation/Pabies-Pfam-description-annotation.tsv")
colnames(annot) <- c("GeneID","Confidence","PfamDomains")

#' Gene raw expression
counts <- read.csv("analysis/kallisto/raw-unormalised-gene-expression-tech-rep-combined_data.csv",
               row.names=1)
colnames(counts) <- gsub("\\.+","-",colnames(counts))

#' Setup graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
#' 
#' # Process
#' Re-order the samples
stopifnot(all(colnames(counts) == samples$ID))

#' ## DE at the gene level
#' A generic function to do the DE
#' default to a 1% FDR and abs(lfc) >=0.5

".de" <- function(counts,samples,design,nam,res.nam=NULL,alpha=0.01,lfc=0.5){
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = design)
  
  # Variance stabilisation
  vst <- varianceStabilizingTransformation(dds,blind=FALSE)
  vst <- assay(vsd)
  vst <- vst - min(vst)
  
  # Validate the VST 
  meanSdPlot(vst[rowSums(vst)>0,])
  
  # Differential Expression
  suppressMessages(dds <- DESeq(dds))
  
  # Dispersion Estimation
  #
  # The dispersion estimation is adequate
  plotDispEsts(dds)
  
  # Obtain the results 
  # resultsNames(dds)
  if(is.null(res.nam)){
    res <- results(dds)
  } else {
    res <- results(dds,name=res.nam)
  }
  
  
  # Plot the Median vs Average
  # There are many genes that appear differentially expressed at a 1% FDR cutoff
  # The log2 fold-change range is relatively broad, with three extreme
  # values
  plotMA(res,alpha)
  
  # Plot the log10 odds (i.e. -log10 FDR) vs. log2 fold change
  #
  # The volcano plot shows the same results as the MA plot; a
  # large number of genes show significant fold-changes
  volcanoPlot(res,alpha=alpha,lfc = lfc)
  
  # take a look at some results
  # res[res$padj < 1e-150 & !is.na(res$padj),]
  
  
  # Plot the adjusted p-value histogram
  #
  # Which is almost evenly distributed, with an enrichment for
  # lower p-values (more significant)
  hist(res$padj,breaks=seq(0,1,.01))
  
  # Select genes below alpha
  #
  # Note the 0.5 cutoff on the log fold change is motivated by the
  # Schurch et al. RNA, 2016 publication
  sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
  message(sprintf("There are %s genes differentially expressed at a %s cutoff",
          sum(sel),alpha))
  
  # Equally distributed between F (-1) and M (1)
  pander(table(sign(res[sel,"log2FoldChange"])))
  
  # Order the annotations
  annot <- annot[match(sub("\\.1$","",rownames(res)),annot$GeneID),]
  
  # Write them out
  outdir=file.path("analysis/DESeq2",nam)
  dir.create(outdir,
             showWarnings=FALSE,recursive=TRUE)
  write.csv(cbind(res,annot),
            file=file.path(outdir,paste0(nam,".csv")))
  
  # write the de selected genes
  de.sel <- rownames(res)[sel]
  write(sub("\\.1$","",de.sel),file=file.path(outdir,paste0(nam,"-DE-geneIDs.txt")))
  
  # #### Cluster the VST expression of the genotype genes
  # The color code for the column next to the row dendogram (genes)
  # is red for significant positive fold change and blue for
  # significant negative fold change
  heatmap.2(as.matrix(vst[rownames(vst) %in% de.sel,]),
            scale="row",labRow=NA,trace="none",
            RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                                 "blue",
                                 "red"),
            col=hpal,cexCol = 0.9, margins = c(6.1,0.1),
            labCol=samples[s.sel,"ID"])
  
  # #### Interactive MDS
  dir.create(file.path(outdir,nam),showWarnings = FALSE)
  res.df <- as.data.frame(res)
  res.df$log10MeanNormCount <- log10(res.df$baseMean)
  idx <- rowSums(counts(dds)) > 0
  res.df <- res.df[idx,]
  res.df$padj[is.na(res.df$padj)] <- 1
  rownames(res.df) <- sub("\\.1$","",rownames(res.df))
  glMDPlot(res.df,
           xval = "log10MeanNormCount",
           yval = "log2FoldChange",
           counts = counts(dds)[idx,],
           anno = annot[idx,],
           groups = as.character(dds$Sex),
           samples = colnames(dds),
           status = sign(res[idx,"log2FoldChange"]) * sel[idx],
           display.columns = c("GeneID", "Confidence","PfamDomains","padj","baseMean"),
           path = outdir,
           folder = nam,
           launch=FALSE)
  
  return(res)
}

#' ### Male vs. Female, 08-01
s.sel <- which(samples$Date == "08-01" & samples$Sex %in% c("F","M") & samples$Type=="WT")
res <- .de(counts[,s.sel],samples[s.sel,"Sex",drop=FALSE],~Sex,nam="M-vs-F_Aug-01")

#' ### Male vs. Veg. lead., 08-01

s.sel <- which(samples$Date == "08-01" & samples$Sex %in% c("M","VL") & samples$Type=="WT")

dds <- DESeqDataSetFromMatrix(
  countData = counts[,s.sel],
  colData = samples[s.sel,"Sex",drop=FALSE],
  design = ~Sex)

#' Variance stabilisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Validate the VST 
meanSdPlot(vst[rowSums(vst)>0,])

#' Differential Expression
suppressMessages(dds <- DESeq(dds))

#' Dispersion Estimation
#'
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Obtain the results 
# resultsNames(dds)
res <- results(dds)

#' assume a 1% FDR and limit to lfc >=0.5
alpha=0.01
lfc=0.5

#' Plot the Median vs Average
#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
plotMA(res,alpha)

#' Plot the log10 odds (i.e. -log10 FDR) vs. log2 fold change
#'
volcanoPlot(res,alpha=alpha,lfc = lfc)

# take a look at some results
# res[res$padj < 1e-150 & !is.na(res$padj),]


#' Plot the adjusted p-value histogram
hist(res$padj,breaks=seq(0,1,.01))

#' Select genes below alpha
#'
#' Note the 0.5 cutoff on the log fold change is motivated by the
#' Schurch et al. RNA, 2016 publication
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
sprintf("There are %s genes differentially expressed at a %s cutoff",
        sum(sel),alpha)

#' M (-1) and VL (1) (is this right??)
pander(table(sign(res[sel,"log2FoldChange"])))

#' Order the annotations
annot <- annot[match(sub("\\.1$","",rownames(res)),annot$GeneID),]

#' Write them out
nam="M-vs-VL_Aug-01"
outdir=file.path("analysis/DESeq2",nam)
dir.create(outdir,
           showWarnings=FALSE,recursive=TRUE)
write.csv(cbind(res,annot),
          file=file.path(outdir,paste0(nam,".csv")))

# write the de selected genes
de.sel <- rownames(res)[sel]
write(sub("\\.1$","",de.sel),file=file.path(outdir,paste0(nam,"-DE-geneIDs.txt")))

#' #### Cluster the VST expression of the genotype genes
#' The color code for the column next to the row dendogram (genes)
#' is yellow for significant positive fold change and darkorange for
#' significant negative fold change
heatmap.2(as.matrix(vst[rownames(vst) %in% de.sel,]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "blue",
                               "red"),
          col=hpal,cexCol = 0.9, margins = c(6.1,0.1),
          labCol=samples[s.sel,"ID"])

#' #### Interactive MDS
dir.create(file.path(outdir,"report"),showWarnings = FALSE)
res.df <- as.data.frame(res)
res.df$log10MeanNormCount <- log10(res.df$baseMean)
idx <- rowSums(counts(dds)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
rownames(res.df) <- sub("\\.1$","",rownames(res.df))
glMDPlot(res.df,
         xval = "log10MeanNormCount",
         yval = "log2FoldChange",
         counts = counts(dds)[idx,],
         anno = annot[idx,],
         groups = as.character(dds$Sex),
         samples = as.character(samples[s.sel,"ID"]),
         status = sign(res[idx,"log2FoldChange"]) * sel[idx],
         display.columns = c("GeneID", "Confidence","PfamDomains","padj","baseMean"),
         path = outdir,
         folder = nam,
         launch=FALSE)

#' ### Veg. lead. 08-19 vs. Female, 08-01
s1.sel <- which(samples$Date == "08-01" & samples$Sex %in% c("F") & samples$Type=="WT")
s2.sel <- which(samples$Date == "08-19" & samples$Sex %in% c("VL") & samples$Type=="WT")
s.sel <- c(s1.sel,s2.sel)

dds <- DESeqDataSetFromMatrix(
  countData = counts[,s.sel],
  colData = samples[s.sel,"Sex",drop=FALSE],
  design = ~Sex)

#' Variance stabilisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)

#' Validate the VST 
meanSdPlot(vst[rowSums(vst)>0,])

#' Differential Expression
suppressMessages(dds <- DESeq(dds))

#' Dispersion Estimation
#'
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Obtain the results 
# resultsNames(dds)
res <- results(dds)

#' assume a 1% FDR and limit to lfc >=0.5
alpha=0.01
lfc=0.5

#' Plot the Median vs Average
#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
plotMA(res,alpha)

#' Plot the log10 odds (i.e. -log10 FDR) vs. log2 fold change
#'
volcanoPlot(res,alpha=alpha,lfc = lfc)

# take a look at some results
# res[res$padj < 1e-150 & !is.na(res$padj),]


#' Plot the adjusted p-value histogram
#'
hist(res$padj,breaks=seq(0,1,.01))

#' Select genes below alpha
#'
#' Note the 0.5 cutoff on the log fold change is motivated by the
#' Schurch et al. RNA, 2016 publication
sel <- res$padj<alpha & !is.na(res$padj) & abs(res$log2FoldChange) >= 0.5
sprintf("There are %s genes differentially expressed at a %s cutoff",
        sum(sel),alpha)

#' F 08-01 (-1) and VL 08-19 (1) ???
pander(table(sign(res[sel,"log2FoldChange"])))

#' Order the annotations
annot <- annot[match(sub("\\.1$","",rownames(res)),annot$GeneID),]

#' Write them out
nam="F_Aug-01-vs-VL_Aug-19"
outdir=file.path("analysis/DESeq2",nam)
dir.create(outdir,
           showWarnings=FALSE,recursive=TRUE)
write.csv(cbind(res,annot),
          file=file.path(outdir,paste0(nam,".csv")))

# write the de selected genes
de.sel <- rownames(res)[sel]
write(sub("\\.1$","",de.sel),file=file.path(outdir,paste0(nam,"-DE-geneIDs.txt")))

#' #### Cluster the VST expression of the genotype genes
#' The color code for the column next to the row dendogram (genes)
#' is yellow for significant positive fold change and darkorange for
#' significant negative fold change
heatmap.2(as.matrix(vst[rownames(vst) %in% de.sel,s.sel]),
          scale="row",labRow=NA,trace="none",
          RowSideColors=ifelse(sign(res[sel,"log2FoldChange"])==-1,
                               "blue",
                               "red"),
          col=hpal,cexCol = 0.9, margins = c(6.1,0.1),
          labCol=samples[s.sel,"ID"])

#' #### Interactive MDS
dir.create(file.path(outdir,"report"),showWarnings = FALSE)
res.df <- as.data.frame(res)
res.df$log10MeanNormCount <- log10(res.df$baseMean)
idx <- rowSums(counts(dds)) > 0
res.df <- res.df[idx,]
res.df$padj[is.na(res.df$padj)] <- 1
rownames(res.df) <- sub("\\.1$","",rownames(res.df))
glMDPlot(res.df,
         xval = "log10MeanNormCount",
         yval = "log2FoldChange",
         counts = counts(dds)[idx,],
         anno = annot[idx,],
         groups = as.character(dds$Sex),
         samples = as.character(samples[s.sel,"ID"]),
         status = sign(res[idx,"log2FoldChange"]) * sel[idx],
         display.columns = c("GeneID", "Confidence","PfamDomains","padj","baseMean"),
         path = outdir,
         folder = nam,
         launch=FALSE)

#' ### Veg. lead. 08-19 vs. Female, 08-01
s.sel <- which(samples$Date %in% c("08-01","08-19") & samples$Sex %in% c("F","VL") & samples$Type=="WT")
res <- .de(counts[,s.sel],samples[s.sel,c("Date","Sex")],~Date*Sex,nam="M-vs-F_Aug-01",res.nam="Date08.19.SexVL")

# dds <- DESeqDataSetFromMatrix(
#   countData = counts[,s.sel],
#   colData = samples[s.sel,c("Date","Sex")],
#   design = ~Date*Sex)
# # we look at the covariates: Date, Sex and Date:Sex (interaction term)
# 
# dds <- DESeq(dds)
# 
# resultsNames(dds)

#' ### Veg. lead. 10-25 vs. acrocona transition shoot, 10-18
# s.sel <- which(samples$Date %in% c("10-25","10-18") & samples$Sex %in% c("TA","VL"))
# res <- .de(counts[,s.sel],samples[s.sel,c("Date","Sex")],~Date*Sex,nam="VL_Oct-25-vs-TA_Oct-18",res.nam="Date08.19.SexVL")
# res <- .de(counts[,s.sel],samples[s.sel,"Sex",drop=FALSE],~Sex,nam="M-vs-F_Aug-01")


#' ```{r empty,echo=FALSE,eval=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
