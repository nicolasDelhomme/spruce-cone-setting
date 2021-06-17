#' ---
#' title: "Spruce Cone development plots"
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
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(vsn))

#' Helpers
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Load saved data
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/biological-samples.csv",row.names = 1)

#' Gene raw expression
counts <- read.csv("analysis/kallisto/raw-unormalised-gene-expression-tech-rep-combined_data.csv",
               row.names=1)
colnames(counts) <- gsub("\\.+","-",colnames(counts))

#' Setup graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
rainbow17 <- rainbow(17)
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' # Process
#' Re-order the samples
stopifnot(all(colnames(counts) == samples$ID))

#' ## DE at the gene level
".plot" <- function(counts,samples,design,fact,pal=pal,legend.x="topright",heatmap=TRUE){
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = samples,
    design = design)
  
  # Variance stabilisation
  vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
  vst <- assay(vsd)
  vst <- vst - min(vst)
  
  # Validate the VST 
  meanSdPlot(vst[rowSums(vst)>0,])
  
  # PCA
  pc <- prcomp(t(vst))
  percent <- round(summary(pc)$importance[2,]*100)
  
  # ### 3 first dimensions
  mar=c(5.1,4.1,4.1,2.1)
  scatterplot3d(pc$x[,1],
                pc$x[,2],
                pc$x[,3],
                xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
                ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
                zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
                color=pal[as.integer(fact)],
                pch=19)
  legend(legend.x,pch=19,
         col=pal[1:nlevels(fact)],
         legend=levels(fact))
  par(mar=mar)
  
  # 1st and 2nd dim
  plot(pc$x[,1],
       pc$x[,2],
       xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
       ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
       col=pal[as.integer(fact)],
       pch=19,
       main="Principal Component Analysis",sub="variance stabilized counts")
  
  legend(legend.x,bty="n",col=pal,levels(fact),pch=19)
  
  # 2nd and 3rd dim
  plot(pc$x[,2],
       pc$x[,3],
       xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
       ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
       col=pal[as.integer(fact)],
       pch=19,
       main="Principal Component Analysis",sub="variance stabilized counts")
  
  legend(legend.x,bty="n",col=pal,levels(fact),pch=19)
  
  # Saturate the expression
  vst.sat <- vst[sel,]
  vst.sat[vst.sat<3] <- 3
  vst.sat[vst.sat>8] <- 8
  

  # heatmap
  if(heatmap){
    sel <- featureSelect(vst,samples$Biol_repl_ID, exp = 5,nrep = 2)
    heatmap.2(vst[sel,],trace = "none",labRow = FALSE, 
              col=hpal, ColSideColors = pal[as.integer(fact)],
              scale="row")
    
    heatmap.2(vst.sat,trace = "none",labRow = FALSE, 
              col=hpal, ColSideColors = pal[as.integer(fact)])
  }
  
  # Sample Hierarchical clustering
  hc <- hclust(dist(t(vst.sat)))
  plot(hc, labels=samples$ID,
       main = "Hierarchical clustering",cex=0.7)
  
}

#' ## 
# select the desired samples
sel <- which(samples$Date %in% c("08-01","08-18","08-19","09-16") & samples$Type == "WT")

# First argument does not change
# Second argument, "ID","Biol_repl_ID" needs to be there, add any variable that needs to be in the design
# Third argument (design) needs to be in the second argument
# Fourth is the factor to color the points
# Fifth is the palette (if <=8, pal, otherwise pal12)
.plot(counts[,sel],
      samples[sel,c("ID","Biol_repl_ID","Date")],
      ~Date,
      fact=factor(as.character(samples[sel,]$Biol_repl_ID)),
      pal=pal12)

#' ## 
# select the desired samples
s1.sel <- which(samples$Date %in% c("08-01","08-19","10-25") & 
                  samples$Sex %in% c("F","VL") & samples$Type == "WT")
s2.sel <- which(samples$Type=="A")
sel <- c(s1.sel,s2.sel)

# First argument does not change
# Second argument, "ID","Biol_repl_ID" needs to be there, add any variable that needs to be in the design
# Third argument (design) needs to be in the second argument
# Fourth is the factor to color the points
# Fifth is the palette (if <=8, pal, otherwise pal12)
.plot(counts[,sel],
      samples[sel,c("ID","Biol_repl_ID","Date")],
      ~Date,
      fact=factor(as.character(samples[sel,]$Biol_repl_ID)),
      pal=pal)

.plot(counts[,sel],
      samples[sel,c("ID","Biol_repl_ID","Date")],
      ~Date,
      fact=factor(as.character(samples[sel,]$Biol_repl_ID)),
      pal=pal,legend.x="top",heatmap = FALSE)

.plot(counts[,sel],
      samples[sel,c("ID","Biol_repl_ID","Date")],
      ~Date,
      fact=factor(as.character(samples[sel,]$Biol_repl_ID)),
      pal=pal,legend.x="topleft",heatmap = FALSE)

# select the desired samples
sel <- which(samples$Date %in% c("08-18","09-12","10-08","04-10","04-29","05-12","06-11"))

# First argument does not change
# Second argument, "ID","Biol_repl_ID" needs to be there, add any variable that needs to be in the design
# Third argument (design) needs to be in the second argument
# Fourth is the factor to color the points
# Fifth is the palette (if <=8, pal, otherwise pal12)
.plot(counts[,sel],
      samples[sel,c("ID","Biol_repl_ID","Date")],
      ~Date,
      fact=factor(as.character(samples[sel,]$Biol_repl_ID)),
      pal=rainbow17)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
