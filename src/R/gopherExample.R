#' ---
#' title: "Gopher example"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Source the helper
source("~/Git/UPSCb/src/R/gopher.R")

## New gopher source
source("~delhomme/Git/UPSCb/src/R/gopher.R")


#' # Run
#' 
#' You can check the source code in "~/Git/UPSCb/src/R/gopher.R"
#' 
#' genes is the list of gene of interest
#' 
#' background is the population in which to look for an enrichment, typically
#' every expressed genes in an experiment
#' 

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

# TODO - HERE you probably want to save the vst objects as an rda file, i.e. use the save function
# so that you can instead of creating them every-time, just load the existing object using the load function

# TODO - HERE YOU do not do any selection on the samples, so you are using all samples, not just the one
# you are interested in. You should do first the selection, then the vst and then the filtering

#' task can be one of 'all', 'go', 'kegg', 'pfam'
#' HERE the first problem is that your IDs are transcript IDs (the .1 at the end). For enrichment, you need
#' gene IDs, so they need to be trimmed of the last `.1`
ID <- scan("~sakhter/Git/UPSCb/projects/spruce-cone-development/doc/GeneID_cluster_F_01_08.csv",what="character")
ID <- sub("\\.\\d+$","",ID)


#' HERE this is to devise an adequate noise cutoff
sels <- sapply(1:10,function(i){
  featureSelect(vst.tissue,samples$Biol_repl_ID,exp=i)
})

plot(sels)

#' FIXME - This won't work because sels is a list. 
#' 
#' population <-rownames(vst)[sels]
#' 
#' From the plot above, select a cutoff that gives you
#' between 20000 and 30000 genes, i.e. replace the X below
sel <- sels[[X]]
population <-rownames(vst)[sel]
population <- sub("\\.\\d+$","",population)

# genes=your list of DEG genes
# background=expressed genes in your samples (look in the Heatmap.R)
# task=go (pfam or kegg) c("go","pfam","kegg")

#' HERE the url was incorrect, it needs to be "pabies" litteraly.
#' NOW I do get enrichments
enrichment <- gopher(genes=ID,port = 12000,task = "go",background = population,url = "pabies")

#' # Results
#' 
#' mpat, mt, napt and nt are the value used for the computation of the 
#' parent-child relationship. 
#' 
#' padj is the Benjamini-Hochberg multiple correction
#' of the pval. 
#' 
#' namespace is the "root": Biological Process, Cellular Component
#' or Molecular function. 
#' 
#' id is the GO ID and def is the definition of the GO ID
#' 
str(test$GO)
str(enrichment$GO)

# write it to disk
write.table(test$GO,file="GeneID_cluster_F_01_08_GO")

#' ###################
#' New code - 20180228
#' ###################
#' This you need to do only once per script
#' 
#' Get the genes associated with a GO ID
gene.go.mapping <- read.delim("/mnt/picea/storage/reference/Picea-abies/v1.0/gopher/gene_to_go.tsv",
                              col.names=c("gene","GO"),header=FALSE)
#' And reverse the mapping GO -> gene
go.gene.mapping <- strsplit(as.character(gene.go.mapping$GO),"\\|")
go.gene.mapping <- data.frame(GO=unlist(go.gene.mapping),
                              gene=rep(as.character(gene.go.mapping$gene),elementNROWS(go.gene.mapping)))

#' This you have to do for every of your list of interest
#' 
#' Subset to the gene in you list of interest
go.gene.subset <- go.gene.mapping[go.gene.mapping$gene %in% ID,]

#' And then lookup a GO ID of interest, I just took one at random.
go.gene.subset[go.gene.subset$GO %in% "GO:0006355","gene"]
#' ###################


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
