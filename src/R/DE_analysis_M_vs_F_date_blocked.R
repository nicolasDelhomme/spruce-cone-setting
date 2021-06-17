#' ---
#' title: "Spruce cone development Differential Expression Analysis"
#' author: "Veronika Nordal & Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/onilsson/cone-development/")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/onilsson/cone-development")
#' ```

library(DESeq2)
library(VennDiagram)
library(RColorBrewer)
source("~/Git/UPSCb/src/R/volcanoPlot.R")
source("~/Git/UPSCb/src/R/plotMA.R")

pal <- brewer.pal(8,"Dark2")

#' # Differential expression
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-cone-development/doc/samples.csv",
                    row.names=1,as.is=TRUE)

#' Load the count table
count.table <- read.csv("analysis/HTSeq/raw-unormalised_tech-rep-combined_data.csv",
                        row.names=1,as.is=TRUE)

#' Block date effect
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=colnames(count.table),
                       sex=samples$Sex,
                       date=samples$Sampling,
                       rep=samples$Sample),
  design = ~date+sex)


dds <- DESeq(dds)
plotDispEsts(dds)

#' Compare male with female (date blocked)
res <- results(dds,contrast=c("sex","M","F"))
alpha=0.01
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=1e-6)
hist(res$padj,breaks=seq(0,1,.01))

#' Number of DE genes, cutoff 1e-6 
sum(res$padj<1e-6,na.rm=TRUE)
#' Number of down- and upregulated genes
table(sign(res$log2FoldChange[res$padj<1e-6]))

#' Compare different cutoffs for padj and log2FoldChange
plot(density(res$log2FoldChange[res$padj<1e-6 & !is.na(res$padj)]))
abline(v=c(-1,1))
plot(density(res$log2FoldChange[res$padj<1e-3 & !is.na(res$padj)]))
abline(v=c(-1,1))
plot(density(res$log2FoldChange[res$padj<1e-3 & !is.na(res$padj) & abs(res$log2FoldChange) > 1]))
hist(res$log2FoldChange[res$padj<1e-3 & !is.na(res$padj) & abs(res$log2FoldChange) > 1])
hist(res$log2FoldChange[res$padj<1e-3 & !is.na(res$padj) & abs(res$log2FoldChange) > 1],breaks=seq(-10,10,0.5))
abline(v=c(-2,2))
hist(res$log2FoldChange[res$padj<1e-3 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],breaks=seq(-10,10,0.5))
hist(res$log2FoldChange[res$padj<1e-6 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],breaks=seq(-10,10,0.5))
hist(res$log2FoldChange[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],breaks=seq(-10,10,0.5))

#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)

sex.effect.date.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])

#' Write files
#' List with P-value and fold change per gene
write.csv(res,file="analysis/DESeq2/differential-expression_M-vs-F_date-blocked.csv")
write.csv2(res,file="analysis/DESeq2/differential-expression_M-vs-F_date-blocked_semi-colon.csv")
#' List of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_M-vs-F_date-blocked_gene-list.txt")

# to change the design e.g for sex
# design(dds) <- ~sex


#' F vs V date blocked
# dds <- DESeqDataSetFromMatrix(
#   countData = count.table,
#   colData = data.frame(condition=colnames(count.table),
#                        sex=samples$Sex,
#                        date=samples$Sampling,
#                        rep=samples$Sample),
#   design = ~date+sex)
# dds <- DESeq(dds)
# plotDispEsts(dds)
resultsNames(dds)
res <- results(dds,contrast=c("sex","F","V"))

#' Compare different cutoffs
alpha=0.01
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
volcanoPlot(res,alpha=0.001)
volcanoPlot(res,alpha=1e-6)
hist(res$padj,breaks=seq(0,1,.01))
sum(res$padj<1e-6,na.rm=TRUE)
sum(res$padj<1e-3,na.rm=TRUE)
table(sign(res$log2FoldChange[res$padj<1e-6]))
table(sign(res$log2FoldChange[res$padj<1e-3]))
plot(density(res$log2FoldChange[res$padj<1e-3 & !is.na(res$padj)]))
abline(v=c(-1,1))
sum(res$padj<1e-2,na.rm=TRUE)
table(sign(res$log2FoldChange[res$padj<1e-2]))
#' Number of DE genes padj cutoff 1e-2 and log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_F-vs-V_date-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_F-vs-V_date-blocked_gene-list.txt")
FvsV.effect.date.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])

#' Compare dates
res <- results(dds,contrast=c("date","10-08","09-12"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
date.effect.sex.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_10-08-vs-09-12_sex-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_10-08-vs-09-12_sex-blocked_gene-list.txt")

#' The complete model
#' We include the genotype.
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=colnames(count.table),
                       sex=samples$Sex,
                       date=samples$Sampling,
                       rep=samples$Sample,
                       genotype=ifelse(samples$Sample %in% c(1,2,4),"G1","G2")),
  design = ~sex+genotype+date)
dds <- DESeq(dds)

#' Compare dates again (sex and genotype blocked)
#' Contrast is 10-08 vs 09-12, so expression is (log2) 10-08 - 09-12 (linear 10-08/09-12).
#' Positive logFC indicate higher expression in 10-08.
res <- results(dds,contrast=c("date","10-08","09-12"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
hist(res$log2FoldChange[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],breaks=seq(-10,10,0.5))
hist(res$log2FoldChange[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1],breaks=seq(-10,10,0.5))
date.effect.genotype.sex.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1)
#' Venn diagram to compare the lists of date DE genes with and without genotype blocked
plot.new()
grid.draw(venn.diagram(list(date.effect.sex.blocked,date.effect.genotype.sex.blocked),
  filename=NULL,
  category.names=c("with G","blocked G"),
  col=pal[1:2]
))

#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_10-08-vs-09-12_sex_and_genotype-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_10-08-vs-09-12_sex_and_genotype-blocked_gene-list.txt")

#' Compare male with female again (date and genotype blocked)
res <- results(dds,contrast=c("sex","M","F"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
sex.effect.genotype.date.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])
sex.effect.genotype.date.blocked.2 <-rownames(res[!is.na(res$padj) & res$padj < 1e-2 & abs(res$log2FoldChange) >= 1.5,])
#' Venn diagram to compare the lists of date DE genes with and without genotype blocked
plot.new()
grid.draw(venn.diagram(list(sex.effect.date.blocked,sex.effect.genotype.date.blocked),
                       filename=NULL,
                       category.names=c("with G","blocked G"),
                       col=pal[1:2]
))
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_M-vs-F_date_and_genotype-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_M-vs-F_date_and_genotype-blocked_gene-list.txt")

#' Compare female vs veg again (date and genotype blocked)
res <- results(dds,contrast=c("sex","F","V"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
FvsV.effect.genotype.date.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])
FvsV.effect.genotype.date.blocked.2 <-rownames(res[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5,])
#' Venn diagram to compare the lists of DE genes with and without genotype blocked
plot.new()
grid.draw(venn.diagram(list(FvsV.effect.date.blocked,FvsV.effect.genotype.date.blocked),
                       filename=NULL,
                       category.names=c("with G","blocked G"),
                       col=pal[1:2]))
                       
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_F-vs-V_date_and_genotype-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
  file="analysis/DESeq2/differential-expression_F-vs-V_date_and_genotype-blocked_gene-list.txt")                      

#' Compare male vs veg (date and genotype blocked)
res <- results(dds,contrast=c("sex","M","V"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
MvsV.effect.genotype.date.blocked <- rownames(res[!is.na(res$padj) & res$padj < alpha,])
MvsV.effect.genotype.date.blocked.2 <-rownames(res[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5,])
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_M-vs-V_date_and_genotype-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_M-vs-V_date_and_genotype-blocked_gene-list.txt")

#' Compare genotype G1 vs G2 (date and sex blocked)
res <- results(dds,contrast=c("genotype","G1","G2"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_G1-vs-G2_date_and_sex-blocked_semi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5],
      file="analysis/DESeq2/differential-expression_G1-vs-G2_date_and_sex-blocked_gene-list.txt")

#' Venn diagram to compare the lists of DE genes comparing M vs F, F vs V and M vs V
plot.new()
grid.draw(venn.diagram(list(sex.effect.genotype.date.blocked,FvsV.effect.genotype.date.blocked,MvsV.effect.genotype.date.blocked),
                       filename=NULL,
                       category.names=c("M vs F","F vs V","M vs V"),
                       col=pal[1:3]))
plot.new()
grid.draw(venn.diagram(list(sex.effect.genotype.date.blocked.2,FvsV.effect.genotype.date.blocked.2,MvsV.effect.genotype.date.blocked.2),
                       filename=NULL,
                       category.names=c("M vs F","F vs V","M vs V"),
                       col=pal[1:3]))
#' List of DE genes common between the three comparisons
common.MvsF.FvsV.MvsV.gene.list <- intersect(sex.effect.genotype.date.blocked.2,intersect(FvsV.effect.genotype.date.blocked.2,MvsV.effect.genotype.date.blocked.2))
write(common.MvsF.FvsV.MvsV.gene.list,
      file="analysis/DESeq2/common_MvsF_FvsV_MvsV_date_and_genotype-blocked_gene-list.txt")
#' List of DE genes present only in MvsF comparison
MvsF.minus.FvsV.gene.list <- setdiff(sex.effect.genotype.date.blocked.2,FvsV.effect.genotype.date.blocked.2)
MvsF.only.gene.list <- setdiff(MvsF.minus.FvsV.gene.list,MvsV.effect.genotype.date.blocked.2)
#' List of DE genes common to MvsF and MvsV but not present in FvsV
common.MvsF.MvsV.gene.list <- intersect(sex.effect.genotype.date.blocked.2,MvsV.effect.genotype.date.blocked.2)
common.MvsF.MvsV.minus.FvsV.gene.list <-setdiff(common.MvsF.MvsV.gene.list,FvsV.effect.genotype.date.blocked.2)
write(common.MvsF.MvsV.minus.FvsV.gene.list,
file="analysis/DESeq2/common_MvsF_MvsV_minus_FvsV_date_and_genotype-blocked_gene-list.txt")
#' List of DE genes common to MvsF and FvsV but not present in MvsV
common.MvsF.FvsV.gene.list <- intersect(sex.effect.genotype.date.blocked.2,FvsV.effect.genotype.date.blocked.2)
common.MvsF.FvsV.minus.MvsV.gene.list <-setdiff(common.MvsF.FvsV.gene.list,MvsV.effect.genotype.date.blocked.2)
write(common.MvsF.FvsV.minus.MvsV.gene.list,
      file="analysis/DESeq2/common_MvsF_FvsV_minus_MvsV_date_and_genotype-blocked_gene-list.txt")
#' List of DE genes common to FvsV and MvsV but not present in MvsF
common.FvsV.MvsV.gene.list <- intersect(FvsV.effect.genotype.date.blocked.2,MvsV.effect.genotype.date.blocked.2)
common.FvsV.MvsV.minus.MvsF.gene.list <-setdiff(common.FvsV.MvsV.gene.list,sex.effect.genotype.date.blocked.2)
write(common.FvsV.MvsV.minus.MvsF.gene.list,
      file="analysis/DESeq2/common_FvsV_MvsV_minus_MvsF_date_and_genotype-blocked_gene-list.txt")

#' To compare the dates in female samples separately
 dds <- DESeqDataSetFromMatrix(
countData = count.table,
colData = data.frame(condition=colnames(count.table),
                     sex=samples$Sex,
                     date=samples$Sampling,
                     rep=samples$Sample,
                     genotype=ifelse(samples$Sample %in% c(1,2,4),"G1","G2"),
                     Fdate=ifelse(samples$Sex %in% c("F"),samples$Sampling,"00-00")),
design = ~genotype+Fdate)
dds <- DESeq(dds)
res <- results(dds,contrast=c("Fdate","10-08","09-12"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_10-08vs09-12_female_only_genotype-blocked_s
 emi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >=1.5],file="analysis/DESeq2/differential-expression_10-08vs09-12_female_only_genotype-blocked_gene-list.txt")

#' To compare the dates in male samples separately
dds <- DESeqDataSetFromMatrix(
  countData = count.table,
  colData = data.frame(condition=colnames(count.table),
                       sex=samples$Sex,
                       date=samples$Sampling,
                       rep=samples$Sample,
                       genotype=ifelse(samples$Sample %in% c(1,2,4),"G1","G2"),
                       Mdate=ifelse(samples$Sex %in% c("M"),samples$Sampling,"00-00")),
  design = ~genotype+Mdate)
dds <- DESeq(dds)
res <- results(dds,contrast=c("Mdate","10-08","09-12"))
BiocGenerics::plotMA(res,alpha=alpha)
volcanoPlot(res,alpha=alpha)
#' Number of DE genes, padj cutoff 1e-2 + log2FoldChange cutoff 1.5
sum(res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >= 1.5)
#' Write files
write.csv2(res,file="analysis/DESeq2/differential-expression_10-08vs09-12_male_only_genotype-blocked_s
 emi-colon.csv")
write(rownames(res)[res$padj<1e-2 & !is.na(res$padj) & abs(res$log2FoldChange) >=1.5],file="analysis/DESeq2/differential-expression_10-08vs09-12_male_only_genotype-blocked_gene-list.txt")
