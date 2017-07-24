args <- commandArgs(TRUE); #arg1 = count matrix  arg2=name converter arg3=cancerType
library(limma)
library(edgeR)
library(RColorBrewer)
library(gplots)
library(reshape2)

counts<-read.table(args[1], header=T, row.names=1)
colnames(counts)<-gsub("...star.star_align.", "", colnames(counts))
colnames(counts)<-gsub(".Aligned.out.bam", "", colnames(counts))
info<-counts[,1:5]
counts<-counts[,c(-1,-2,-3,-4,-5)]
dge<-DGEList(counts, genes=info)
dge<-(dge[,which(colSums(dge$counts) > 1000)])
dge<-calcNormFactors(dge)
v<-voom(dge)
#note to check this name converter stuff.  
name_converter<-read.table(args[2], header=F)
rownames(name_converter)<-name_converter[,1]
colnames(v$E)<-gsub("\\.", "-", colnames(v$E))
colnames(v$E)<-name_converter[colnames(v$E),2]
print(dim(v$E))
Gene_Expr<-melt(as.matrix(v$E))
colnames(Gene_Expr)<-c("Gene", "Sample", "logCPM")
save(Gene_Expr, file=paste(args[3], "_gencode_tmm_cpm.RData", sep=""))
#savehistory(paste(args[3], "_ucsc_tmm_cpm.Rhistory", sep=""))
