#Differential binding analysis of ChIP-seq data
#R Bioconductor package 'diffbind'

#Reference
#https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
#https://hbctraining.github.io/Intro-to-ChIPseq/lessons/08_diffbind_differential_peaks.html

#load libraries
library(DiffBind)
library(tidyverse)
library(RColorBrewer)

#read in samples
#input - BAM files (mapped only) and .narrowPeak files from MACS2
samples1 <- read.csv('sample-h3k4me3.csv')
dbObj1 <- dba(sampleSheet=samples1)

#compute count data for each consensus peak region/site
#consensus sites - identified as peaks in at least two samples
dbObj1 <- dba.count(dbObj1, bUseSummarizeOverlaps = TRUE, bParallel = FALSE)
#view FRiP (fraction of reads in peaks) values 
dbObj1
#check normalisation - default seq. depth
norm = dba.normalize(dbObj1)

#plot heatmap/PCA, based on count data
dba.plotPCA(dbObj,  attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
dev.off()
dba.plotHeatmap(dbObj, correlations = TRUE, colSideCols = c("blue", "red", "green"), col = colorRampPalette(brewer.pal(9,"Blues"))(256))

#assign pair-wise contrasts for differential binding analysis via deseq2
#e.g. WT-Mal comparison
dbmodel1 = dba.contrast(dbObj1,contrast=c("Factor","H3K4me3.Mal","H3K4me3.WT"), reorderMeta=list(Factor="H3K4me3.WT"))
#differential binding analysis
dbmodel1 <- dba.analyze(dbmodel1, method=DBA_DESEQ2)
dba.show(dbmodel1, bContrasts=T)	
#volcano plot
dba.plotVolcano(dbmodel1, contrast=1,method=DBA_DESEQ2, bUsePval = TRUE, th = 0.1, fold = 0.58496250072)

#generate reports
res_deseq1 <- dba.report(dbmodel1, method=DBA_DESEQ2, contrast = 1, th=1)

#write to file
out1 <- as.data.frame(res_deseq1)
write.table(out1, file="H3K4me3-WT-Mal.csv", sep=",", quote=F, row.names=F)



