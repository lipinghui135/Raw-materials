dir.create('analysis/06_genomic_landscape/')
########################################################
load("raw_datas/Preprocessed/tcga_cli.RData")
###############################################
wgcna.degs.intersect.genes=read.csv('analysis/02_WGCNA/wgcna.degs.intersect.genes.csv')
head(wgcna.degs.intersect.genes)
dim(wgcna.degs.intersect.genes) ## 101
wgcna.degs.intersect.genes=as.character(wgcna.degs.intersect.genes$x)

library(maftools)
load(file = "GDCdata/TCGA-SARC_SNP.Rdata")
maf.sarc <- data
maf.sarc$Tumor_Sample_Barcode=substr(maf.sarc$Tumor_Sample_Barcode,1,15)
maf <- read.maf(maf.sarc)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

maf_sub=subsetMaf(maf,tsb=tcga_cli$SampleID)
pdf("analysis/06_genomic_landscape/oncoplot.pdf",height = 6,width = 6)
oncoplot(maf_sub,genes=wgcna.degs.intersect.genes)
dev.off()


