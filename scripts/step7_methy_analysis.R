dir.create('analysis/07_methy_analysis/')

load("raw_datas/Preprocessed/tcga_cli.RData")
smp=tcga_cli$SampleID
smp

anno <- read.delim("raw_datas/illuminaMethyl450_hg38_GDC")  
head(anno)

# samples 269
# version 07-20-2019
# type of data DNA methylation
# unit beta value
# platform Illumina Human Methylation 450

data=getTCGAMethyCpGByCode(code = 'SARC')
methy=data$M450k
dim(methy)
# [1] 395967    269


methy[c('cg00000029','cg00000165'),c('TCGA-X6-A8C5-01','TCGA-DX-A8BG-01')]
table(substr(colnames(methy),14,15))
# 01  02  06  11 
# 261   3   1   4


methy_tumor=methy[,match(smp,colnames(methy))]
methy_normal=methy[,which(substr(colnames(methy),14,15)==11)]

beta=cbind(methy_tumor,methy_normal)
group=rep(c("Tumor","Normal"),c(ncol(methy_tumor),ncol(methy_normal)))

library(minfi)
# Find differentially methylated positions
res=dmpFinder(beta,pheno = group, type = "categorical")
dim(res)
head(res)
table(res$qval<0.05)
# FALSE   TRUE 
# 391236   4731


################ limma 
beta=cbind(methy_tumor,methy_normal)
group=rep(c("Tumor","Normal"),c(ncol(methy_tumor),ncol(methy_normal)))

group <- factor(group)
design <- model.matrix(~group)
beta <- beta[complete.cases(beta), ]
mValues <- log2(beta / (1 - beta))
fit <- lmFit(mValues, design)
fit <- eBayes(fit)
results <- topTable(fit, coef="groupTumor", number=nrow(mValues))
head(results)

table(results$adj.P.Val<0.05)

###############################################
wgcna.degs.intersect.genes=read.csv('analysis/02_WGCNA/wgcna.degs.intersect.genes.csv')
head(wgcna.degs.intersect.genes)
dim(wgcna.degs.intersect.genes) ## 101
wgcna.degs.intersect.genes=as.character(wgcna.degs.intersect.genes$x)

dt=anno[anno$gene %in% wgcna.degs.intersect.genes,]
dim(dt) ## 994
head(dt)

DMP_sig=rownames(results)[which(results$adj.P.Val<0.05)]
dt=anno[anno$X.id %in% DMP_sig,]
dim(dt)
head(dt)
DMG_sig=intersect(as.character(dt$gene),wgcna.degs.intersect.genes)
# [1] "ZWINT"  "CDKN3"  "STMN1"  "CCNA2"  "TPX2"   "CCNB1"  "ARID5B" "CDK1"
length(DMG_sig)
DMG_sig
# [1] "CDKN3"  "STMN1"  "ARID5B"



library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

annotation <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

gene_beta <- aggregate(beta, by = list(annotation$UCSC_RefGene_Name), FUN = mean)

