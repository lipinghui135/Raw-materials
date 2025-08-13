setwd("Z:/users/lishuang/Work1/Work/20240524_Liposarcomas")
source('Z:/projects/codes/mg_base.R')
##################
library(IOBR)
load("raw_datas/Preprocessed/gse71118.gse21050.gse30929.expr_adjusted.RData")
expr=expr_adjusted
dim(expr) ## 21655   107
range(expr)
#### ESTIMATE
geo.exp.estimate<-deconvo_estimate(eset=expr)
save(geo.exp.estimate,file='raw_datas/TME/geo.exp.estimate.RData')
### CIBERSORT
geo.exp.cibersort<-deconvo_cibersort(eset=expr,arrays=T)
save(geo.exp.cibersort,file='raw_datas/TME/geo.exp.cibersort.RData')
############ TIMER 
geo.exp.timer<-deconvo_timer(eset=as.matrix(expr),indications=rep('sarc',ncol(expr)))
save(geo.exp.timer,file='raw_datas/TME/geo.exp.timer.RData')
############ EPIC 
geo.exp.epic<-deconvo_epic(eset=as.matrix(expr),tumor = TRUE)
save(geo.exp.epic,file='raw_datas/TME/geo.exp.epic.RData')
############ MCP-counter 
geo.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(expr))
save(geo.exp.mcp,file='raw_datas/TME/geo.exp.mcp.RData')
############ Xcell
geo.exp.xcell=deconvo_xcell(eset=expr,project = 'Liposarcoma',arrays = T)
save(geo.exp.xcell,file='raw_datas/TME/geo.exp.xcell.RData')

################################ ssGSEA
geneset=readxl::read_excel('Z:/users/lishuang/public/immu_signature/29_immu_signature_pmid_30594216.xlsx',col_names = T,skip =1)
geneset=data.frame(geneset,check.names = F)
dim(geneset)
apply(geneset, 2, function(x){length(na.omit(x))})

geneList=apply(geneset, 2, function(x){as.character(na.omit(x))})

geo.immu.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = expr,genelist = geneList)
dim(geo.immu.ssgsea)
save(geo.immu.ssgsea,file = 'raw_datas/TME/geo.immu.ssgsea.RData')
#####################
geneset=read.table('Z:/users/lishuang/public/immu_signature/28_immu_signature.PMID28052254.txt',sep = "\t",header = T,stringsAsFactors = F)
geneList=split(geneset$Gene,geneset$SetName)

geo.immu.ssgsea_28=ssGSEAScore_by_muti_group_genes(gene.exp = expr,genelist = geneList)
dim(geo.immu.ssgsea_28)
save(geo.immu.ssgsea_28,file = 'raw_datas/TME/geo.immu.ssgsea_28.RData')

######################## GSE30929
library(IOBR)
load("raw_datas/Preprocessed/gse30929_exprs.RData")
expr=gse30929_exprs
dim(expr) ## 21655   107
range(expr)
#### ESTIMATE
gse30929.exp.estimate<-deconvo_estimate(eset=expr)
save(gse30929.exp.estimate,file='raw_datas/TME/gse30929.exp.estimate.RData')
### CIBERSORT
gse30929.exp.cibersort<-deconvo_cibersort(eset=expr,arrays=T)
save(gse30929.exp.cibersort,file='raw_datas/TME/gse30929.exp.cibersort.RData')
############ TIMER 
gse30929.exp.timer<-deconvo_timer(eset=as.matrix(expr),indications=rep('sarc',ncol(expr)))
save(gse30929.exp.timer,file='raw_datas/TME/gse30929.exp.timer.RData')
############ EPIC 
gse30929.exp.epic<-deconvo_epic(eset=as.matrix(expr),tumor = TRUE)
save(gse30929.exp.epic,file='raw_datas/TME/gse30929.exp.epic.RData')
############ MCP-counter 
gse30929.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(expr))
save(gse30929.exp.mcp,file='raw_datas/TME/gse30929.exp.mcp.RData')
############ Xcell
gse30929.exp.xcell=deconvo_xcell(eset=expr,project = 'Liposarcoma',arrays = T)
save(gse30929.exp.xcell,file='raw_datas/TME/gse30929.exp.xcell.RData')

################################ ssGSEA
geneset=readxl::read_excel('Z:/users/lishuang/public/immu_signature/29_immu_signature_pmid_30594216.xlsx',col_names = T,skip =1)
geneset=data.frame(geneset,check.names = F)
dim(geneset)
apply(geneset, 2, function(x){length(na.omit(x))})

geneList=apply(geneset, 2, function(x){as.character(na.omit(x))})

gse30929.immu.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = expr,genelist = geneList)
dim(gse30929.immu.ssgsea)
save(gse30929.immu.ssgsea,file = 'raw_datas/TME/gse30929.immu.ssgsea.RData')
#####################
geneset=read.table('Z:/users/lishuang/public/immu_signature/28_immu_signature.PMID28052254.txt',sep = "\t",header = T,stringsAsFactors = F)
geneList=split(geneset$Gene,geneset$SetName)

gse30929.immu.ssgsea_28=ssGSEAScore_by_muti_group_genes(gene.exp = expr,genelist = geneList)
dim(gse30929.immu.ssgsea_28)
save(gse30929.immu.ssgsea_28,file = 'raw_datas/TME/gse30929.immu.ssgsea_28.RData')

######################## TCGA
library(IOBR)
load("raw_datas/Preprocessed/tcga_exprs.RData")
expr=tcga_exprs
dim(expr) ## 21655   107
range(expr)
#### ESTIMATE
tcga.exp.estimate<-deconvo_estimate(eset=expr)
save(tcga.exp.estimate,file='raw_datas/TME/tcga.exp.estimate.RData')
### CIBERSORT
tcga.exp.cibersort<-deconvo_cibersort(eset=expr,arrays=F)
save(tcga.exp.cibersort,file='raw_datas/TME/tcga.exp.cibersort.RData')
############ TIMER 
tcga.exp.timer<-deconvo_timer(eset=as.matrix(expr),indications=rep('sarc',ncol(expr)))
save(tcga.exp.timer,file='raw_datas/TME/tcga.exp.timer.RData')
############ EPIC 
tcga.exp.epic<-deconvo_epic(eset=as.matrix(expr),tumor = TRUE)
save(tcga.exp.epic,file='raw_datas/TME/tcga.exp.epic.RData')
############ MCP-counter 
tcga.exp.mcp<-deconvo_mcpcounter(eset=as.matrix(expr))
save(tcga.exp.mcp,file='raw_datas/TME/tcga.exp.mcp.RData')
############ Xcell
tcga.exp.xcell=deconvo_xcell(eset=expr,project = 'Liposarcoma',arrays = F)
save(tcga.exp.xcell,file='raw_datas/TME/tcga.exp.xcell.RData')

################################ ssGSEA
geneset=readxl::read_excel('Z:/users/lishuang/public/immu_signature/29_immu_signature_pmid_30594216.xlsx',col_names = T,skip =1)
geneset=data.frame(geneset,check.names = F)
dim(geneset)
apply(geneset, 2, function(x){length(na.omit(x))})

geneList=apply(geneset, 2, function(x){as.character(na.omit(x))})

tcga.immu.ssgsea=ssGSEAScore_by_muti_group_genes(gene.exp = expr,genelist = geneList)
dim(tcga.immu.ssgsea)
save(tcga.immu.ssgsea,file = 'raw_datas/TME/tcga.immu.ssgsea.RData')
#####################
geneset=read.table('Z:/users/lishuang/public/immu_signature/28_immu_signature.PMID28052254.txt',sep = "\t",header = T,stringsAsFactors = F)
geneList=split(geneset$Gene,geneset$SetName)

tcga.immu.ssgsea_28=ssGSEAScore_by_muti_group_genes(gene.exp = expr,genelist = geneList)
dim(tcga.immu.ssgsea_28)
save(tcga.immu.ssgsea_28,file = 'raw_datas/TME/tcga.immu.ssgsea_28.RData')

