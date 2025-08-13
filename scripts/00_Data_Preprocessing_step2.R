
library(dplyr)
library(stringr)
################# GSE21122
load("raw_datas/Preprocessed/gse21122.cli.RData")
load("raw_datas/Preprocessed/gse21122.exp.RData")

selected=c("geo_accession","disease status:ch1")
clinical_data=gse21122.cli[,selected]
colnames(clinical_data)=c("geo_accession","disease_status")
table(clinical_data$disease_status)
# Leiomyosarcoma:N/A Liposarcoma:Dedifferentiated        Liposarcoma:Myxoid/RC 
# 26                           46                           20 
# Liposarcoma:Pleomorphic         MFH:Myxofibrosarcoma              MFH:Pleomorphic 
# 23                           31                            3 
# Normal control 
# 9 
unique(clinical_data$disease_status)
# "Liposarcoma:Dedifferentiated": 
# "Leiomyosarcoma:N/A": 
# "MFH:Myxofibrosarcoma": 
# "MFH:Pleomorphic": 
# "Liposarcoma:Myxoid/RC": 
# "Liposarcoma:Pleomorphic": 
# "Normal control": 

clinical_data=clinical_data %>% mutate(disease_status=case_when(disease_status=="Liposarcoma:Dedifferentiated" ~ "Liposarcoma:Dedifferentiated liposarcoma",
                                          disease_status=="Leiomyosarcoma" ~ "Leiomyosarcoma",
                                          disease_status=="MFH:Myxofibrosarcoma" ~ "Myxofibrosarcoma:Myxofibrosarcoma",
                                          disease_status=="MFH:Pleomorphic" ~ "Myxofibrosarcoma:Pleomorphic MFH",
                                          disease_status=="Liposarcoma:Myxoid/RC" ~ "Liposarcoma:Myxoid/round cell liposarcoma",
                                          disease_status=="Liposarcoma:Pleomorphic" ~ "Liposarcoma:Pleomorphic liposarcoma",
                                          disease_status=="Normal control" ~ "Normal control"))
head(clinical_data)
clinical_data= clinical_data %>% tidyr::separate(disease_status, c("sample_type","histology"),sep = ":")
table(clinical_data$sample_type,clinical_data$histology)


gse21122_cli=clinical_data %>% filter(sample_type %in% c("Liposarcoma","Normal control"))
gse21122_exprs=gse21122.exp %>% dplyr::select(gse21122_cli$geo_accession)
dim(gse21122_cli)
dim(gse21122_exprs)

################# GSE159659
load("raw_datas/Preprocessed/gse159659.cli.RData")
load("raw_datas/Preprocessed/gse159659.exp.RData")

selected=c("geo_accession","tissue:ch1")
clinical_data=gse159659.cli[,selected]
colnames(clinical_data)=c("geo_accession","disease_status")
table(clinical_data$disease_status)
# adipose tissue    dedifferentiated liposarcoma well differentiated liposarcoma 
# 15                              15                              15 
unique(clinical_data$disease_status)
# "dedifferentiated liposarcoma"
# "well differentiated liposarcoma"
# "adipose tissue":

clinical_data=clinical_data %>% mutate(disease_status=case_when(disease_status=="dedifferentiated liposarcoma" ~ "Liposarcoma:Dedifferentiated liposarcoma",
                                                                disease_status=="well differentiated liposarcoma" ~ "Liposarcoma:Well-differentiated liposarcoma",
                                                                disease_status=="adipose tissue" ~ "Normal control"))
head(clinical_data)
clinical_data= clinical_data %>% tidyr::separate(disease_status, c("sample_type","histology"),sep = ":")
table(clinical_data$sample_type)
table(clinical_data$sample_type,clinical_data$histology)


gse159659_cli=clinical_data %>% filter(sample_type %in% c("Liposarcoma","Normal control"))
gse159659_exprs=gse159659.exp %>% dplyr::select(gse159659_cli$geo_accession)
dim(gse159659_cli)
dim(gse159659_exprs)


################# GSE30929
load("raw_datas/Preprocessed/gse30929.cli.RData")
load("raw_datas/Preprocessed/gse30929.exp.RData")

selected=c("geo_accession","tissue:ch1","subtype:ch1","drfs:ch1","tt.drfs:ch1")
clinical_data=gse30929.cli[,selected]
# Distant recurrence-free survival (DRFS)
colnames(clinical_data)=c("geo_accession","tissue","histology","DRFS","DRFS.time")
table(clinical_data$tissue)
# primary liposarcoma 
# 140
table(clinical_data$histology)
# dedifferentiated              myxoid   myxoid/round cell         pleomorphic well-differentiated 
# 40                  17                  11                  20                  52
unique(clinical_data$histology)

clinical_data=clinical_data %>% mutate(histology=case_when(histology=="dedifferentiated" ~ "Liposarcoma:Dedifferentiated liposarcoma",
                                                           histology=="myxoid" ~ "Liposarcoma:Myxoid",
                                                           histology=="myxoid/round cell" ~ "Liposarcoma:Myxoid/round cell liposarcoma",
                                                           histology=="pleomorphic" ~ "Liposarcoma:Pleomorphic liposarcoma",
                                                           histology=="well-differentiated" ~ "Liposarcoma:Well-differentiated liposarcoma"))
clinical_data$DRFS[which(clinical_data$DRFS=="FALSE")]=0
clinical_data$DRFS[which(clinical_data$DRFS=="TRUE")]=1
head(clinical_data)
clinical_data= clinical_data %>% tidyr::separate(histology, c("sample_type","histology"),sep = ":")
table(clinical_data$sample_type,clinical_data$histology)


gse30929_cli=clinical_data %>% filter(sample_type %in% c("Liposarcoma","Normal control"))
gse30929_exprs=gse30929.exp %>% dplyr::select(gse30929_cli$geo_accession)
dim(gse30929_cli)
dim(gse30929_exprs)
# [1] 12548   140

gse30929_cli$DRFS.time=gse30929_cli$DRFS.time/12

################# GSE21050
load("raw_datas/Preprocessed/gse21050.cli.RData")
load("raw_datas/Preprocessed/gse21050.exp.RData")

selected=c("geo_accession","tissue:ch1","diagnosis:ch1","metastasis:ch1","time:ch1")
clinical_data=gse21050.cli[,selected]
# Distant recurrence-free survival (DRFS): 
colnames(clinical_data)=c("geo_accession","tissue","histology","DRFS","DRFS.time")
table(clinical_data$tissue)
# Extremities  Head and neck Internal trunk     Trunk wall 
# 187              3             66             54 
# "Internal trunk": 
# "Extremities": 
# "Trunk wall": 
# "Head and neck":
table(clinical_data$histology)
# Leiomyosarcoma Liposarcoma - dedifferentiated                          Other 
# 85                             62                             27 
# Undifferentiated sarcoma 
# 136 
unique(clinical_data$histology)
# "Liposarcoma - dedifferentiated"
# "Leiomyosarcoma": 
# "Undifferentiated sarcoma"
# "Other"

clinical_data=clinical_data %>% mutate(histology=case_when(histology=="Liposarcoma - dedifferentiated" ~ "Liposarcoma:Dedifferentiated liposarcoma",
                                                           histology=="Leiomyosarcoma" ~ "Leiomyosarcoma",
                                                           histology=="Undifferentiated sarcoma" ~ "Undifferentiated sarcoma",
                                                           histology=="Other" ~ "Other"))
clinical_data$DRFS[which(clinical_data$DRFS=="no")]=0
clinical_data$DRFS[which(clinical_data$DRFS=="yes")]=1
head(clinical_data)
clinical_data= clinical_data %>% tidyr::separate(histology, c("sample_type","histology"),sep = ":")
table(clinical_data$sample_type,clinical_data$histology)


gse21050_cli=clinical_data %>% filter(sample_type %in% c("Liposarcoma","Normal control"))
gse21050_exprs=gse21050.exp %>% dplyr::select(gse21050_cli$geo_accession)
dim(gse21050_cli)
dim(gse21050_exprs)
# [1] 12548   62

################# GSE71118
load("raw_datas/Preprocessed/gse71118.cli.RData")
load("raw_datas/Preprocessed/gse71118.exp.RData")

selected=c("geo_accession","source_name_ch1","cinsarc:ch1","metastasis:ch1","time:ch1")
clinical_data=gse71118.cli[,selected]
# Distant recurrence-free survival (DRFS): 
colnames(clinical_data)=c("geo_accession","histology","cinsarc","DRFS","DRFS.time")
table(clinical_data$histology)
# Dedifferentiated liposarcoma               Leiomyosarcoma                  Liposarcoma 
# 44                           89                            1 
# Myxofibrosarcoma                        Other     Undifferentiated sarcoma 
# 43                           46                           89 
unique(clinical_data$histology)
# "Dedifferentiated liposarcoma": 
# "Myxofibrosarcoma":
# "Undifferentiated sarcoma": 
# "Leiomyosarcoma": 
# "Other"
# "Liposarcoma": 


clinical_data=clinical_data %>% dplyr::mutate(histology=case_when(histology=="Dedifferentiated liposarcoma" ~ "Liposarcoma:Dedifferentiated liposarcoma",
                                                           histology=="Myxofibrosarcoma" ~ "Myxofibrosarcoma",
                                                           histology=="Undifferentiated sarcoma" ~ "Undifferentiated sarcoma",
                                                           histology=="Leiomyosarcoma" ~ "Leiomyosarcoma",
                                                           histology=="Other" ~ "Other",
                                                           # histology=="Liposarcoma" ~ "Liposarcoma:Liposarcoma"
                                                           histology=="Liposarcoma" ~ "Other"))
clinical_data$DRFS[which(clinical_data$DRFS=="No")]=0
clinical_data$DRFS[which(clinical_data$DRFS=="Yes")]=1
head(clinical_data)
clinical_data= clinical_data %>% tidyr::separate(histology, c("sample_type","histology"),sep = ":")
table(clinical_data$sample_type,clinical_data$histology)


gse71118_cli=clinical_data %>% filter(sample_type %in% c("Liposarcoma","Normal control"))
gse71118_exprs=gse71118.exp %>% dplyr::select(gse71118_cli$geo_accession)
dim(gse71118_cli)
dim(gse71118_exprs)
# [1] 21655   44

save(gse21122_cli,file = "raw_datas/Preprocessed/gse21122_cli.RData")
save(gse21122_exprs,file = "raw_datas/Preprocessed/gse21122_exprs.RData")
save(gse159659_cli,file = "raw_datas/Preprocessed/gse159659_cli.RData")
save(gse159659_exprs,file = "raw_datas/Preprocessed/gse159659_exprs.RData")
save(gse30929_cli,file = "raw_datas/Preprocessed/gse30929_cli.RData")
save(gse30929_exprs,file = "raw_datas/Preprocessed/gse30929_exprs.RData")
save(gse21050_cli,file = "raw_datas/Preprocessed/gse21050_cli.RData")
save(gse21050_exprs,file = "raw_datas/Preprocessed/gse21050_exprs.RData")
save(gse71118_cli,file = "raw_datas/Preprocessed/gse71118_cli.RData")
save(gse71118_exprs,file = "raw_datas/Preprocessed/gse71118_exprs.RData")

##################################### TCGA
load("raw_datas/Preprocessed/tcga.expr.tpm.RData")
load("raw_datas/Preprocessed/tcga.cli.RData")

tcga.exprs=log2(expr+1)
tcga.cli=tcga.cli
dim(tcga.exprs)
colnames(tcga.exprs)=substr(colnames(tcga.exprs),1,15)
table(substr(colnames(tcga.exprs),14,15))
# 01  02  06  11 
# 259   3   1   2

clinical_data= tcga.cli %>% tidyr::separate(histology, c("sample_type","histology"),sep = ":")
table(clinical_data$sample_type,clinical_data$histology)


tcga_cli=clinical_data %>% filter(sample_type %in% c("Liposarcoma","Normal control"))
tcga_cli$SampleID=paste0(tcga_cli$SampleID,"-01")
rownames(tcga_cli)=tcga_cli$SampleID
dim(tcga_cli)

smp=intersect(tcga_cli$SampleID,colnames(tcga.exprs))
length(smp)
# [1] 58

tcga_cli=tcga_cli[match(smp,tcga_cli$SampleID),]
tcga_exprs=tcga.exprs %>% dplyr::select(smp)
dim(tcga_exprs)
# [1] 19938    58

all(colnames(tcga_exprs)==tcga_cli$SampleID)

save(tcga_cli,file = "raw_datas/Preprocessed/tcga_cli.RData")
save(tcga_exprs,file = "raw_datas/Preprocessed/tcga_exprs.RData")


# GSE21122
# GSE159659
# GSE30929
# GSE21050

##################### combat
library("FactoMineR")
library("factoextra")
pca.plot = function(dat,col,title ="PCA - Biplot"){
  df.pca <- PCA(t(dat), graph = FALSE)
  fviz_pca_ind(df.pca,
               geom.ind = "point",
               col.ind = col ,
               addEllipses = TRUE,
               legend.title = "Groups",
               title=title)
}

dim(gse21122_exprs)
dim(gse159659_exprs)
dim(gse30929_exprs)
dim(gse21050_exprs)

common.genes=Reduce(intersect,list(rownames(gse21122_exprs),
                                   rownames(gse159659_exprs)))
length(common.genes)
# [1] 11771

data=cbind(gse21122_exprs[common.genes,],
           gse159659_exprs[common.genes,])
data=as.matrix(data)
dim(data)

########
group=rep(c('GSE21122','GSE159659'),
          times=c(ncol(gse21122_exprs),ncol(gse159659_exprs)))
group

sampleInfor=data.frame(ID=colnames(data),batch=group,stringsAsFactors = F)
rownames(sampleInfor)=sampleInfor$ID
head(sampleInfor)

clinical_data=rbind(gse21122_cli,gse159659_cli)
clinical_data$batch=sampleInfor$batch
head(clinical_data)

exprs=data
exprs=as.matrix(exprs)
batch=as.character(sampleInfor[colnames(exprs),'batch'])

pdf("raw_datas/Preprocessed/ad.pca.before.pdf",height = 4,width = 4)
pca.plot(dat = exprs,col = batch)
dev.off()


library(sva)
expr_adjusted=ComBat(dat = exprs
                     , batch = batch
                     , par.prior = T)

pdf("raw_datas/Preprocessed/ad.pca.pdf",height = 4,width = 4)
pca.plot(dat = expr_adjusted,col = batch)
dev.off()

save(expr_adjusted,file = 'raw_datas/Preprocessed/gse21122.gse159659.expr_adjusted.RData')
save(clinical_data,file = 'raw_datas/Preprocessed/gse21122.gse159659.clinical_data.RData')

################################## GSE30929 和 GSE21050
load("raw_datas/Preprocessed/gse71118_cli.RData")
load("raw_datas/Preprocessed/gse71118_exprs.RData")
load("raw_datas/Preprocessed/gse21050_cli.RData")
load("raw_datas/Preprocessed/gse21050_exprs.RData")

dim(gse71118_exprs)
dim(gse21050_exprs)

common.genes=Reduce(intersect,list(rownames(gse71118_exprs),
                                   rownames(gse21050_exprs)))
length(common.genes)
# [1] 21655

data=cbind(gse71118_exprs[common.genes,],
           gse21050_exprs[common.genes,])
data=as.matrix(data)
dim(data)

########
group=rep(c('GSE71118','GSE21050'),
          times=c(ncol(gse71118_exprs),ncol(gse21050_exprs)))
group

sampleInfor=data.frame(ID=colnames(data),batch=group,stringsAsFactors = F)
rownames(sampleInfor)=sampleInfor$ID
head(sampleInfor)

clinical_data=bind_rows(gse71118_cli,gse21050_cli)
clinical_data$batch=sampleInfor$batch
head(clinical_data)

exprs=data
exprs=as.matrix(exprs)
batch=as.character(sampleInfor[colnames(exprs),'batch'])

pdf("raw_datas/Preprocessed/ad.pca.before_2.pdf",height = 4,width = 4)
pca.plot(dat = exprs,col = batch)
dev.off()


library(sva)
expr_adjusted=ComBat(dat = exprs
                     , batch = batch
                     , par.prior = T)

pdf("raw_datas/Preprocessed/ad.pca_2.pdf",height = 4,width = 4)
pca.plot(dat = expr_adjusted,col = batch)
dev.off()

save(expr_adjusted,file = 'raw_datas/Preprocessed/gse71118.gse21050.expr_adjusted.RData')
save(clinical_data,file = 'raw_datas/Preprocessed/gse71118.gse21050.clinical_data.RData')

load("raw_datas/Preprocessed/gse71118_cli.RData")
load("raw_datas/Preprocessed/gse71118_exprs.RData")
load("raw_datas/Preprocessed/gse21050_cli.RData")
load("raw_datas/Preprocessed/gse21050_exprs.RData")
load("raw_datas/Preprocessed/gse30929_cli.RData")
load("raw_datas/Preprocessed/gse30929_exprs.RData")

dim(gse71118_exprs)
dim(gse21050_exprs)

common.genes=Reduce(intersect,list(rownames(gse71118_exprs),
                                   rownames(gse21050_exprs),
                                   rownames(gse30929_exprs)))
length(common.genes)
# [1] 12548

data=cbind(gse71118_exprs[common.genes,],
           gse21050_exprs[common.genes,],
           gse30929_exprs[common.genes,])
data=as.matrix(data)
dim(data)

########
group=rep(c('GSE71118','GSE21050','GSE30929'),
          times=c(ncol(gse71118_exprs),ncol(gse21050_exprs),ncol(gse30929_exprs)))
group

sampleInfor=data.frame(ID=colnames(data),batch=group,stringsAsFactors = F)
rownames(sampleInfor)=sampleInfor$ID
head(sampleInfor)

clinical_data=bind_rows(gse71118_cli,gse21050_cli,gse30929_cli)
clinical_data$batch=sampleInfor$batch
head(clinical_data)

exprs=data
exprs=as.matrix(exprs)
batch=as.character(sampleInfor[colnames(exprs),'batch'])

pdf("raw_datas/Preprocessed/ad.pca.before_3.pdf",height = 4,width = 4)
pca.plot(dat = exprs,col = batch)
dev.off()


library(sva)
expr_adjusted=ComBat(dat = exprs
                     , batch = batch
                     , par.prior = T)

pdf("raw_datas/Preprocessed/ad.pca_3.pdf",height = 4,width = 4)
pca.plot(dat = expr_adjusted,col = batch)
dev.off()

save(expr_adjusted,file = 'raw_datas/Preprocessed/gse71118.gse21050.gse30929.expr_adjusted.RData')
save(clinical_data,file = 'raw_datas/Preprocessed/gse71118.gse21050.gse30929.clinical_data.RData')

dt=data.frame(exprs=(t(exprs_data)[,seleceted_genes]),group=clinical_sub$sample_type)
head(dt)

library(ggplot2)
library(gghalves)
library(ggpubr)
mycol = c("#ff7b00","#006fbf")

p=ggplot(data = dt, aes(x = group, y = exprs, fill = group)) +
  geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
  geom_point(aes(y = exprs, color = group), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
  labs(y = paste("the expression of TXNDC5","\n"), x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = mycol) +
  scale_colour_manual(values = mycol) +
  theme_classic()
p=p+ggpubr::stat_compare_means(comparisons = list(c(1,2)),method = 'wilcox.test',label= "p.format")
