
library(GEOquery)
library(dplyr)
################################## GSE21122
gset=getGEO('GSE21122')
gse21122.cli=pData(gset$GSE21122_series_matrix.txt.gz)
table(gse21122.cli$`disease status:ch1`)
# Leiomyosarcoma:N/A Liposarcoma:Dedifferentiated        Liposarcoma:Myxoid/RC 
# 26                           46                           20 
# Liposarcoma:Pleomorphic         MFH:Myxofibrosarcoma              MFH:Pleomorphic 
# 23                           31                            3 
# Normal control 
# 9


library(stringr)
gpl=fData(gset$GSE21122_series_matrix.txt.gz)
dt=data.frame(ID=gpl$ID,Symbol=gpl$`Gene Symbol`)
colnames(dt)[2]='Symbol'
head(dt)

probe<-dt %>% 
  select("ID","Symbol") %>%
  filter(Symbol!="") %>%
  filter(!grepl("///",Symbol))
dim(probe)
# [1] 19820     2
head(probe)


library(affy)
####
celpath="raw_datas/GEO/GSE21122_RAW/"
## read all files in that dir
celfiles <- ReadAffy(celfile.path=celpath)
# boxplot(celfiles)
# hist(celfiles)
eset <- rma(celfiles)
# Background correcting
# Normalizing
# Calculating Expression
exprData=exprs(eset)
colnames(exprData)=gsub("(GSM[0-9]+).*","\\1",colnames(exprData))
dim(exprData)
exprData[1:4,1:5]
range(exprData)
# [1]  3.675598 14.956062
class(exprData)
exprData=data.frame(exprData)
exprData$ID = rownames(exprData)

exprData = merge(exprData, probe, by = "ID")
head(exprData)
dim(exprData)
# [1] 19820   160

exprData = exprData %>% 
  dplyr::select(-ID) %>% 
  dplyr::select(Symbol, everything())
dim(exprData)
# [1] 19820   160

exprData[1:3,1:5]
table(duplicated(exprData$Symbol))

# FALSE  TRUE 
# 12548  7272

exprData <- exprData %>% 
  group_by(Symbol) %>%
  summarise_all(mean)
dim(exprData)
# [1] 12548   159

exprData=data.frame(exprData,check.names=F)
exprData[1:3,1:5]
rownames(exprData)=exprData$Symbol
exprData=exprData[,-1]
exprData[1:3,1:5]

gse21122.exp=exprData
dim(gse21122.exp)
# [1] 12548   158

save(gse21122.exp,file = 'raw_datas/Preprocessed/gse21122.exp.RData')
save(gse21122.cli,file = 'raw_datas/Preprocessed/gse21122.cli.RData')

################################## GSE159659
gset=getGEO('GSE159659')
gse159659.cli=pData(gset$GSE159659_series_matrix.txt.gz)
table(gse159659.cli$`tissue:ch1`)
# adipose tissue    dedifferentiated liposarcoma well differentiated liposarcoma 
# 15                              15                              15 


library(stringr)
gpl=fData(gset$GSE159659_series_matrix.txt.gz)
Symbol=str_split(gpl$SPOT_ID.1,pattern = " // ",simplify = T)[,3]
pattern <- ".*\\((?<ID>[A-Za-z0-9]*)\\),.*"
res <- stringr::str_match(string = Symbol, pattern = pattern)
Symbol <- res[,2]
dt=data.frame(ID=gpl$ID,Symbol)
head(dt)

probe<-dt %>% 
  select("ID","Symbol") %>%
  filter(!is.na(Symbol)) %>%
  filter(Symbol!="") %>%
  filter(!grepl("///",Symbol))
dim(probe)
# [1] 18838     2
head(probe)


library(affy)
####
celpath="raw_datas/GEO/GSE159659_RAW/"
## read all files in that dir
files=list.files(celpath,full.names = T)
files
celfiles <- oligo::read.celfiles(filenames = files)
# boxplot(celfiles)
# hist(celfiles)
eset <- oligo::rma(celfiles)
# Background correcting
# Normalizing
# Calculating Expression
exprData=oligo::exprs(eset)
colnames(exprData)=gsub("(GSM[0-9]+).*","\\1",colnames(exprData))
dim(exprData)
exprData[1:4,1:5]
range(exprData)
# [1]  1.985236 13.994767

class(exprData)
exprData=data.frame(exprData)
exprData$ID = rownames(exprData)

exprData = merge(exprData, probe, by = "ID")
head(exprData)
dim(exprData)
# [1] 18838    47

exprData = exprData %>% 
  dplyr::select(-ID) %>% 
  dplyr::select(Symbol, everything())
dim(exprData)
# [1] 18838    46

exprData[1:3,1:5]
table(duplicated(exprData$Symbol))

# FALSE 
# 18838

exprData <- exprData %>% 
  group_by(Symbol) %>%
  summarise_all(mean)
dim(exprData)
# [1] 18838    46

exprData=data.frame(exprData,check.names=F)
exprData[1:3,1:5]
rownames(exprData)=exprData$Symbol
exprData=exprData[,-1]
exprData[1:3,1:5]

gse159659.exp=exprData
dim(gse159659.exp)
# [1] 19360   45

save(gse159659.exp,file = 'raw_datas/Preprocessed/gse159659.exp.RData')
save(gse159659.cli,file = 'raw_datas/Preprocessed/gse159659.cli.RData')

################################## GSE30929
gset=getGEO('GSE30929')
gse30929.cli=pData(gset$GSE30929_series_matrix.txt.gz)
table(gse30929.cli$`tissue:ch1`)
# primary liposarcoma 
# 140 


library(stringr)
gpl=fData(gset$GSE30929_series_matrix.txt.gz)
dt=data.frame(ID=gpl$ID,Symbol=gpl$`Gene Symbol`)
colnames(dt)[2]='Symbol'
head(dt)

probe<-dt %>% 
  select("ID","Symbol") %>%
  filter(Symbol!="") %>%
  filter(!grepl("///",Symbol))
dim(probe)
# [1] 19820     2
head(probe)


library(affy)
####
celpath="raw_datas/GEO/GSE30929_RAW/"
## read all files in that dir
celfiles <- ReadAffy(celfile.path=celpath)
# boxplot(celfiles)
# hist(celfiles)
eset <- affy::rma(celfiles)
# Background correcting
# Normalizing
# Calculating Expression
exprData=exprs(eset)
colnames(exprData)=gsub("(GSM[0-9]+).*","\\1",colnames(exprData))
dim(exprData)
exprData[1:4,1:5]
range(exprData)
# [1]  2.243086 13.302879
class(exprData)
exprData=data.frame(exprData)
exprData$ID = rownames(exprData)

exprData = merge(exprData, probe, by = "ID")
head(exprData)
dim(exprData)
# [1] 19820   142

exprData = exprData %>% 
  dplyr::select(-ID) %>% 
  dplyr::select(Symbol, everything())
dim(exprData)
# [1] 19820   141

exprData[1:3,1:5]
table(duplicated(exprData$Symbol))

# FALSE  TRUE 
# 12548  7272

exprData <- exprData %>% 
  group_by(Symbol) %>%
  summarise_all(mean)
dim(exprData)
# [1] 12548   141

exprData=data.frame(exprData,check.names=F)
exprData[1:3,1:5]
rownames(exprData)=exprData$Symbol
exprData=exprData[,-1]
exprData[1:3,1:5]

gse30929.exp=exprData
dim(gse30929.exp)
# [1] 12548   140

save(gse30929.exp,file = 'raw_datas/Preprocessed/gse30929.exp.RData')
save(gse30929.cli,file = 'raw_datas/Preprocessed/gse30929.cli.RData')

################################## GSE21050
gset=getGEO('GSE21050')
gse21050.cli=pData(gset$GSE21050_series_matrix.txt.gz)
table(gse21050.cli$`diagnosis:ch1`)
# Leiomyosarcoma Liposarcoma - dedifferentiated                          Other 
# 85                             62                             27 
# Undifferentiated sarcoma 
# 136 


library(stringr)
gpl=fData(gset$GSE21050_series_matrix.txt.gz)
dt=data.frame(ID=gpl$ID,Symbol=gpl$`Gene Symbol`)
colnames(dt)[2]='Symbol'
head(dt)

probe<-dt %>% 
  select("ID","Symbol") %>%
  filter(Symbol!="") %>%
  filter(!grepl("///",Symbol))
dim(probe)
# [1] 42976     2
head(probe)


library(affy)
####
celpath="raw_datas/GEO/GSE21050_RAW/"
## read all files in that dir
celfiles <- ReadAffy(celfile.path=celpath)
# boxplot(celfiles)
# hist(celfiles)
eset <- affy::rma(celfiles)
# Background correcting
# Normalizing
# Calculating Expression
exprData=exprs(eset)
colnames(exprData)=gsub("(GSM[0-9]+).*","\\1",colnames(exprData))
dim(exprData)
exprData[1:4,1:5]
range(exprData)
# [1]  2.387562 14.912034

class(exprData)
exprData=data.frame(exprData)
exprData$ID = rownames(exprData)

exprData = merge(exprData, probe, by = "ID")
head(exprData)
dim(exprData)
# [1] 42976   312

exprData = exprData %>% 
  dplyr::select(-ID) %>% 
  dplyr::select(Symbol, everything())
dim(exprData)
# [1] 42976   311

exprData[1:3,1:5]
table(duplicated(exprData$Symbol))

# FALSE  TRUE 
# 21655 21321

exprData <- exprData %>% 
  group_by(Symbol) %>%
  summarise_all(mean)
dim(exprData)
# [1] 21655   311

exprData=data.frame(exprData,check.names=F)
exprData[1:3,1:5]
rownames(exprData)=exprData$Symbol
exprData=exprData[,-1]
exprData[1:3,1:5]

gse21050.exp=exprData
dim(gse21050.exp)
# [1] 21655   311

save(gse21050.exp,file = 'raw_datas/Preprocessed/gse21050.exp.RData')
save(gse21050.cli,file = 'raw_datas/Preprocessed/gse21050.cli.RData')

################################## GSE71118
gset=getGEO('GSE71118')
gse71118.cli=pData(gset$GSE71118_series_matrix.txt.gz)
table(gse71118.cli$source_name_ch1)
# Dedifferentiated liposarcoma               Leiomyosarcoma                  Liposarcoma 
# 44                           89                            1 
# Myxofibrosarcoma                        Other     Undifferentiated sarcoma 
# 43                           46                           89


library(stringr)
gpl=fData(gset$GSE71118_series_matrix.txt.gz)
dt=data.frame(ID=gpl$ID,Symbol=gpl$`Gene Symbol`)
colnames(dt)[2]='Symbol'
head(dt)

probe<-dt %>% 
  select("ID","Symbol") %>%
  filter(Symbol!="") %>%
  filter(!grepl("///",Symbol))
dim(probe)
# [1] 42976     2
head(probe)


library(affy)
####
celpath="raw_datas/GEO/GSE71118_RAW/"
## read all files in that dir
celfiles <- ReadAffy(celfile.path=celpath)
# boxplot(celfiles)
# hist(celfiles)
eset <- affy::rma(celfiles)
# Background correcting
# Normalizing
# Calculating Expression
exprData=exprs(eset)
colnames(exprData)=gsub("(GSM[0-9]+).*","\\1",colnames(exprData))
dim(exprData)
exprData[1:4,1:5]
range(exprData)
# [1]  2.387562 14.912034

class(exprData)
exprData=data.frame(exprData)
exprData$ID = rownames(exprData)

exprData = merge(exprData, probe, by = "ID")
head(exprData)
dim(exprData)
# [1] 42976   312

exprData = exprData %>% 
  dplyr::select(-ID) %>% 
  dplyr::select(Symbol, everything())
dim(exprData)
# [1] 42976   311

exprData[1:3,1:5]
table(duplicated(exprData$Symbol))

# FALSE  TRUE 
# 21655 21321

exprData <- exprData %>% 
  group_by(Symbol) %>%
  summarise_all(mean)
dim(exprData)
# [1] 21655   313

exprData=data.frame(exprData,check.names=F)
exprData[1:3,1:5]
rownames(exprData)=exprData$Symbol
exprData=exprData[,-1]
exprData[1:3,1:5]

gse71118.exp=exprData
dim(gse71118.exp)
# [1] 21655   311

save(gse71118.exp,file = 'raw_datas/Preprocessed/gse71118.exp.RData')
save(gse71118.cli,file = 'raw_datas/Preprocessed/gse71118.cli.RData')


setwd('Z:/users/Work1/Work/20240524_Liposarcomas')


library(tidyverse)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)


# # getGDCprojects()$project_id

# query <- GDCquery(
#   project = c("TCGA-SARC"),
#   data.category = "Transcriptome Profiling",
#   data.type = "Gene Expression Quantification",
#   workflow.type = "STAR - Counts"
# )

# GDCdownload(query)

# data = GDCprepare(query)

# saveRDS(data, "GDCdata/TCGA-SARC.rds")


data = readRDS("GDCdata/TCGA-SARC.rds")

############ 
assayNames(data)
# [1] "unstranded"       "stranded_first"   "stranded_second"  "tpm_unstrand"    
# [5] "fpkm_unstrand"    "fpkm_uq_unstrand"

exprs_fpkm <- assay(data, i = 'fpkm_unstrand') %>% as.data.frame()
dim(exprs_fpkm) ## 60660   177
。
gene_anno <- rowData(data)
# sample_infor <- colData(data)
# saveRDS(sample_infor,"GDCdata/TCGA-SARC_clinical.rds")

anno_info = gene_anno %>% 
  as.data.frame() %>% 
  subset(gene_type == "protein_coding") %>% 
  dplyr::select(gene_id, gene_name)
head(anno_info)
dim(anno_info) ## 19962*2

exprs_fpkm[1:4,1:5]
exprs_fpkm$gene_id = rownames(exprs_fpkm)

#################
expr = merge(exprs_fpkm, anno_info, by = "gene_id")
expr[1:4,1:4]

expr = expr %>% 
  dplyr::select(-gene_id) %>% 
  dplyr::select(gene_name, everything())
dim(expr) ## 19962   178
expr[1:3,1:5]
table(duplicated(expr$gene_name))

expr <- expr %>% 
  group_by(gene_name) %>%
  summarise_all(mean)
dim(expr) # 19938   178
expr=data.frame(expr,check.names=F)
expr[1:3,1:5]
rownames(expr)=expr$gene_name
expr=expr[,-1]
expr[1:3,1:5]
write.csv(expr, "GDCdata/TCGA-SARC_expr_fpkm.csv", row.names = T)
save(expr, file = "raw_datas/Preprocessed/tcga.expr.fpkm.RData")


library(dplyr)
# query <- GDCquery(project = "TCGA-SARC",
#                   data.category = "Clinical",
#                   data.type="Clinical Supplement",
#                   data.format = "BCR Biotab")
# GDCdownload(query)
# data = GDCprepare(query)
# saveRDS(data, "GDCdata/TCGA-SARC_Clinical.rds")
data=readRDS('GDCdata/TCGA-SARC_Clinical.rds')
patient_infor=data$clinical_patient_sarc
clinical_infor=data$clinical_follow_up_v4.0_sarc
colnames(patient_infor)=patient_infor[1,]
colnames(clinical_infor)=clinical_infor[1,]
patient_infor=patient_infor[-c(1,2),]
clinical_infor=clinical_infor[-c(1,2),]

#########################
head(patient_infor)
head(clinical_infor)
dim(patient_infor)
# [1] 261  68
dim(clinical_infor)
# [1] 327  14

all(patient_infor$bcr_patient_barcode==clinical_infor$bcr_patient_barcode)
# [1] FALSE

#####################################
clinical_infor=read.delim2('GDCdata/clinical.project-tcga-sarc.2024-05-27/clinical.tsv',stringsAsFactors = F)
dim(clinical_infor)
# [1] 522 210
# bcr_patient_barcode, days_to_last_followup,vital_status

surv_data = dplyr::select(clinical_infor, c(case_submitter_id, days_to_death,days_to_last_follow_up,vital_status)) %>% 
  subset(!is.na(days_to_death) & !is.na(vital_status)) %>%
  mutate(vital_status=case_when(vital_status=="Dead" ~ 1,
                                vital_status=="Alive" ~ 0,
                                TRUE ~ NA)) %>% data.frame()
head(surv_data)
surv_data[which(surv_data=="'--",arr.ind = T)]=NA
surv_data[which(surv_data=="[Not Available]",arr.ind = T)]=NA

mode(surv_data$days_to_death)="numeric"
mode(surv_data$days_to_last_follow_up)="numeric"

surv_data=surv_data %>% 
  rowwise() %>% 
  mutate(OS.time=max(days_to_death,days_to_last_follow_up,na.rm = T)) %>%
  select(case_submitter_id,OS.time,vital_status) %>% data.frame()

colnames(surv_data)=c('SampleID','OS.time','OS')
surv_data$OS.time[which(surv_data$OS.time=="'--")]=NA
surv_data$OS.time[which(surv_data$OS.time=="[Not Available]")]=NA
mode(surv_data$OS.time)="numeric"
mode(surv_data$OS)="numeric"

table(duplicated(surv_data$SampleID))

surv_data=surv_data %>%
  group_by(SampleID) %>%
  slice_max(OS.time,n=1) %>%
  ungroup() %>% data.frame()

surv_data=surv_data[!duplicated(surv_data),]
dim(surv_data)
# histological_type
# tumor_depth
# local_disease_recurrence
# metastatic_diagnosis
# tumor_tissue_site

# c(bcr_patient_barcode,age_at_initial_pathologic_diagnosis,gender,tumor_tissue_site,histological_type,tumor_depth,local_disease_recurrence,metastatic_diagnosis)

meta_data = dplyr::select(patient_infor, c(bcr_patient_barcode,age_at_initial_pathologic_diagnosis,gender,tumor_tissue_site,histological_type,tumor_depth,local_disease_recurrence,metastatic_diagnosis,person_neoplasm_cancer_status))
meta_data = data.frame(meta_data)
colnames(meta_data)=c('SampleID','age','gender','primary_site','histology','tumor_depth','local_disease_recurrence','metastatic_diagnosis','disease_free_status')


meta_data[which(meta_data=="[Not Available]",arr.ind = T)]=NA
meta_data[which(meta_data=="[Discrepancy]",arr.ind = T)]=NA
meta_data[which(meta_data=="NA",arr.ind = T)]=NA

unique(meta_data$histology)
names(table(meta_data$histology))

# [1] "Dedifferentiated liposarcoma"                                            :Liposarcoma:Dedifferentiated liposarcoma
# [2] "Desmoid Tumor"                                                           :Desmoid/Aggressive fibromatosis
# [3] "Giant cell 'MFH' / Undifferentiated pleomorphic sarcoma with giant cells":Undifferentiated Pleomorphic Sarcoma
# [4] "Leiomyosarcoma (LMS)"                                                    :Leiomyosarcoma
# [5] "Malignant Peripheral Nerve Sheath Tumors (MPNST)"                        :Malignant Peripheral Nerve Sheath Tumor
# [6] "Myxofibrosarcoma"                                                        :Myxofibrosarcoma
# [7] "Pleomorphic 'MFH' / Undifferentiated pleomorphic sarcoma"                :Undifferentiated Pleomorphic Sarcoma
# [8] "Sarcoma; synovial; poorly differentiated"                                :Synovial Sarcoma   
# [9] "Synovial Sarcoma - Biphasic"                                             :Synovial Sarcoma
# [10] "Synovial Sarcoma - Monophasic"                                          :Synovial Sarcoma
# [11] "Undifferentiated Pleomorphic Sarcoma (UPS)"                             :Undifferentiated Pleomorphic Sarcoma


meta_data=meta_data %>% dplyr::mutate(histology=case_when(histology=="Dedifferentiated liposarcoma" ~ "Liposarcoma:Dedifferentiated liposarcoma",
                                                          histology=="Desmoid Tumor" ~ "Desmoid/Aggressive fibromatosis",
                                                          histology %in% c("Giant cell 'MFH' / Undifferentiated pleomorphic sarcoma with giant cells","Pleomorphic 'MFH' / Undifferentiated pleomorphic sarcoma","Undifferentiated Pleomorphic Sarcoma (UPS)") ~ "Undifferentiated Pleomorphic Sarcoma",
                                                          histology=="Leiomyosarcoma (LMS)" ~ "Leiomyosarcoma",
                                                          histology=="Malignant Peripheral Nerve Sheath Tumors (MPNST)" ~ "Malignant Peripheral Nerve Sheath Tumor",
                                                          histology=="Myxofibrosarcoma" ~ "Myxofibrosarcoma",
                                                          histology %in% c("Sarcoma; synovial; poorly differentiated","Synovial Sarcoma - Biphasic","Synovial Sarcoma - Monophasic") ~ "Synovial sarcoma"))



setdiff(meta_data$SampleID,surv_data$SampleID)
setdiff(surv_data$SampleID,meta_data$SampleID)
dim(meta_data)
dim(surv_data)

tcga.cli=dplyr::left_join(meta_data,surv_data,by='SampleID')
tcga.cli=tcga.cli %>% select(!c(primary_site,local_disease_recurrence,metastatic_diagnosis,disease_free_status))
head(tcga.cli)
dim(tcga.cli) 
# [1] 261   7

write.csv(tcga.cli,file = 'raw_datas/Preprocessed/TCGA-SARC_Clinical.csv',row.names = F)
save(tcga.cli, file = "raw_datas/Preprocessed/tcga.cli.RData")


library(TCGAbiolinks)
query <- GDCquery(
  project = "TCGA-SARC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

GDCdownload(query)

GDCprepare(query, save = T,save.filename = "TCGA-SARC_SNP.Rdata")

query <- GDCquery(
  project = "TCGA-SARC",
  data.category = "DNA Methylation",
  platform = "Illumina human methylation 450",
  legacy = TRUE
)

GDCdownload(query)
data <- GDCprepare(query)


