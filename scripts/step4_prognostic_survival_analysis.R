dir.create('analysis/04_Prognostic_Model/',recursive = T)
################## GSE30929
# Well-differentiated liposarcoma: WDLPS
# Dedifferentiated liposarcoma: DDLPS
# Myxoid: MLPS
# Myxoid/round cell liposarcoma: MRCLS, MRCLS usually grows in the arms and legs. These tumors grow slowly, and they can spread to other parts of the body.
# Pleomorphic liposarcoma: PLPS
# myxoid pleomorphic liposarcoma: MPLPS
###################### TCGA
# dedifferentiated liposarcoma (DDLPS)
# leiomyosarcoma (arising in both gynecologic and soft tissue sites) (LMS)
# undifferentiated pleomorphic sarcoma (UPS)
# myxofibrosarcoma (MFS)
# malignant peripheral nerve sheath tumor (MPNST)

load("raw_datas/Preprocessed/gse30929_cli.RData")
load("raw_datas/Preprocessed/gse30929_exprs.RData")

load("raw_datas/Preprocessed/tcga_cli.RData")
load("raw_datas/Preprocessed/tcga_exprs.RData")


mode(gse30929_cli$DRFS)="numeric"
mode(gse30929_cli$DRFS.time)="numeric"

mode(tcga_cli$OS.time)="numeric"
mode(tcga_cli$OS)="numeric"

########################################################
wgcna.degs.intersect.genes=read.csv('analysis/02_WGCNA/wgcna.degs.intersect.genes.csv')
head(wgcna.degs.intersect.genes)
dim(wgcna.degs.intersect.genes) ## 101

###################
library(dplyr)
geo.immu.deg.sig=read.csv("analysis/03_Immune_Subtype/gse30929.immu.deg.sig.csv",row.names = 1)
dim(geo.immu.deg.sig)
# [1] 828   6

#################
geneset=read.table

genes_1=wgcna.degs.intersect.genes$x
genes_2=geo.immu.deg.sig %>% rownames()
genes_3=unique(geneset$Gene)
length(genes_1) ## 101
length(genes_2) ## 828
length(genes_3) ## 782

############################
geo.deg.sig=read.csv("analysis/01_DEGs/geo.deg.sig.csv",row.names = 1)
geo.deg.genes=rownames(geo.deg.sig)
head(geo.deg.genes)

turquoise_gene=read.csv(file = 'analysis/02_WGCNA/turquoise.genes.csv',stringsAsFactors = F)
turquoise_gene=as.character(turquoise_gene$x)
###################
myList=list(WGCNA=turquoise_gene,'DEGs (Tumor vs Normal)'=geo.deg.genes,`DEGs (Immunity_H vs Immunity_L)`=genes_2,`Immune-related signature`=genes_3)

mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

p=plot(eulerr::venn(myList),labels = list(col = "gray20", font = 2), 
       edges = list(col = "gray60", lex = 1),
       fills = list(fill = c("#11A579","#3969AC","#F2B701","#E73F74"), alpha = 1),
       quantities = list(cex = .8, col = 'gray20'))
p
savePDF('analysis/04_Prognostic_Model/genes.venn.pdf',p,height = 5,width = 6)

optimal_gene=Reduce(intersect,list(genes_1,genes_2,genes_3))
optimal_gene
length(optimal_gene) ## 10

#################################### GSE30929
exprs_data=gse30929_exprs
surv_data=gse30929_cli[,c('geo_accession','DRFS','DRFS.time')]
colnames(surv_data)=c('sample','status','time')
##########################
library(dplyr)
exprs_sub = subset(exprs_data, rownames(exprs_data) %in% c(optimal_gene)) %>% 
  t() %>% as.data.frame()
dim(exprs_sub) ## 140*10
exprs_sub=tibble::rownames_to_column(exprs_sub,var = 'sample')
head(exprs_sub)

cox_data=merge(surv_data,exprs_sub,by.x = 'sample',by.y ='sample',sort = F)
rownames(cox_data)=cox_data$sample
cox_data=cox_data[,-1]
cox_data=cox_data[which(cox_data$time>0),]
dim(cox_data) ## 140*12
cox_data[1:3,1:4]

cox.res=cox_batch(t(scale(cox_data[,optimal_gene]))
                   ,time = cox_data$time
                   ,event = cox_data$status)
dim(cox.res)
table(cox.res$p.value<0.05) ## 10
table(cox.res$p.value<0.01) ## 10

geo.cox.res=cox.res
geo_cox_data=cox_data


#################################### GSE30929
dim(cox_data)
train=data.frame(time=cox_data$time*365,status=cox_data$status,cox_data[,-c(1:2)])
train[1:4,1:5]
# RSF_genes=RSF.fs(train = train)
# CoxBoost_genes=CoxBoost.fs(train = train)
# stepwiseCox_genes=stepwiseCox.both.fs(train)
# Lasso_genes=Lasso.fs(train)
# genes=Reduce(intersect,list(Lasso_genes,RSF_genes))
# genes
# # [1] "CDKN3" "KNTC1" "PRC1"
# train_2=data.frame(time=cox_data$time*365,status=cox_data$status,cox_data[,genes])
# stepwiseCox_genes_2=stepwiseCox.both.fs(train_2)
# stepwiseCox_genes_2
# # [1] "KNTC1" "PRC1"

options(ggrepel.max.lihcerlaps = Inf)
mg_lasso_cox_use=function(dat,time,event,nfolds=3,lambda.min=T,show_text=T,figLabels=c('A','B')){
  library("glmnet") 
  library('survival')
  t.inds=which(!is.na(time)&!is.na(event)&time>0)
  dat=dat[t.inds,]
  time=as.numeric(time[t.inds])
  event=as.numeric(event[t.inds])
  y=Surv(time,event)
  set.seed(456789123)
  # set.seed(123456)
  fit1_cv = cv.glmnet(as.matrix(dat), y, family = "cox", nfolds=nfolds
                      ,nlambda=100, alpha=1,type.measure="mse")
  fit<-glmnet(dat, y, family = "cox")
  if(lambda.min){
    lambda=fit1_cv$lambda.min
  }else{
    lambda=fit1_cv$lambda.1se
  }
  coefficients<-coef(fit,s=lambda)
  Active.Index<-which(coefficients[,1]!=0)
  genes=row.names(coefficients)[Active.Index]
  Active.coefficients<-coefficients[Active.Index]  
  g=mg_plot_lasso(fit,fit1_cv,lambda = lambda,show_text=show_text,figLabels=figLabels)
  return(list(Mode1=fit,Model2=fit1_cv,Genes=genes,Coef=Active.coefficients,lambda=lambda,plot=g))
}

length(optimal_gene)
p.value=0.05
cox_data=geo_cox_data
dim(cox_data)
cox.sig=rownames(geo.cox.res)[which(geo.cox.res$p.value<p.value)]
length(cox.sig) ## 10

lasso.res=mg_lasso_cox_use((cox_data[,optimal_gene])
                           , time = cox_data$time
                           , event = cox_data$status
                           , nfolds = 10
                           , lambda.min = T
                           , figLabels = c('B', 'C'))
lasso.res$Genes
lasso.res$lambda

lasso.res$plot
savePDF('analysis/04_Prognostic_Model/gse30929.lasso.res.pdf',lasso.res$plot,height = 12,width = 6)

lst.modl=createCoxModel_use((cox_data[,lasso.res$Genes])
                            , time = cox_data$time
                            , event = cox_data$status
                            , isStep = T)
lst.modl$Cox
lst.modl$Genes
lst.modl$fmla

lst.modl.Coef=lst.modl$Coef
names(lst.modl.Coef)=lst.modl$Genes
lst.modl.Coef

gse30929.risk.score=lst.modl$Score
gse30929.risk.score=scale(gse30929.risk.score)
range(gse30929.risk.score)
length(gse30929.risk.score)

#######################
dt=data.frame(cox_data[,c('time','status')],cox_data[,lst.modl$Genes])
fmla <- as.formula(paste0("Surv(time, status) ~",paste0(lst.modl$Genes,collapse = '+')))
cox <- coxph(fmla, data =dt)
cox
pdf('analysis/04_Prognostic_Model/gse30929.model.ggforest.pdf', height = 5, width = 6,onefile = F)
survminer::ggforest(cox,data=dt,noDigits = 3)
dev.off()

################## ROC 
dt=data.frame(cox_data[,c('time','status')],RiskScore=gse30929.risk.score)
dim(dt)
head(dt)

gse30929.TimeROC=ggplotTimeROC_use(time = dt$time,status = dt$status,score = dt$RiskScore,mks = c(1,3,5))
gse30929.TimeROC
############
ROC=timeROC::timeROC(T=dt$time,
                     delta = dt$status,
                     marker = dt$RiskScore,
                     cause = 1, weighting = "marginal",
                     times = c(1,3,5),
                     iid = TRUE)
ROC

auc_data=data.frame(times=ROC$times,AUC=ROC$AUC,confint(ROC,level = 0.95,na.rm=T)$CI_AUC/100)
auc_data=round(auc_data,2)
auc_data
# times  AUC X2.5. X97.5.
# t=1     1 0.81  0.72   0.91
# t=3     3 0.80  0.71   0.89
# t=5     5 0.70  0.57   0.84

mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

lbs=paste0(paste0(auc_data$times,"-Years"),',AUC=',auc_data$AUC,',95%CI(',paste(auc_data$X2.5.,auc_data$X97.5., sep = "-"),")")

pdf("analysis/04_Prognostic_Model/gse30929.timeROC.pdf",height = 7,width = 7)
plot(ROC, time=1, col="#E73F74", add=FALSE, lwd=2, title = "")   
# plot(ROC, time=2 , col="#80BA5A", add=TRUE, lwd=2)    
plot(ROC, time=3 , col="#008695", add=TRUE, lwd=2)    
plot(ROC, time=5, col="#F2B701", add=TRUE, lwd=2)

legend("bottomright",
       legend=lbs,
       col=c("#E73F74", "#008695","#F2B701"),
       lty=1, lwd=2,bty = "n")
dev.off()

################
dt=data.frame(cox_data[,c('time','status')],RiskScore=gse30929.risk.score)
dt=dt %>% mutate(Group=ifelse(dt$RiskScore>0,"High","Low"))
rownames(dt)=rownames(cox_data)
head(dt)
writeMatrix(dt,outpath = 'analysis/04_Prognostic_Model/gse30929.group.txt')

cutoff <- survminer::surv_cutpoint(data.frame(time=dt$time,
                                              event = dt$status,
                                              risk = dt$RiskScore), time = "time", event = "event",
                                   variables = c("risk"))
cutoff=cutoff$cutpoint$cutpoint
cutoff

risk.group.color=c("#E73F74","#F2B701")
names(risk.group.color)=c('High','Low')

gse30929.km=ggplotKMCox(dt[,c("time","status","Group")]
                    , title = 'RiskType'
                    , pal = unname(risk.group.color)
                    , labs = c('High', 'Low')
                    , add_text = '')
gse30929.km



load("raw_datas/Preprocessed/gse71118.gse21050.expr_adjusted.RData")
load("raw_datas/Preprocessed/gse71118.gse21050.clinical_data.RData")  dedifferentiated liposarcoma

mode(clinical_data$DRFS.time)="numeric"
mode(clinical_data$DRFS)="numeric"

geo.model.dat=expr_adjusted[match(lst.modl$Genes,row.names(expr_adjusted)),]
dim(geo.model.dat)
# [1]   2 106

lst.vd.mod1=createCoxModel((t(geo.model.dat))
                           ,time=clinical_data$DRFS.time
                           ,event = clinical_data$DRFS)
geo.risk.score=lst.vd.mod1$Score
geo.risk.score=predictScoreByCoxModel(coxModel = lst.modl,dat =t(geo.model.dat))
geo.risk.score=scale(geo.risk.score)

lst.modl$fmla
# [1] "RiskScore=+0.79*KNTC1+0.366*PRC1"
lst.vd.mod1$fmla
# [1] "RiskScore=+0.311*KNTC1-0.038*PRC1"

dt=data.frame(clinical_data[,c('DRFS.time','DRFS')],RiskScore=geo.risk.score)
colnames(dt)[1:2]=c('time','status')
head(dt)

geo.TimeROC=ggplotTimeROC_use(time = dt$time,status = dt$status,score = dt$RiskScore,mks = c(1,3,5))
geo.TimeROC

######################## GSE21050
# GSE21050=getGEOExpData('GSE21050')
# gse21050.exprs=GSE21050$Exp$GPL570_54613_Data_col1
# range(gse21050.exprs)
# gse21050.exprs=log2(gse21050.exprs)
# gse21050.exprs[1:4,1:4]
# gse21050.exprs=exp_probe2symbol_v2(datExpr = gse21050.exprs,GPL = "GPL570")
# gse21050_exprs=gse21050.exprs %>% select(gse21050_cli$geo_accession)
load("raw_datas/Preprocessed/gse21050_exprs.RData")
load("raw_datas/Preprocessed/gse21050_cli.RData")
mode(gse21050_cli$DRFS.time)="numeric"
mode(gse21050_cli$DRFS)="numeric"

gse21050.model.dat=gse21050_exprs[match(lst.modl$Genes,row.names(gse21050_exprs)),]
dim(gse21050.model.dat)
# [1]   2 62

lst.vd.mod1=createCoxModel((t(gse21050.model.dat))
                           ,time=gse21050_cli$DRFS.time
                           ,event = gse21050_cli$DRFS)
gse21050.risk.score=lst.vd.mod1$Score
gse21050.risk.score=predictScoreByCoxModel(coxModel = lst.modl,dat =scale(t(gse21050.model.dat)))
gse21050.risk.score=scale(gse21050.risk.score)

lst.modl$fmla
# [1] "RiskScore=+0.79*KNTC1+0.366*PRC1"
lst.vd.mod1$fmla
# [1] "RiskScore=+0.296*KNTC1+0.06*PRC1"

dt=data.frame(gse21050_cli[,c('DRFS.time','DRFS')],RiskScore=gse21050.risk.score)
colnames(dt)[1:2]=c('time','status')
head(dt)

gse21050.TimeROC=ggplotTimeROC_use(time = dt$time,status = dt$status,score = dt$RiskScore,mks = c(1,3,5))
gse21050.TimeROC

######################## GSE71118
load("raw_datas/Preprocessed/gse71118_exprs.RData")
load("raw_datas/Preprocessed/gse71118_cli.RData")
mode(gse71118_cli$DRFS.time)="numeric"
mode(gse71118_cli$DRFS)="numeric"

gse71118.model.dat=gse71118_exprs[match(lst.modl$Genes,row.names(gse71118_exprs)),]
dim(gse71118.model.dat)
# [1]   2 44

lst.vd.mod1=createCoxModel((t(gse71118.model.dat))
                           ,time=gse71118_cli$DRFS.time
                           ,event = gse71118_cli$DRFS)
# gse71118.risk.score=lst.vd.mod1$Score
gse71118.risk.score=predictScoreByCoxModel(coxModel = lst.modl,dat =(t(gse71118.model.dat)))
gse71118.risk.score=scale(gse71118.risk.score)

lst.modl$fmla
# [1] "RiskScore=+0.79*KNTC1+0.366*PRC1"
lst.vd.mod1$fmla
# [1] "RiskScore=+0.291*KNTC1-0.201*PRC1"

dt=data.frame(gse71118_cli[,c('DRFS.time','DRFS')],RiskScore=gse71118.risk.score)
colnames(dt)[1:2]=c('time','status')
head(dt)

gse71118.TimeROC=ggplotTimeROC_use(time = dt$time,status = dt$status,score = dt$RiskScore,mks = c(1,3,5))
gse71118.TimeROC

############################# TCGA
load("raw_datas/Preprocessed/tcga_cli.RData")
load("raw_datas/Preprocessed/tcga_exprs.RData")
mode(tcga_cli$OS.time)="numeric"
mode(tcga_cli$OS)="numeric"


# tcga_cli_supp=getTCGAOSBySamples(samples = tcga_cli$SampleID)
# all(tcga_cli_supp$OS.time==tcga_cli$OS.time)
# dt=cbind(tcga_cli,tcga_cli_supp)
# dt[tcga_cli_supp$OS.time!=tcga_cli$OS.time,]
# tcga_cli=cbind(tcga_cli,tcga_cli_supp)
############
dim(tcga_cli)
dim(tcga_exprs)

tcga.model.dat=tcga_exprs[match(lst.modl$Genes,row.names(tcga_exprs)),]
dim(tcga.model.dat)
# [1]   2 58

lst.vd.mod2=createCoxModel(t(tcga.model.dat)
                           ,time=tcga_cli$OS.time
                           ,event = tcga_cli$OS)
# tcga.risk.score=lst.vd.mod2$Score
tcga.risk.score=predictScoreByCoxModel(coxModel = lst.modl,dat =scale(t(tcga.model.dat)))
tcga.risk.score=scale(tcga.risk.score)

lst.modl$fmla
# [1] "RiskScore=+0.79*KNTC1+0.366*PRC1"
lst.vd.mod2$fmla
# [1] "RiskScore=+0.352*KNTC1+0.254*PRC1"

dt=data.frame(tcga_cli[,c('OS.time','OS')],RiskScore=tcga.risk.score)
colnames(dt)[1:2]=c('time','status')
head(dt)

tcga.TimeROC=ggplotTimeROC_use(time = dt$time,status = dt$status,score = dt$RiskScore,mks = c(1,3,5))
tcga.TimeROC

############
ROC=timeROC::timeROC(T=dt$time/365,
                     delta = dt$status,
                     marker = dt$RiskScore,
                     cause = 1, weighting = "marginal",
                     times = c(1,3,5),
                     iid = TRUE)
ROC

auc_data=data.frame(times=ROC$times,AUC=ROC$AUC,confint(ROC,level = 0.95,na.rm=T)$CI_AUC/100)
auc_data=round(auc_data,2)
auc_data
# times  AUC X2.5. X97.5.
# t=1     1 0.71  0.48   0.94
# t=3     3 0.69  0.52   0.87
# t=5     5 0.79  0.63   0.94
mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

lbs=paste0(paste0(auc_data$times,"-Years"),',AUC=',auc_data$AUC,',95%CI(',paste(auc_data$X2.5.,auc_data$X97.5., sep = "-"),")")

pdf("analysis/04_Prognostic_Model/tcga.timeROC.pdf",height = 7,width = 7)
plot(ROC, time=1, col="#E73F74", add=FALSE, lwd=2, title = "")   
# plot(ROC, time=2 , col="#80BA5A", add=TRUE, lwd=2)    
plot(ROC, time=3 , col="#008695", add=TRUE, lwd=2)    
plot(ROC, time=5, col="#F2B701", add=TRUE, lwd=2)

legend("bottomright",
       legend=lbs,
       col=c("#E73F74", "#008695","#F2B701"),
       lty=1, lwd=2,bty = "n")
dev.off()

################
dt=data.frame(tcga_cli[,c('OS.time','OS')],RiskScore=tcga.risk.score)
dt=dt %>% mutate(Group=ifelse(dt$RiskScore>0,"High","Low"))
colnames(dt)[1:2]=c("time","status")
dt$time=dt$time/365
head(dt)

writeMatrix(dt,outpath = 'analysis/04_Prognostic_Model/tcga.group.txt')

cutoff <- survminer::surv_cutpoint(data.frame(time=dt$time/365,
                                              event = dt$status,
                                              risk = dt$RiskScore), time = "time", event = "event",
                                   variables = c("risk"))
cutoff=cutoff$cutpoint$cutpoint
cutoff

tcga.km=ggplotKMCox(dt[,c("time","status","Group")]
                    , title = 'RiskType'
                    , pal = unname(risk.group.color)
                    , labs = c('High', 'Low')
                    , add_text = '')
tcga.km

###########################
plot1=mg_merge_plot(gse30929.TimeROC,tcga.TimeROC,nrow = 1,ncol = 2,labels = c("GSE30929",'TCGA'))
plot2=mg_merge_plot(gse30929.km,tcga.km,nrow = 1,ncol = 2,labels = c("GSE30929",'TCGA'))

plot=mg_merge_plot(gse30929.km,tcga.km,
                   nrow = 1,ncol = 2,labels = c("GSE30929",'TCGA'))
plot

savePDF('analysis/04_Prognostic_Model/model.timeROC.pdf',plot1,height = 6,width = 12)
savePDF('analysis/04_Prognostic_Model/model.km.pdf',plot2,height = 6,width = 12)


prognostic.model=data.frame(gene=lst.modl$Genes,coef=lst.modl$Coef)
rownames(prognostic.model)=NULL
head(prognostic.model)
write.csv(prognostic.model,file = 'analysis/04_Prognostic_Model/prognostic.model.csv')


plot_sankey=function(df_m){
  library(ggalluvial)
  library(ggplot2)
  library(dplyr)
  corLodes=to_lodes_form(df_m, axes = 1:ncol(df_m), id = "Cohort")
  mycol <- rep(c("#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
  # mycol=c(pal_npg()(9)[c(1,2)],pal_lancet()(9)[c(7,4,3)])
  mycol=mycolors
  ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
    scale_x_discrete(expand = c(0, 0)) +  

    geom_flow(width = 2/10,aes.flow = "forward") + 
    geom_stratum(alpha = .9,width = 2/10) +
    scale_fill_manual(values = mycol) +
    #size =
    geom_text(stat = "stratum", size = 2,color="black") +
    xlab("") + ylab("") + theme_bw() + 
    theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
    theme(panel.grid =element_blank()) + 
    theme(panel.border = element_blank()) + 
    ggtitle("") + guides(fill = FALSE)   
}

clin.color=mycolors[1:5]

####################### GSE30929
load("raw_datas/Preprocessed/gse30929_cli.RData")
gse30929.group=readMatrix('analysis/04_Prognostic_Model/gse30929.group.txt')
gse30929.subtype=readMatrix('analysis/03_Immune_Subtype/gse30929.subtype.txt',header = F)
all(gse30929_cli$geo_accession==rownames(gse30929.group))
all(gse30929_cli$geo_accession==rownames(gse30929.subtype))
all(rownames(gse30929.group)==rownames(gse30929.subtype))

dt=data.frame(gse30929_cli,RiskScore=gse30929.group$RiskScore,RiskType=gse30929.group$Group,Cluster=gse30929.subtype$V2)
# Distant recurrence-free survival (DRFS)
dt$Metastasis=NA
dt$Metastasis[which(dt$DRFS==1)]="Yes"
dt$Metastasis[which(dt$DRFS==0)]="No"
dt$Metastasis=factor(dt$Metastasis,levels = c("Yes","No"))
head(dt)
# "Dedifferentiated liposarcoma","Myxoid/round cell liposarcoma", "Myxoid","Pleomorphic liposarcoma","Well-differentiated liposarcoma"
ordered=c("Pleomorphic liposarcoma","Dedifferentiated liposarcoma","Myxoid/round cell liposarcoma","Myxoid","Well-differentiated liposarcoma")
dt$histology=factor(dt$histology,levels =ordered)
head(dt)

p.all=list()
for(i in c('histology','Metastasis','Cluster')){
  print(i)
  dt1=dt[,c(i,'RiskScore')]
  colnames(dt1)=c("feature","RiskScore")
  dt1=data.frame(dt1)
  dt1=na.omit(dt1)
  dt1=na.omit(dt1)
  head(dt1)
  library(gghalves)
  p=ggplot(data = dt1, aes(x = feature, y = RiskScore, fill = feature)) +
    geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
    geom_point(aes(y = RiskScore, color = feature), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
    labs(y="RiskScore",x = stringr::str_to_title(i)) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values =clin.color ) +
    scale_colour_manual(values = clin.color) +
    theme_niwot()
  p
  cmp=data.frame(combn(c(1:length(unique(dt1$feature))),2))
  cmp
  p=p+ggpubr::stat_compare_means(comparisons = as.list(cmp),method = 'wilcox.test',label= "p.signif")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
p.all[[1]]=p.all[[1]]+theme(axis.text.x = element_text(angle = 30, hjust = 1))

head(dt)
p=plot_sankey(dt[,c('RiskType','Cluster',"histology","Metastasis")])
p.all=c(p.all,list(p))
length(p.all)

plot1=mg_merge_plot(p.all,nrow = 2,ncol = 2,labels = LETTERS[1:4])
plot1
savePDF('analysis/04_Prognostic_Model/gse30929.model.clinical.pdf',plot1,height = 10,width = 10)

####################### TCGA
load("raw_datas/Preprocessed/tcga_cli.RData")
tcga.group=readMatrix('analysis/04_Prognostic_Model/tcga.group.txt')
all(tcga_cli$SampleID==rownames(tcga.group))

dt=data.frame(tcga_cli,RiskScore=tcga.group$RiskScore,RiskType=tcga.group$Group)
dt$Status=NA
dt$Status[which(dt$OS==1)]="Dead"
dt$Status[which(dt$OS==0)]="Alive"
dt$Status=factor(dt$Status,levels = c("Dead","Alive"))

mode(dt$age)="numeric"
dt$Age=ifelse(dt$age>60,">60","<=60")
head(dt)

p.all=list()
for(i in c('Age','gender','Status')){
  print(i)
  dt1=dt[,c(i,'RiskScore')]
  colnames(dt1)=c("feature","RiskScore")
  dt1=data.frame(dt1)
  dt1=na.omit(dt1)
  dt1=na.omit(dt1)
  head(dt1)
  library(gghalves)
  p=ggplot(data = dt1, aes(x = feature, y = RiskScore, fill = feature)) +
    geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
    geom_point(aes(y = RiskScore, color = feature), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
    labs(y="RiskScore",x = stringr::str_to_title(i)) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values =clin.color ) +
    scale_colour_manual(values = clin.color) +
    theme_niwot()
  p
  cmp=data.frame(combn(c(1:length(unique(dt1$feature))),2))
  cmp
  p=p+ggpubr::stat_compare_means(comparisons = as.list(cmp),method = 't.test',label= "p.format")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
plot2=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = LETTERS[1])
plot2

cor.test(x=dt$age,y = dt$RiskScore,method = 'spearman')
# rho 
# -0.01050476


load("raw_datas/TME/gse30929.immu.ssgsea_28.RData")
geo.immu.ssgsea=t(rbind(gse30929.immu.ssgsea_28))
all(rownames(geo.immu.ssgsea)==gse30929_cli$geo_accession)

load("raw_datas/TME/tcga.immu.ssgsea_28.RData")
tcga.immu.ssgsea=t(rbind(tcga.immu.ssgsea_28))
all(rownames(tcga.immu.ssgsea)==tcga_cli$SampleID)


p1=mg_PlotMutiBoxplot(geo.immu.ssgsea
                      , group = gse30929.group$Group
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , fill = T
                      , binwidth = 0.1
                      , group_cols = risk.group.color
                      , test_method = 't.test') + labs(fill = "RiskType")

p1

p2=mg_PlotMutiBoxplot(tcga.immu.ssgsea
                      , group = tcga.group$Group
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "ssGSEA score"
                      , fill = T
                      , binwidth = 0.1
                      , group_cols = risk.group.color
                      , test_method = 't.test') + labs(fill = "RiskType")

p2

plot1=mg_merge_plot(p1,p2,nrow = 2,ncol = 1,labels = LETTERS[1:2])
plot1

############################# GSE30929
load("raw_datas/TME/gse30929.exp.estimate.RData")
rownames(gse30929.exp.estimate)=gse30929.exp.estimate$ID
gse30929.exp.estimate=gse30929.exp.estimate[,-1]
colnames(gse30929.exp.estimate)=gsub("_estimate","",colnames(gse30929.exp.estimate))

all(rownames(gse30929.group)==gse30929_cli$geo_accession)
gse30929_cli_all=data.frame(gse30929_cli,Groups=gse30929.group$Group,stringsAsFactors = F)

dt=cbind(gse30929_cli_all,gse30929.exp.estimate[gse30929_cli_all$geo_accession,])
head(dt)

selected=colnames(gse30929.exp.estimate)[-3]
p.all=list()
for(i in selected){
  print(i)
  dt1=dt[,c('Groups',i)]
  colnames(dt1)=c("feature","RiskScore")
  dt1=data.frame(dt1)
  dt1=na.omit(dt1)
  dt1=na.omit(dt1)
  library(gghalves)
  p=ggplot(data = dt1, aes(x = feature, y = RiskScore, fill = feature)) +
    geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
    geom_point(aes(y = RiskScore, color = feature), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
    labs(x="",y = stringr::str_to_title(i)) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values =risk.group.color ) +
    scale_colour_manual(values = risk.group.color) +
    theme_niwot()
  p
  cmp=data.frame(combn(c(1:length(unique(dt1$feature))),2))
  cmp
  p=p+ggpubr::stat_compare_means(comparisons = as.list(cmp),method = 'wilcox.test',label= "p.signif")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
plot2=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = LETTERS[3])
plot2

############################# TCGA
load("raw_datas/TME/tcga.exp.estimate.RData")
rownames(tcga.exp.estimate)=tcga.exp.estimate$ID
tcga.exp.estimate=tcga.exp.estimate[,-1]
colnames(tcga.exp.estimate)=gsub("_estimate","",colnames(tcga.exp.estimate))

all(rownames(tcga.group)==tcga_cli$SampleID)
tcga_cli_all=data.frame(tcga_cli,Groups=tcga.group$Group,stringsAsFactors = F)

dt=cbind(tcga_cli_all,tcga.exp.estimate[tcga_cli_all$SampleID,])
head(dt)

selected=colnames(tcga.exp.estimate)[-3]
p.all=list()
for(i in selected){
  print(i)
  dt1=dt[,c('Groups',i)]
  colnames(dt1)=c("feature","RiskScore")
  dt1=data.frame(dt1)
  dt1=na.omit(dt1)
  dt1=na.omit(dt1)
  library(gghalves)
  p=ggplot(data = dt1, aes(x = feature, y = RiskScore, fill = feature)) +
    geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
    geom_point(aes(y = RiskScore, color = feature), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
    labs(x="",y = stringr::str_to_title(i)) +
    guides(fill = FALSE, color = FALSE) +
    scale_fill_manual(values =risk.group.color ) +
    scale_colour_manual(values = risk.group.color) +
    theme_niwot()
  p
  cmp=data.frame(combn(c(1:length(unique(dt1$feature))),2))
  cmp
  p=p+ggpubr::stat_compare_means(comparisons = as.list(cmp),method = 'wilcox.test',label= "p.signif")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
plot3=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = LETTERS[4])
plot3

plot23=mg_merge_plot(plot2,plot3,nrow = 2,ncol = 1)
plot23

figure=mg_merge_plot(plot1,plot23,nrow = 1,ncol = 2,widths = c(1.5,1))
figure
savePDF('analysis/04_Prognostic_Model/model.immu.boxplot.pdf',figure,height = 10,width = 12)

gse30929.tme=cbind(geo.immu.ssgsea,gse30929.exp.estimate[,-3])
gse30929.tme[1:4,1:5]

cr=psych::corr.test(y=data.frame(gse30929.tme[rownames(gse30929.group),],check.names = F),
                    x=data.frame(t(gse30929_exprs)[rownames(gse30929.group),lst.modl$Genes],RiskScore=gse30929.group$RiskScore),
                    method = 'spearman')

df_cor=cr$r
df_pval=cr$p.adj
df_cor=round(df_cor,2)
head(df_cor)

### Pivot data from wide to long
library(tidyverse)
g = pivot_longer(data=rownames_to_column(as.data.frame(df_cor),var = "from"),
                 cols = 2:(ncol(df_cor)+1), ## Columns to pivot into longer format
                 names_to = "to",
                 values_to = "cor")
gp = pivot_longer(data=rownames_to_column(as.data.frame(df_pval)),
                  cols = 2:(ncol(df_pval)+1),
                  names_to = "gene",
                  values_to = "p")
all(g$from==gp$rowname & g$to==gp$gene)
g$p.adj = gp$p

###################
df=g
df <- df %>%
  mutate(col = cut(cor, breaks = c(-1, 0, 1),
                   labels = c("negative", "positive")),
         p.signif = cut(p.adj, breaks = c(0,0.0001, 0.001, 0.01, 0.05,1),
                        labels = c("****", "**", "**","*",""),
                        right = FALSE, include.lowest = TRUE))
df=data.frame(df)
head(df)
mode(df$p.adj)="numeric"
mode(df$cor)="numeric"

# df=df[which(df$p.adj<0.05 & abs(df$cor)>0.4),]
length(unique(df$to))
head(df)
writeMatrix(df,outpath = 'analysis/04_Prognostic_Model/gse30929.model.immu.cor.txt')

corr.mat=pivot_wider(df[,c(1,2,3)],names_from  ="from",values_from ='cor')
p.mat=pivot_wider(df[,c(1,2,4)],names_from  ="from",values_from ='p.adj')
corr.mat=data.frame(corr.mat)
p.mat=data.frame(p.mat)
head(corr.mat)

rownames(corr.mat)=corr.mat$to
corr.mat=corr.mat[-1]
rownames(p.mat)=p.mat$to
p.mat=p.mat[-1]
dim(corr.mat)
# [1] 31  3
####################
library(corrplot)
pdf('analysis/04_Prognostic_Model/model.TME.corrplot.pdf',height = 10,width = 12,onefile = F)
corrplot(corr = as.matrix(t(corr.mat)),
         p.mat = as.matrix(t(p.mat)),
         mar = c(0,0,0,0),
         col=colorRampPalette(c('#3969AC', 'white','#E73F74'))(50),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.8,cl.cex = 0.8,
         addgrid.col = 'white',
         method = "pie",
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()

############################## GSEA analysis
R version 4.3.2


T.exhausted=c("PDCD1","CD274","LAG3","HAVCR2","TIGIT")
setdiff(T.exhausted,rownames(gse30929_exprs))
rownames(gse30929_exprs)[grep("^WUCAM",rownames(gse30929_exprs))]

dt=t(gse30929_exprs[T.exhausted,rownames(gse30929.group)])
dt=data.frame(dt)
dt=reshape::rename(dt,c("PDCD1"="PD-1","CD274"="PD-L1","LAG3"="LAG-3","HAVCR2"="TIM-3"))
head(dt)

mg_PlotMutiBoxplot(dt
                   , group = gse30929.group$Group
                   , legend.pos = "top"
                   , add = 'boxplot'
                   , ylab = "log2(TPM+1)"
                   , fill = T
                   , binwidth = 0.1
                   , group_cols = risk.group.color
                   , test_method = 'kruskal.test') + 
  labs(title = "T cell exhaustion markers",fill = "RiskType")

#################### TCGA
setdiff(T.exhausted,rownames(tcga_exprs))
dt=t(tcga_exprs[T.exhausted,rownames(tcga.group)])
dt=data.frame(dt)
dt=reshape::rename(dt,c("PDCD1"="PD-1","CD274"="PD-L1","LAG3"="LAG-3","HAVCR2"="TIM-3"))
head(dt)

mg_PlotMutiBoxplot(dt
                   , group = tcga.group$Group
                   , legend.pos = "top"
                   , add = 'boxplot'
                   , ylab = "log2(TPM+1)"
                   , fill = T
                   , binwidth = 0.1
                   , group_cols = risk.group.color
                   , test_method = 'kruskal.test') + 
  labs(title = "T cell exhaustion markers",fill = "RiskType")

################# IPS
library(IOBR)
############################## GSE30929
gse30929.ips=deconvo_ips(eset = gse30929_exprs,project = 'GSE30929',plot=F)
gse30929.ips=cbind(gse30929.ips,RiskType=gse30929.group$Group)
rownames(gse30929.ips)=gse30929.ips$ID
dim(gse30929.ips)
head(gse30929.ips)
# writeMatrix(tcga.ips,outpath = "analysis/05_Model_TME/tcga.ips.txt")

dt=gse30929.ips[,c('RiskType','IPS_IPS')]
colnames(dt)[2]='feature'

p=ggplot(data = dt, aes(x = RiskType, y = feature, fill = RiskType)) +
  geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
  geom_point(aes(y = feature, color = RiskType), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
  labs(y = paste('Immune Phenotype Score',""), x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = risk.group.color) +
  scale_colour_manual(values = risk.group.color) +
  theme_niwot()
p=p+ggpubr::stat_compare_means(comparisons = list(c(1,2)),method = 'wilcox.test',label= "p.signif")
p

############################## TCGA
tcga.ips=deconvo_ips(eset = tcga_exprs,project = 'TCGA-SARC',plot=F)
tcga.ips=cbind(tcga.ips,RiskType=tcga.group$Group)
rownames(tcga.ips)=tcga.ips$ID
dim(tcga.ips)
head(tcga.ips)
# writeMatrix(tcga.ips,outpath = "analysis/05_Model_TME/tcga.ips.txt")

dt=tcga.ips[,c('RiskType','IPS_IPS')]
colnames(dt)[2]='feature'

p=ggplot(data = dt, aes(x = RiskType, y = feature, fill = RiskType)) +
  geom_half_violin(side='R',position = position_nudge(x = 0.2, y = 0), alpha = 0.6) +
  geom_point(aes(y = feature, color = RiskType), position = position_jitter(width = 0.15), size = 1, alpha = 0.6) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.6) +
  labs(y = paste('Immune Phenotype Score',""), x = NULL) +
  guides(fill = FALSE, color = FALSE) +
  scale_fill_manual(values = risk.group.color) +
  scale_colour_manual(values = risk.group.color) +
  theme_niwot()
p=p+ggpubr::stat_compare_means(comparisons = list(c(1,2)),method = 'wilcox.test',label= "p.signif")
p

##################### IMvigor210CoreBiologies ######
library("IMvigor210CoreBiologies")
data(cds)

imvigor.cli<-pData(cds)
imvigor.cli=rownames_to_column(imvigor.cli,var = 'SampleID')
imvigor.cli=imvigor.cli[,c('SampleID','os','censOS','binaryResponse')]
colnames(imvigor.cli)=c('SampleID','OS.time','OS','Response')
imvigor.cli$OS.time=imvigor.cli$OS.time*30
rownames(imvigor.cli)=imvigor.cli$SampleID

expr=mg_get_immu_pd1_treament_exp()
imvigor.exp=expr$fpkm
range(imvigor.exp)
imvigor.exp<-log2(imvigor.exp+1)
dim(imvigor.exp)
range(imvigor.exp)

#########
setdiff(lst.modl$Genes,rownames(imvigor.exp))
imvigor.model.dat=imvigor.exp[lst.modl$Genes,]
dim(imvigor.model.dat)
range(imvigor.model.dat)

all(colnames(imvigor.exp)==rownames(imvigor.cli))

vd.mod1=createCoxModel((t(imvigor.model.dat))
                       ,time=imvigor.cli$OS.time/365
                       ,event = imvigor.cli$OS)

imvigor.risk.score=vd.mod1$Score
# imvigor.risk.score=predictScoreByCoxModel(coxModel = lst.modl,dat = t(imvigor.model.dat))
imvigor.risk.score=scale(imvigor.risk.score)

lst.modl$fmla
# [1] "RiskScore=+0.79*KNTC1+0.366*PRC1"
vd.mod1$fmla
# [1] "RiskScore=-0.3*KNTC1+0.028*PRC1"

dt=data.frame(imvigor.cli[,c('OS.time','OS',"Response")],RiskScore=imvigor.risk.score,stringsAsFactors = F)
dt=dt %>% mutate(Group=ifelse(dt$RiskScore>0,"High","Low"))
dt$OS.time=dt$OS.time/365
head(dt)

imvigor.km=ggplotKMCox(dt[,c("OS.time","OS","Group")]
                       , title = 'RiskType'
                       , pal = unname(risk.group.color)
                       , labs = c('High', 'Low')
                       , add_text = '')
imvigor.km

results1=prop.table(table(dt$Response,dt$Group),margin=2)
results1=reshape2::melt(results1)
colnames(results1)<-c("type","Risk","Percentage")
results1$Percentage<-round(results1$Percentage,digits=2)
results1

ggplot(results1,aes(x=Risk,y=Percentage,fill=type))+
  geom_bar(position = "fill",stat="identity")+
  scale_fill_manual(values = pal_npg(palette = 'nrc',alpha = 0.9)(9)[c(1,2,3)])+
  theme_bw()+theme(legend.position = 'top')+
  geom_text(aes(label = Percentage),position=position_stack(vjust =0.5),size = 5)+
  labs(fill="Response")

mg_violin(data = na.omit(dt[,c('Response','RiskScore')]),melt=T)+labs(y="RiskScore")

########################### GSE78220
gse78220.cli=getGEOSampleData('GSE78220')
table(gse78220.cli$`vital status`)

gse78220.cli=gse78220.cli[,c('Acc','Title','overall survival (days)','vital status','anti-pd-1 response','Desc')]
colnames(gse78220.cli)=c('SampleID','Patient','OS.time','OS','Response','Type')
gse78220.cli$OS[which(gse78220.cli$OS=='Alive')]=0
gse78220.cli$OS[which(gse78220.cli$OS=='Dead')]=1
gse78220.cli$Type=stringr::str_split(gse78220.cli$Type,pattern = ", ",simplify =T)[,2]
table(gse78220.cli$Response)
gse78220.cli$Response[which(gse78220.cli$Response=="Complete Response")]='CR'
gse78220.cli$Response[which(gse78220.cli$Response=="Partial Response")]='PR'
gse78220.cli$Response[which(gse78220.cli$Response=="Progressive Disease")]='PD'
head(gse78220.cli)

table(gse78220.cli$Type)
gse78220.cli=gse78220.cli[which(gse78220.cli$Type=="pre anti-PD-1 treatment"),]
gse78220.cli=gse78220.cli[which(gse78220.cli$OS.time>0),]
mode(gse78220.cli$OS.time)="numeric"
mode(gse78220.cli$OS)="numeric"
dim(gse78220.cli)

gse78220.exp=readxl::read_xlsx('raw_datas/GSE78220_PatientFPKM.xlsx')
gse78220.exp=data.frame(gse78220.exp)
rownames(gse78220.exp)=gse78220.exp[,1]
gse78220.exp=gse78220.exp[,-1]
colnames(gse78220.exp)=gsub("(.*?)\\..*","\\1",colnames(gse78220.exp))
colnames(gse78220.exp)=gse78220.cli$SampleID[match(colnames(gse78220.exp),gse78220.cli$Patient)]
gse78220.exp[1:4,1:5]
range(gse78220.exp)
gse78220.exp=log2(gse78220.exp+1)
dim(gse78220.exp)

gse78220.exp=gse78220.exp[,gse78220.cli$SampleID]
match(lst.modl$Genes,rownames(gse78220.exp))

gse78220.model.dat=gse78220.exp[match(lst.modl$Genes,row.names(gse78220.exp)),]
dim(gse78220.model.dat)
# [1]   4 26
range(gse78220.model.dat)

vd.mod2=createCoxModel(scale(t(gse78220.model.dat))
                       ,time=gse78220.cli$OS.time/365
                       ,event = gse78220.cli$OS)
gse78220.risk.score=vd.mod2$Score
gse78220.risk.score=scale(gse78220.risk.score)
# gse78220.risk.score=predictScoreByCoxModel(coxModel = lst.modl,dat =t(gse78220.model.dat))
# gse78220.risk.score=scale(gse78220.risk.score)

vd.mod2$fmla
# "RiskScore=+0.348*KNTC1-0.197*PRC1"

dt=data.frame(gse78220.cli[,c('OS.time','OS',"Response")],RiskScore=gse78220.risk.score,stringsAsFactors = F)
dt=dt %>% mutate(Group=ifelse(dt$RiskScore>0,"High","Low"))
dt$OS.time=dt$OS.time/365
head(dt)

dt=data.frame(imvigor.cli[,c('OS.time','OS',"Response")],RiskScore=imvigor.risk.score,stringsAsFactors = F)
dt=dt %>% mutate(Group=ifelse(dt$RiskScore>0,"High","Low"))
dt$OS.time=dt$OS.time/365
head(dt)

gse78220.km=ggplotKMCox(dt[,c("OS.time","OS","Group")]
                       , title = 'RiskType'
                       , pal = unname(risk.group.color)
                       , labs = c('High', 'Low')
                       , add_text = '')
gse78220.km

p2=ggplotTimeROC_use(time = dt$OS.time,status = dt$OS,score = dt$RiskScore,mks = c(0.5,1,2))
p2

mg_violin(data = na.omit(dt[,c('Response','RiskScore')]),melt=T)

