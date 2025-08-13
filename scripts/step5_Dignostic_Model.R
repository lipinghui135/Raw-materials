
load("raw_datas/Preprocessed/gse21122_cli.RData")
load("raw_datas/Preprocessed/gse21122_exprs.RData")
load("raw_datas/Preprocessed/gse159659_cli.RData")
load("raw_datas/Preprocessed/gse159659_exprs.RData")
load("raw_datas/Preprocessed/gse21122.gse159659.expr_adjusted.RData")
load("raw_datas/Preprocessed/gse21122.gse159659.clinical_data.RData")

######################################
dir.create('analysis/05_Dignostic_Model/')
prognostic.model=read.csv(file = 'analysis/04_Prognostic_Model/prognostic.model.csv',row.names = 1,stringsAsFactors = F)
prognostic.model
# gene      coef
# 1 KNTC1 0.7901395
# 2  PRC1 0.3661083

final.genes=prognostic.model$gene
final.genes
# [1] "KNTC1" "PRC1" 
final.genes=c("KNTC1","PRC1")

exprs_1=gse21122_exprs
exprs_1=data.frame(t(exprs_1[final.genes,]))
exprs_1$Group=gse21122_cli$sample_type
head(exprs_1)

exprs_2=gse159659_exprs
exprs_2=data.frame(t(exprs_2[final.genes,]))
exprs_2$Group=gse159659_cli$sample_type
head(exprs_2)

exprs_3=expr_adjusted
exprs_3=data.frame(t(exprs_3[final.genes,]))
exprs_3$Group=clinical_data$sample_type
head(exprs_3)

mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

sampleType.color=c("#11A579","#E73F74")
# cli.color=c("#11A579","#3969AC","#F2B701","#E73F74")

p1=mg_PlotMutiBoxplot(exprs_1[,final.genes]
                      , group = exprs_1$Group
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "log2(TPM+1)"
                      , fill = T
                      , binwidth = 0.1
                      , group_cols = sampleType.color
                      , test_method = 't.test') + labs(fill = "Tissue",title = 'GSE21122')

p1

p2=mg_PlotMutiBoxplot(exprs_2[,final.genes]
                      , group = exprs_2$Group
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "log2(TPM+1)"
                      , fill = T
                      , binwidth = 0.1
                      , group_cols = sampleType.color
                      , test_method = 't.test') + labs(fill = "Tissue",title = 'GSE159659')

p2

p3=mg_PlotMutiBoxplot(exprs_3[,final.genes]
                      , group = exprs_3$Group
                      , legend.pos = "top"
                      , add = 'boxplot'
                      , ylab = "Expresion"
                      , fill = T
                      , binwidth = 0.1
                      , group_cols = sampleType.color
                      , test_method = 't.test') + labs(fill = "Tissue",title = 'GSE21122 and GSE159659)')

p3

plot1=mg_merge_plot(p1,p2,p3,nrow = 1,ncol = 3,labels = 'A',common.legend = T)
plot1

############################### ROC
model_data=data.frame(tissue=exprs_1$Group,exprs_1[,final.genes])
dim(model_data)

test_data_1=data.frame(tissue=exprs_2$Group,exprs_2[,final.genes])
dim(test_data_1)

test_data_2=data.frame(tissue=exprs_3$Group,exprs_3[,final.genes])
dim(test_data_2)


gene.colors=mycolors[1:5]
names(gene.colors)=final.genes

roc.list <- roc(tissue ~ ., data = model_data)

auc.res=sapply(roc.list,function(x){round(x$auc,3)},simplify=T)
auc.res
data.labels=paste0(names(auc.res)," , AUC = ",paste(round(auc.res,2)))

p1=ggroc(roc.list, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=8)+
  scale_colour_manual(labels=data.labels,values = gene.colors)+
  theme_bw()+labs(title = "GSE21122")+
  theme(legend.position = c(0.8, 0.2))
p1

roc.list <- roc(tissue ~ ., data = test_data_1)
auc.res=sapply(roc.list,function(x){round(x$auc,3)},simplify=T)
auc.res
data.labels=paste0(names(auc.res)," , AUC = ",paste(round(auc.res,2)))

p2=ggroc(roc.list, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=8)+
  scale_colour_manual(labels=data.labels,values = gene.colors)+
  theme_bw()+labs(title = "GSE159659")+
  theme(legend.position = c(0.8, 0.2))
p2

roc.list <- roc(tissue ~ ., data = test_data_2)
auc.res=sapply(roc.list,function(x){round(x$auc,3)},simplify=T)
auc.res
data.labels=paste0(names(auc.res)," , AUC = ",paste(round(auc.res,2)))

p3=ggroc(roc.list, legacy.axes = TRUE )+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               color="darkgrey", linetype=8)+
  scale_colour_manual(labels=data.labels,values = gene.colors)+
  theme_bw()+labs(title = "GSE21122 and GSE159659")+
  theme(legend.position = c(0.8, 0.2))
p3

plot2=mg_merge_plot(p1,p2,p3,nrow = 1,ncol = 3,labels = 'B')
plot2

genes.ggroc=mg_merge_plot(plot1,plot2,nrow = 2,ncol = 1)
genes.ggroc
savePDF('analysis/05_Dignostic_Model/model.genes.ggroc.pdf',genes.ggroc,height = 10,width = 15)

# library(pROC)
# library(e1071)
# paste0(final.genes,collapse = '+')
# #########
# model =svm(tissue ~ . ,model_data, probability=TRUE)
# summary(model)

################### GSE21122 
fmla <- as.formula(paste0("tissue ~ ",paste(final.genes , collapse= "+")))
model <- glm(fmla, data=model_data, family=binomial)
summary(model)

########## ICGC
pred<-predict(object =model,newdata = model_data[,final.genes])
pred <- as.ordered(pred)
model.roc <- roc(model_data$tissue,pred)

pdf('analysis/05_Dignostic_Model/gse21122.roc.pdf',height =5,width = 5)
plot(model.roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), 
     grid.col=c("blue", "red"),
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE,main='GSE21122')
dev.off()

########## GSE159659
pred<-predict(object =model,newdata = test_data_1[,final.genes])
pred <- as.ordered(pred)
model.roc <- roc(test_data_1$tissue,pred)

pdf('analysis/05_Dignostic_Model/gse159659.roc.pdf',height =5,width = 5)
plot(model.roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), 
     grid.col=c("blue", "red"),
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE,main='GSE159659')
dev.off()

########## GSE21122 and GSE159659
pred<-predict(object =model,newdata = test_data_2[,final.genes])
pred <- as.ordered(pred)
model.roc <- roc(test_data_2$tissue,pred)

pdf('analysis/05_Dignostic_Model/geo.roc.pdf',height =5,width = 5)
plot(model.roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.5, 0.2), 
     grid.col=c("blue", "red"),
     max.auc.polygon=T, auc.polygon.col="skyblue", print.thres=TRUE,main='GSE21122 and GSE159659')
dev.off()

save.image(file = '20240524_Liposarcomas_final.RData')



