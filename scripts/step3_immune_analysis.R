source('scripts/my_funcitions.R')

dir.create('analysis/03_Immune_Subtype/')
subtype.color=mycolors[c(5,3)]
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Times"),
          axis.text = element_text(size = 12), 
          axis.title = element_text(size = 14),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          # plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 14, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}
##################################### GSE30929
load("raw_datas/Preprocessed/gse30929_exprs.RData")
load("raw_datas/Preprocessed/gse30929_cli.RData")

exprs_data=gse30929_exprs
surv_data=gse30929_cli[,c('geo_accession','DRFS','DRFS.time')]
colnames(surv_data)=c('sample','status','time')
mode(surv_data$time)="numeric"
mode(surv_data$status)="numeric"
dim(surv_data) ## 140   3
head(surv_data)

load("raw_datas/TME/gse30929.immu.ssgsea.RData")
load("raw_datas/TME/gse30929.immu.ssgsea_28.RData")

geo.immu.ssgsea=t(rbind(gse30929.immu.ssgsea_28))
all(rownames(geo.immu.ssgsea)==surv_data$sample)

##############
immu_data=geo.immu.ssgsea
dim(immu_data)
##############
cox_data=cbind(surv_data[,c('status','time')],immu_data)
cox_data=cox_data[which(cox_data$time>0),]
dim(cox_data) ## 140*30
cox_data[1:4,1:4]

cox.res=cox_batch(t(scale(cox_data[,-c(1,2)]))
                  ,time = cox_data$time
                  ,event = cox_data$status)
dim(cox.res)
table(cox.res$p.value<0.05)
table(cox.res$p.value<0.01)
cox.res[which(cox.res$p.value<0.05),]
# p.value        HR Low 95%CI High 95%CI
# Activated CD4 T cell     0.012591232 1.3723156 1.0702409  1.7596507
# Eosinophil               0.005724637 0.6239327 0.4465030  0.8718687
# Mast cell                0.016534187 0.6829679 0.5000212  0.9328507
# Memory B cell            0.010299372 1.4162090 1.0856147  1.8474766
# T follicular helper cell 0.002085433 0.6372212 0.4782560  0.8490241


library(ConsensusClusterPlus)
clusterAlg_use=c('pam','hc','km','kmdist')[2]
distance_use=c('pearson','spearman','euclidean','canberra','maximum','minkowski')[1]
# "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans".
#########
clusterAlg_use
distance_use
#######
# df_exp=t(immu_data)[which(cox.res$p.value<0.05),]
df_exp=t(immu_data)
# df_exp=sweep(df_exp,1,apply(df_exp, 1, mean))
# df_exp=sweep(df_exp,1,apply(df_exp, 1, median))
df_exp=t(scale(t(df_exp)))
# df_exp=as.dist(1-cor(df_exp,method = 'spearman'))
# df_exp=dist(t(df_exp),method = 'canberra')
dim(df_exp)

clust_subtype = ConsensusClusterPlus(df_exp
                                     , maxK = 10, reps = 500
                                     , pItem = 0.8, pFeature = 1
                                     , title = "GSE30929_subtype"
                                     , clusterAlg = clusterAlg_use
                                     , seed = 123456
                                     , distance = distance_use
                                     , innerLinkage='complete'
                                     , finalLinkage="complete"
                                     , plot = "pdf", writeTable = T)

k1=2
gse30929.subtype=data.frame(clust_subtype[[k1]]$consensusClass)
colnames(gse30929.subtype)=c('Cluster')
gse30929.subtype$Cluster=paste0('C',gse30929.subtype$Cluster)
table(gse30929.subtype)
gse30929.subtype$Cluster[which(gse30929.subtype$Cluster=="C1")]="Immunity_H"
gse30929.subtype$Cluster[which(gse30929.subtype$Cluster=="C2")]="Immunity_L"

write.table(gse30929.subtype,'analysis/03_Immune_Subtype/gse30929.subtype.txt',sep = "\t",quote = F)

library(survcomp)
p=ggplotKMCox(data.frame(time = surv_data$time
                       , event = surv_data$status
                       , groups = gse30929.subtype$Cluster)
            , title='Groups'
            , pal = subtype.color
            , labs = c('Immunity_H','Immunity_L')
            , add_text = '')
p
savePDF('analysis/03_Immune_Subtype/gse30929.subtype.km.pdf',p,height = 6,width = 6)

ht=Heatmap(as.matrix(t(scale((immu_data))))
        , name = "ssGSEA"
        # , row_split = tme.type$Cluster
        , cluster_rows = F
        , cluster_row_slices = T
        # , row_title_gp = gpar(fill = mycolor)
        # , column_km = 3
        # , column_km_repeats=500
        , show_row_dend = F
        , column_split = gse30929.subtype[rownames(immu_data), 'Cluster']
        , cluster_columns = F
        , cluster_column_slices = T
        , show_column_dend = F
        , show_column_names = F
        , col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))
        , column_title_gp = gpar(fill = subtype.color)
        , border = TRUE)

pdf("analysis/03_Immune_Subtype/gse30929.immu.heatmap.pdf",height = 6,width = 8)
ht
dev.off()

# library(ComplexHeatmap)
# fh = function(x) fastcluster::hclust(dist(x))
# set.seed(20240601)
# ht=Heatmap(as.matrix(t(scale((immu_data))))
#         , name = "ssGSEA"
#         # , cluster_rows = fh, 
#         , cluster_columns = fh
#         # , row_split = tme.type$Cluster
#         , cluster_rows = F
#         , cluster_row_slices = T
#         # , row_title_gp = gpar(fill = mycolor)
#         , column_km = 2
#         , column_km_repeats=500
#         , show_row_dend = F
#         # , column_split = geo.subtype[rownames(immu_data), 'Cluster']
#         # , cluster_columns = F
#         , cluster_column_slices = T
#         , show_column_dend = F
#         , show_column_names = F
#         , col = circlize::colorRamp2(c(-3, 0, 3), c('#3B4992FF', 'white', '#EE0000FF'))
#         # , column_title_gp = gpar(fill = subtype.color)
#         , border = TRUE)
# vectorList=column_order(ht)
# names(vectorList)=paste0("C",names(vectorList))
# vectorList
# geo.subtype=rbind(data.frame(Cluster="C1",SampleID=vectorList[[1]]),
#                   data.frame(Cluster="C2",SampleID=vectorList[[2]])
#                   # data.frame(Cluster="C3",SampleID=vectorList[[3]])
#                   )
# head(geo.subtype)
# geo.subtype$SampleID=rownames(immu_data)[geo.subtype$SampleID]
# rownames(geo.subtype)=geo.subtype$SampleID
# geo.subtype=geo.subtype[surv_data$sample,]
# ggplotKMCox(data.frame(time = surv_data$time
#                        , event = surv_data$status
#                        , groups = geo.subtype$Cluster)
#             , title='Subtype'
#             # , pal = subtype.color
#             # , labs = c('C1','C2','C3','C4')
#             , add_text = '')
# 
# ########################
# set.seed(20240601)
# dim(immu_data)
# immu_data.scaled=as.matrix(scale(immu_data))
# km=kmeans(immu_data.scaled,centers = 3)
# km
# #### between_SS / total_SS =  77.4 % 
# library(fpc)
# plotcluster(immu_data.scaled,km$cluster)
# 
# geo.subtype_2=data.frame(km$cluster)
# colnames(geo.subtype_2)[1]="Cluster"
# geo.subtype$Cluster=paste0("C",geo.subtype_2$Cluster)
# head(geo.subtype_2)
# ggplotKMCox(data.frame(time = surv_data$time
#                        , event = surv_data$status
#                        , groups = geo.subtype_2$Cluster)
#             , title='Subtype'
#             # , pal = subtype.color
#             # , labs = c('C1','C2','C3','C4')
#             , add_text = '')
# all(rownames(geo.subtype)==rownames(geo.subtype_2))

library(scatterplot3d)
pca<-prcomp(immu_data, scale=T)
pca.result<-as.data.frame(pca$x)
pca.result$Subtype<-gse30929.subtype$Cluster


group=levels(factor(pca.result$Subtype))
bioCol=subtype.color
col= bioCol[match(pca.result$Subtype,group)]
#######
pca.result$colour= col

# pdf('PDFs/Fig2B.pdf',height = 6,width = 6)
scatterplot3d(pca.result[,1:3],color=pca.result$colour,
              pch = 16,angle=30,
              # cex.symbols = 2,
              box=T,type="p",
              main = "3D PCA Plot")
legend("bottom", legend =group,pch = 16, inset = -0.2, box.col="white",xpd = TRUE, horiz = TRUE,col=bioCol[1:length(group)])
dev.off()


library(ggbiplot)
p=ggbiplot(pca, scale=1, groups = gse30929.subtype$Cluster,
                  ellipse = T,ellipse.prob=0.5, circle = T,var.axes=F) +
  scale_color_manual(values = subtype.color) +
  # xlim(-3, 3) + ylim(-3, 3) +
  theme_classic() +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  xlab('PCA1') + ylab('PCA2')
p
savePDF('analysis/03_Immune_Subtype/gse30929.subtype.pca.pdf',p,height = 6,width = 6)


dat=apply(immu_data, 2, function(x){ifelse(x>median(x),"High","Low")})
dat=cbind(surv_data[,c('status','time')],dat)

res=c()
for(i in 3:ncol(dat)){
  dt=dat[,c(1,2,i)]
  colnames(dt)[3]="Groups"
  fit=survdiff(Surv(time, status) ~Groups,data = dt)
  print(fit$pvalue)
  res=rbind(res,data.frame(immu=colnames(dat)[i],pval=fit$pvalue,stringsAsFactors = F))
}
head(res)
immu_sig=res$immu[which(res$pval<0.05)]

p.all=list()
for(i in immu_sig){
  dt=dat[,c('status','time',i)]
  colnames(dt)[3]="Groups"
  # print(head(dt))
  p=ggplotKMCox(data.frame(time = dt$time
                           , event = dt$status
                           , groups = dt$Groups)
                , title=i
                , pal = subtype.color
                , labs = c('High','Low')
                , add_text = '')
  p.all=c(p.all,list(p))
}
length(p.all)
plot1=mg_merge_plot(p.all,nrow = 1,ncol = 4)
plot1

load("raw_datas/TME/gse30929.exp.estimate.RData")
dim(gse30929.exp.estimate)
head(gse30929.exp.estimate)
rownames(gse30929.exp.estimate)=gse30929.exp.estimate$ID
gse30929.exp.estimate=gse30929.exp.estimate[,-1]
colnames(gse30929.exp.estimate)=gsub("_estimate","",colnames(gse30929.exp.estimate))

all(rownames(gse30929.subtype)==gse30929_cli$geo_accession)
gse30929_cli_all=data.frame(gse30929_cli,Groups=gse30929.subtype$Cluster,stringsAsFactors = F)

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
    scale_fill_manual(values =subtype.color ) +
    scale_colour_manual(values = subtype.color) +
    theme_niwot()
  p
  cmp=data.frame(combn(c(1:length(unique(dt1$feature))),2))
  cmp
  p=p+ggpubr::stat_compare_means(comparisons = as.list(cmp),method = 'wilcox.test',label= "p.signif")
  p
  p.all=c(p.all,list(p))
}
length(p.all)
plot2=mg_merge_plot(p.all,nrow = 1,ncol = 3,labels = LETTERS[1])
plot2

figure=mg_merge_plot(plot1,plot2,nrow = 2,ncol = 1)
savePDF('analysis/03_Immune_Subtype/gse30929.subtype.tme.pdf',figure,height = 12,width = 18)


res=mg_limma_DEG(exp = exprs_data,group=gse30929.subtype$Cluster,ulab = 'Immunity_H',dlab = 'Immunity_L')
res$Summary
# 1.2-fold    1.3-fold   1.5-fold  2-fold  
# p<0.05   "1413|1496" "1001|800" "547|282" "158|63"
# p<0.01   "1316|1474" "959|794"  "538|281" "157|63"
# FDR<0.05 "1378|1488" "988|797"  "546|282" "158|63"
# FDR<0.01 "1253|1434" "939|783"  "535|276" "157|61"
gse30929.immu.deg.res=res$DEG
write.csv(gse30929.immu.deg.res,file = 'analysis/03_Immune_Subtype/gse30929.immu.deg.res.csv')

cutFC=log2(1.5)
gse30929.immu.deg.sig=gse30929.immu.deg.res %>% filter(abs(logFC)>cutFC & adj.P.Val<0.05)
dim(gse30929.immu.deg.sig)
# [1] 828   6
write.csv(gse30929.immu.deg.sig,'analysis/03_Immune_Subtype/gse30929.immu.deg.sig.csv')

###############################
Up.genes=rownames(gse30929.immu.deg.sig)[which(gse30929.immu.deg.sig$logFC>0)]
Down.genes=rownames(gse30929.immu.deg.sig)[which(gse30929.immu.deg.sig$logFC<0)]
length(c(Up.genes,Down.genes))

deg.up.enrich=mg_clusterProfiler(genes =Up.genes)
deg.down.enrich=mg_clusterProfiler(genes =Down.genes)

#########################################################
DB.color=c("#F97B72","#F2B701","#80BA5A","#3969AC")

df=deg.up.enrich$Enrich_tab
write.csv(df,file = 'analysis/03_Immune_Subtype/gse30929.immu.up.enrich.res.csv',row.names = F)

dt=df %>% filter(FDR<0.05) %>%
  dplyr::select(description,FDR,DB) %>%
  mutate(DB=case_when(DB=='geneontology_Biological_Process' ~ "GO_BP",
                      DB=='geneontology_Cellular_Component' ~ "GO_CC",
                      DB=='geneontology_Molecular_Function' ~ "GO_MF",
                      DB=='pathway_KEGG' ~ "KEGG")) %>%
  group_by(DB) %>%
  slice_head(n=10) %>% data.frame()

dt$FDR=-log10(dt$FDR)
dt$ID=1:nrow(dt)
head(dt)

up.enrich.barplot=ggpubr::ggbarplot(dt, x='ID', y="FDR", fill = "DB", color = "white", 
                                    orientation = "horiz",   
                                    palette = DB.color,    
                                    legend = "right",   
                                    sort.val = "asc",    
                                    sort.by.groc1s = TRUE) +    
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(labels=dt$description,expand=c(0,0))+
  labs(title = "Up-regulated",x=NULL,y="-log10(p.adjust)")


df=deg.down.enrich$Enrich_tab
write.csv(df,file = 'analysis/03_Immune_Subtype/gse30929.immu.down.enrich.res.csv',row.names = F)

dt=df %>% filter(FDR<0.05) %>%
  dplyr::select(description,FDR,DB) %>%
  mutate(DB=case_when(DB=='geneontology_Biological_Process' ~ "GO_BP",
                      DB=='geneontology_Cellular_Component' ~ "GO_CC",
                      DB=='geneontology_Molecular_Function' ~ "GO_MF",
                      DB=='pathway_KEGG' ~ "KEGG")) %>%
  group_by(DB) %>%
  slice_head(n=10) %>% data.frame()

dt$FDR=-log10(dt$FDR)
dt$ID=1:nrow(dt)
head(dt)

down.enrich.barplot=ggpubr::ggbarplot(dt, x='ID', y="FDR", fill = "DB", color = "white", 
                                      orientation = "horiz",   
                                      palette = DB.color,    
                                      legend = "right",   
                                      sort.val = "asc",    
                                      sort.by.groc1s = TRUE) +  
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(labels=stringr::str_wrap(dt$description, width = 60),expand=c(0,0))+
  labs(title = "Down-regulated",x=NULL,y="-log10(p.adjust)")
down.enrich.barplot
plot=mg_merge_plot(up.enrich.barplot,down.enrich.barplot,labels = LETTERS[3:4])
plot
savePDF('analysis/03_Immune_Subtype/gse30929.immmu.degs.enrich.barplot.pdf',plot,height = 10,width = 16)
