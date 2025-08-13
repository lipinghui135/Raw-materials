setwd("/home/pub252/users/Work1/Work/20240524_Liposarcomas/")
source('/home/pub252/projects/codes/mg_base.R')

#######################
load("raw_datas/Preprocessed/gse21122_cli.RData")
load("raw_datas/Preprocessed/gse21122_exprs.RData")
load("raw_datas/Preprocessed/gse159659_cli.RData")
load("raw_datas/Preprocessed/gse159659_exprs.RData")

#################################################### 1-
dir.create('analysis/01_DEGs',recursive = T)
####################################### GSE21122
cutFC=log2(1.5)
###########
res=mg_limma_DEG(exp = gse21122_exprs,group=gse21122_cli$sample_type,ulab = 'Liposarcoma',dlab = 'Normal control')
res$Summary
# 1.2-fold    1.3-fold    1.5-fold  2-fold   
# p<0.05   "2097|1612" "1515|1256" "802|857" "248|431"
# p<0.01   "1615|1355" "1280|1113" "710|794" "229|424"
# FDR<0.05 "1747|1416" "1341|1153" "733|813" "236|427"
# FDR<0.01 "1256|1157" "1065|966"  "629|730" "208|406"
gse21122.deg.res=res$DEG
head(gse21122.deg.res)
# writeMatrix(gse21122.deg.res,'analysis/01_DEGs/gse21122.deg.res.txt')

library(dplyr)
gse21122.deg.sig=gse21122.deg.res %>% filter(abs(logFC)>cutFC & adj.P.Val<0.05)
head(gse21122.deg.sig)

####################################### GSE159659
res=mg_limma_DEG(exp = gse159659_exprs,group=gse159659_cli$sample_type,ulab = 'Liposarcoma',dlab = 'Normal control')
res$Summary
# 1.2-fold    1.3-fold    1.5-fold  2-fold  
# p<0.05   "1834|1881" "1084|1402" "447|774" "79|254"
# p<0.01   "1094|1223" "767|1079"  "349|702" "75|250"
# FDR<0.05 "841|1029"  "635|938"   "310|655" "73|243"
# FDR<0.01 "259|459"   "238|454"   "151|396" "51|194"
gse159659.deg.res=res$DEG
head(gse159659.deg.res)
# writeMatrix(gse159659.deg.res,'analysis/01_DEGs/gse159659.deg.res.txt')

library(dplyr)
gse159659.deg.sig=gse159659.deg.res %>% filter(abs(logFC)>cutFC & adj.P.Val<0.05)
head(gse159659.deg.sig)


p1=mg_volcano(logfc = gse21122.deg.res$logFC,pvalue = gse21122.deg.res$adj.P.Val,
              cutFC = cutFC,cutPvalue = 0.05)+labs(title = 'GSE21122')
p2=mg_volcano(logfc = gse159659.deg.res$logFC,pvalue = gse159659.deg.res$adj.P.Val,
              cutFC = cutFC,cutPvalue = 0.05)+labs(title = 'GSE159659')
figure1A=mg_merge_plot(p1,p2,nrow = 1,ncol = 2,labels = LETTERS[1:2])
figure1A
# savePDF('analysis/01_DEGs/Figure1A.pdf',figure1A,height = 5,width = 10)

gse21122.up.genes=rownames(gse21122.deg.sig)[which(gse21122.deg.sig$logFC>cutFC)]
gse21122.down.genes=rownames(gse21122.deg.sig)[which(gse21122.deg.sig$logFC<(-cutFC))]

gse159659.up.genes=rownames(gse159659.deg.sig)[which(gse159659.deg.sig$logFC>cutFC)]
gse159659.down.genes=rownames(gse159659.deg.sig)[which(gse159659.deg.sig$logFC<(-cutFC))]

mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

library(ggsci)
myList=list(A=gse21122.up.genes,B=gse159659.up.genes,C=gse21122.down.genes,D=gse159659.down.genes)
names(myList)=c("Up-regulated of GSE21122","Up-regulated of GSE159659","Down-regulated of GSE21122","Down-regulated of GSE159659")
dev.off()
p=plot(eulerr::venn(myList),labels = list(col = "gray20", font = 2), 
       edges = list(col = "gray60", lex = 1),
       fills = list(fill = c("#11A579","#3969AC","#E68310","#F97B72"), alpha = 1),
       quantities = list(cex = .8, col = 'gray20'))
p

# savePDF('analysis/01_DEGs/Figure1B.pdf',p,height = 5,width = 5)

################
Up.genes=intersect(gse21122.up.genes,gse159659.up.genes)
length(Up.genes)
# [1] 116

Down.genes=intersect(gse21122.down.genes,gse159659.down.genes)
length(Down.genes)
# [1] 324

dt=rbind(data.frame(Gene=Up.genes,type="Up-regulated"),
      data.frame(Gene=Down.genes,type="Down-regulated"))
head(dt)
# write.csv(dt,file = 'analysis/01_DEGs/geo.degs.limma.csv',row.names = F)

deg.up.enrich=mg_clusterProfiler(genes =Up.genes)
deg.down.enrich=mg_clusterProfiler(genes =Down.genes)

#########################################################
DB.color=c("#F97B72","#F2B701","#80BA5A","#3969AC")

df=deg.up.enrich$Enrich_tab
# write.csv(df,file = 'analysis/01_DEGs/geo.up.enrich.res.csv',row.names = F)

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
# write.csv(df,file = 'analysis/01_DEGs/geo.down.enrich.res.csv',row.names = F)

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
figure1C=mg_merge_plot(up.enrich.barplot,down.enrich.barplot,labels = LETTERS[2:3])
figure1C
# savePDF('analysis/01_DEGs/Figure1C.pdf',figure1C,height = 10,width = 15)


load("raw_datas/Preprocessed/gse21122.gse159659.expr_adjusted.RData")
load("raw_datas/Preprocessed/gse21122.gse159659.clinical_data.RData")
dim(expr_adjusted)
# [1] 11771   143
expr_adjusted[1:4,1:4]
range(expr_adjusted)
# [1]  1.498117 15.467317
head(clinical_data)


expr=expr_adjusted
clinicalData=clinical_data
#############
cutFC=log2(1.5)

res=mg_limma_DEG(exp = expr,group=clinicalData$sample_type,ulab = 'Liposarcoma',dlab = 'Normal control')
res$Summary
# 1.2-fold    1.3-fold   1.5-fold  2-fold  
# p<0.05   "1469|1455" "805|1010" "300|562" "68|234"
# p<0.01   "1258|1310" "724|957"  "282|557" "68|234"
# FDR<0.05 "1340|1358" "758|973"  "292|560" "68|234"
# FDR<0.01 "1026|1154" "637|893"  "260|541" "68|231"
geo.deg.res=res$DEG
head(geo.deg.res)
writeMatrix(geo.deg.res,'analysis/01_DEGs/geo.deg.res.txt')
##########################
library(dplyr)
cutFC=log2(1.5)
geo.deg.res=read.table('analysis/01_DEGs/geo.deg.res.txt',sep = "\t",row.names = 1,header = T)
head(geo.deg.res)
geo.deg.sig=geo.deg.res %>% filter(abs(logFC)>cutFC & adj.P.Val<0.05)
dim(geo.deg.sig)
# [1] 852   6
write.csv(geo.deg.sig,'analysis/01_DEGs/geo.deg.sig.csv')


################
geo.up.genes=rownames(geo.deg.sig)[which(geo.deg.sig$logFC>0)]
geo.down.genes=rownames(geo.deg.sig)[which(geo.deg.sig$logFC<0)]
length(geo.up.genes)
# [1] 292
length(geo.down.genes)
# [1] 560

Up.genes=geo.up.genes
length(Up.genes)
Down.genes=geo.down.genes
length(Down.genes)

dt=rbind(data.frame(Gene=Up.genes,type="Up-regulated"),
         data.frame(Gene=Down.genes,type="Down-regulated"))
head(dt)
# write.csv(dt,file = 'analysis/01_DEGs/geo.degs.limma.csv',row.names = F)

########################## 
p1=mg_volcano(logfc = geo.deg.res$logFC,pvalue = geo.deg.res$adj.P.Val,
              cutFC = cutFC,cutPvalue = 0.05,colors=c("#223D6C","grey","#E0367A"))+labs(title = '')

p1
savePDF('analysis/01_DEGs/geo.degs.volcano.pdf',p1,height = 6,width = 6)
##########################
up_50=geo.deg.sig %>% arrange(desc(logFC)) %>% dplyr::top_n(50,logFC) %>% rownames()
down_50=geo.deg.sig %>% arrange((logFC)) %>% dplyr::top_n(-50,logFC) %>% rownames()

library(ComplexHeatmap)
dt=expr[c(up_50,down_50),rownames(clinicalData)]
dim(dt)

pdf('analysis/01_DEGs/geo.degs.heatmap.pdf',height = 6,width = 6)
Heatmap(as.matrix(t(scale(t(dt))))
        , name = "Expression"
        , row_split = c(rep('Up',length(up_50)),rep("Down",length(down_50)))
        , cluster_rows = T
        , cluster_row_slices = T
        , row_title_gp = gpar(fill = c("#E0367A","#223D6C"))
        , show_row_dend = F
        , column_split = clinicalData$sample_type
        , column_title_gp = gpar(fill = mycolors[c(2,3)])
        , cluster_columns = T
        , cluster_column_slices = T
        , show_column_dend = F
        , show_column_names = F
        , show_row_names = F
        , row_names_gp = gpar(fontsize = 10)
        , col = circlize::colorRamp2(c(-4, 0, 4), c('#3B4992FF', 'white', '#EE0000FF'))
        , border = TRUE)
dev.off()


deg.up.enrich=mg_clusterProfiler(genes =Up.genes)
deg.down.enrich=mg_clusterProfiler(genes =Down.genes)

#########################################################
DB.color=c("#F97B72","#F2B701","#80BA5A","#3969AC")

df=deg.up.enrich$Enrich_tab
write.csv(df,file = 'analysis/01_DEGs/geo.up.enrich.res.csv',row.names = F)

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
write.csv(df,file = 'analysis/01_DEGs/geo.down.enrich.res.csv',row.names = F)

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
savePDF('analysis/01_DEGs/geo.degs.enrich.barplot.pdf',plot,height = 10,width = 16)

