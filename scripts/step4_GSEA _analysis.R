

library(GseaVis)
library(clusterProfiler)
library(org.Hs.eg.db)

mg_limma_DEG=function(exp,group,ulab,dlab){
  library(limma)
  ind1=which(group==ulab)
  ind2=which(group==dlab)
  
  sml <- c(rep('G1',length(ind1)),rep('G0',length(ind2)))    # set group names
  eset=exp[,c(ind1,ind2)]
  fl <- as.factor(sml)
  
  design <- model.matrix(~fl+0)
  colnames(design) <- levels(fl)
  cont.matrix<-makeContrasts(contrasts='G1-G0',levels=design)
  #print(head(eset))
  fit<-lmFit (eset,design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  #print(sml)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(eset))
  regulated=ifelse(tT$logFC>0,'Up','Down')
  lfcs=c(log2(1.2),log2(1.3),log2(1.5),1)
  all.deg.cnt=cbind()
  for(lfc in lfcs){
    deg1=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.05)]
    deg2=regulated[which(abs(tT$logFC)>lfc&tT$P.Value<0.01)]
    deg3=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.05)]
    deg4=regulated[which(abs(tT$logFC)>lfc&tT$adj.P.Val<0.01)]
    all.deg.cnt=cbind(all.deg.cnt,c(paste0(sum(deg1=='Up'),'|',sum(deg1=='Down'))
                                    ,paste0(sum(deg2=='Up'),'|',sum(deg2=='Down'))
                                    ,paste0(sum(deg3=='Up'),'|',sum(deg3=='Down'))
                                    ,paste0(sum(deg4=='Up'),'|',sum(deg4=='Down'))))
  }
  row.names(all.deg.cnt)=c('p<0.05','p<0.01','FDR<0.05','FDR<0.01')
  colnames(all.deg.cnt)=paste0(c('1.2','1.3','1.5','2'),'-fold')
  return(list(Exp=eset,Group=group[c(ind1,ind2)],DEG=tT,Summary=all.deg.cnt))
}


load("raw_datas/Preprocessed/gse30929_exprs.RData")
gse30929.group=read.table('analysis/04_Prognostic_Model/gse30929.group.txt',sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
head(gse30929.group)

# enrichment
res=mg_limma_DEG(exp = gse30929_exprs,group=gse30929.group$Group,ulab = 'High',dlab = 'Low')
res$Summary
dim(res$DEG)
dt=res$DEG %>% 
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene,logFC) %>%
  dplyr::arrange(desc(logFC))


geneList=dt$logFC
names(geneList)=dt$Gene
head(geneList)

gseaRes <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "BP",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

# bacth plot
terms <- gseaRes@result$ID[1:8]

# plot
lapply(terms, function(x){
  gseaNb(object = gseaRes,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.6,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
library(cowplot)
figure=cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')
ggsave2('analysis/04_Prognostic_Model/gse30929.gseaRes.pdf',figure,height = 16,width = 12)
save.image(file = '20240524_Liposarcomas_GSEA.RData')


gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813')

# retain curve,
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       subPlot = 1)
# retain curve and heatmap, 
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       subPlot = 2)

# wrap the term title, termWidth
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       subPlot = 2,
       termWidth = 30)


# add gene in specific pathway
mygene <- c("TOP2A","PRC1","NUSAP1","ASPM","CDC20","CCNB1","BUB1B","CDK1")

# plot
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       subPlot = 2,
       addGene = mygene)

# new style GSEA
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T)
# new style GSEA remove point, 
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       addPoint = F)
# change heatmap color,
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       addPoint = F,
       newHtCol = c("blue","white", "red"))
# new style GSEA with gene name, 
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       addGene = mygene)

gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       rmSegment = T,
       addGene = mygene)
# remove heatmap, 
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       rmSegment = T,
       rmHt = T,
       addGene = mygene)

# add pvalue and NES, 
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       addGene = mygene,
       addPval = T)

# control label ajustment, 
gseaNb(object = gseaRes,
       geneSetID = 'GO:0098813',
       newGsea = T,
       addGene = mygene,
       addPval = T,
       pvalX = 0.75,pvalY = 0.8,
       pCol = 'black',
       pHjust = 0)

# bacth plot
terms <- gseaRes@result$ID[1:8]

# plot
lapply(terms, function(x){
  gseaNb(object = gseaRes,
         geneSetID = x,
         addPval = T,
         pvalX = 0.6,pvalY = 0.6,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
library(cowplot)
figure=cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')
ggsave2('analysis/04_Prognostic_Model/gse30929.gsea.pdf',figure,height = 16,width = 12)


# all plot
gseaNb(object = gseaRes,
       geneSetID = terms[1:3])

