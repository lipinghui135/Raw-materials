# load("raw_datas/Preprocessed/gse21122_cli.RData")
# load("raw_datas/Preprocessed/gse21122_exprs.RData")
# load("raw_datas/Preprocessed/gse159659_cli.RData")
# load("raw_datas/Preprocessed/gse159659_exprs.RData")

#################################################### 2-WGCNA
dir.create('analysis/02_WGCNA',recursive = T)

load("raw_datas/Preprocessed/gse21122.gse159659.expr_adjusted.RData")
load("raw_datas/Preprocessed/gse21122.gse159659.clinical_data.RData")
dim(expr_adjusted)
expr_adjusted[1:4,1:4]
range(expr_adjusted)
# [1]  1.498117 15.467317
head(clinical_data)

##
expr=expr_adjusted
clinicalData=clinical_data
##################
library(WGCNA)
dim(expr) ## 15548   98
datExpr = expr %>% as.matrix()
dim(datExpr) # 15548   98
range(datExpr)


geneFilter <- apply(datExpr, 1, function(x) mean(x) > 1)
length(geneFilter) ## 11771

geneMADs <- apply(datExpr, 1, mad)

topGenes <- names(sort(geneMADs, decreasing = TRUE))[1:(length(geneMADs) * 0.9)]
length(topGenes) ## 10593
##########
datExpr <- datExpr[,]
dim(datExpr)
# [1] 11771   143

datExpr = t(datExpr) %>% 
  as.data.frame() %>% 
  na.omit()
dim(datExpr) ## 143 11771
range(datExpr)
# [1]  1.498117 15.467317

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
# [1] TRUE

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower = sft$powerEstimate
softPower ## 6

######################################################################

pdf("analysis/02_WGCNA/wgcna.softPower.pdf",height = 5,width = 10)
par(mfrow = c(1, 2))
cex1 <- 0.85

（scale-free topology fit index）
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red") 
mean connectivity）
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

######################################################################
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = 60
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 1, minClusterSize = minModuleSize,
                            method = "tree")
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
table(table(dynamicColors)>30)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# （Module Eigengenes）
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes

MEDiss <- 1 - cor(MEs)


METree <- hclust(as.dist(MEDiss), method = "average")


mergeCutHeight <- 0.25  
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = mergeCutHeight, verbose = 3)

mergedColors <- merge$colors
mergedMEs <- merge$newMEs
##########################################################################
pdf("analysis/02_WGCNA/wgcna.geneTree.pdf",height = 3,width = 5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()



moduleEigengenes = moduleEigengenes(datExpr, colors = mergedColors)$eigengenes
MEs = orderMEs(moduleEigengenes)
# all(MEs==merge$newMEs)

#
clinicalData$sample_type=factor(clinicalData$sample_type,levels = c("Normal control","Liposarcoma"))
clinicalData$sample_type=as.numeric(clinicalData$sample_type)

moduleTraitCor = cor(MEs, clinicalData$sample_type, use = "pairwise.complete.obs")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))



textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
dev.off()
pdf("analysis/02_WGCNA/labeledHeatmap.pdf",height = 5,width = 5)
par(mar = c(1, 10, 1, 1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinicalData)[2],
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               xLabelsAngle=0,
               xLabelsAdj =0.5,
               main = paste("Module-trait relationships"))
dev.off()

################
Modules=data.frame(Gene=colnames(datExpr),Module=mergedColors,stringsAsFactors = F)
rownames(Modules)=Modules$Gene
head(Modules)
write.csv(Modules,file = 'analysis/02_WGCNA/geo.wgcna.modules.genes.csv')

table(Modules$Module)


## Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

all(rownames(MEs)==rownames(datExpr))
geneModuleMembership <- as.data.frame(signedKME(datExpr, data.frame(MEs), outputColumnName = ""))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))


all(rownames(datExpr)==rownames(clinicalData))
geneTraitSignificance <- as.data.frame(cor(datExpr, clinicalData, use = 'pairwise.complete.obs'))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

####################################
selected_modules=c("turquoise")

table(Modules$Module)[selected_modules]
# turquoise 
# 334

turquoise_gene=Modules$Gene[which(Modules$Module %in% selected_modules)]
head(turquoise_gene)
length(turquoise_gene)
write.csv(turquoise_gene,file = 'analysis/02_WGCNA/turquoise.genes.csv',row.names = F)
################################## hub genes
geneModuleMembership[1:4,1:5]
geneTraitSignificance[1:4,]

module="turquoise"
pdf('analysis/02_WGCNA/turquoise.verboseScatterplot.pdf',height = 5,width = 5)
verboseScatterplot(abs(geneModuleMembership[turquoise_gene, 'turquoise']),
                   abs(geneTraitSignificance[turquoise_gene, 'sample_type']),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for ",'sample_type'),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
# abline(v = 0.4, col = "red", lwd = 2, lty = 2)
# abline(h = 0.4, col = "red", lwd = 2, lty = 2)
dev.off()

inds1=abs(geneModuleMembership[turquoise_gene, 'turquoise'])>0.4 & MMPvalue[turquoise_gene, 'turquoise']<0.05
inds2=abs(geneTraitSignificance[turquoise_gene, 'sample_type'])>0.4 & GSPvalue[turquoise_gene, 'sample_type']<0.05
all(rownames(geneModuleMembership)==rownames(geneTraitSignificance))
length(inds1)
length(inds2)

hub.genes=turquoise_gene[inds1 & inds2]
length(hub.genes)
# [1] 93

##################################
library(clusterProfiler)
# gene2entrez=bitr(turquoise_gene, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
# head(gene2entrez)


turquoise_bp = enrichGO(turquoise_gene,
                    OrgDb = "org.Hs.eg.db",
                    keyType = "SYMBOL",
                    ont = "BP",
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.2)
results=turquoise_bp@result
results=results[which(results$qvalue<0.05),]
dim(results)
# [1] 433   9
write.csv(results,file = 'analysis/02_WGCNA/turquoise_GO_BP_all.csv')

dt=results %>% arrange(p.adjust) %>% top_n(n=10)
dt$p.adjust=format(dt$p.adjust, digits =3,scientific =TRUE)
mode(dt$p.adjust)="numeric"

p=ggplot(dt, aes(Count,reorder(Description,Count), fill = p.adjust)) +
  geom_col() +
  theme_bw() +
  scale_fill_gradient(low = "darkorange", high = "lightcyan") +
  labs(title = "Biological process of MEturquoise", y = "")
p
savePDF('analysis/02_WGCNA/wgcna.module.bp.barplot.pdf',p,height = 5,width = 8)

##############
geo.deg.sig=read.csv("analysis/01_DEGs/geo.deg.sig.csv",row.names = 1)
geo.deg.genes=rownames(geo.deg.sig)
head(geo.deg.genes)
length(geo.deg.genes) ## 852

turquoise_gene=read.csv(file = 'analysis/02_WGCNA/turquoise.genes.csv',stringsAsFactors = F)
turquoise_gene=as.character(turquoise_gene$x)

intersect(turquoise_gene,geo.deg.genes)
length(intersect(turquoise_gene,geo.deg.genes)) ## 101

myList=list(WGCNA=turquoise_gene,DEGs=geo.deg.genes)

mycolors=c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#F97B72","#4B4B8F","#A5AA99")

p=plot(eulerr::venn(myList),labels = list(col = "gray20", font = 2), 
              edges = list(col = "gray60", lex = 1),
              fills = list(fill = c("#11A579","#3969AC","#F2B701","#E73F74"), alpha = 1),
              quantities = list(cex = .8, col = 'gray20'))
p
savePDF('analysis/02_WGCNA/wgcna.degs.venn.pdf',p,height = 5,width = 6)

wgcna.degs.intersect.genes=Reduce(intersect,list(WGCNA=turquoise_gene,DEGs=geo.deg.genes))
length(wgcna.degs.intersect.genes)

write.csv(wgcna.degs.intersect.genes,file = 'analysis/02_WGCNA/wgcna.degs.intersect.genes.csv',row.names = F)

