rm(list=ls())
Sys.setenv(LANGUAGE = "en")

load("差异分析_input.Rdata")
getwd()
dir.create("WGCNA_Sev_top5000")
setwd("E:/麻醉/WGCNA_Sev_top5000")

fpkm_sd = apply(exp,1,sd)#1是对每一行，2是对每一列
## 使用标准差对基因进行降序排序
fpkm_sd_sorted = order(fpkm_sd, decreasing = T)
## 选择前5000个标准差较大的基因
fpkm_num = fpkm_sd_sorted[1:5000]
## 从表达矩阵提取基因
fpkm_filter = exp[fpkm_num,]
## 对表达矩阵进行转置
datExpr<-as.data.frame(t(fpkm_filter))
datTraits<-cli
library(PerformanceAnalytics)
library(WGCNA)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)

# 设定软阈值范围
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# 获得各个阈值下的R方和平均连接度，RsquaredCut为期望的R^2阈值
sft = pickSoftThreshold(datExpr,powerVector = powers,RsquaredCut = 0.9,verbose = 5)
# 作图：
pdf(file = "Fig1A.pdf",width = 9,height = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
softPower<-sft$powerEstimate
softPower

cor <- WGCNA::cor
net = blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "unsigned", 
  minModuleSize = 100,#设置每个基因模块最少的基因数目为10。
  reassignThreshold = 0, mergeCutHeight = 0.2,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  verbose = 3
)
table(net$colors)

cor<-stats::cor
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
pdf(file = "dendrogram.pdf",width = 6,height = 5)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
# 计算基因之间的相异度
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
minModuleSize = 100
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
mergedColors = labels2colors(net$colors)
table(mergedColors)
pdf("Fig1B.pdf",width = 5,height = 4)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

plotTOM = dissTOM^softPower
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
#pdf("Fig1C_network_heatmap.pdf",width = 6,height = 6)
# 这一部分比较耗时，行和列同时做层级聚类
TOMplot(plotTOM, geneTree, mergedColors, col=colorRampPalette(colors = c("red","yellow","#FFFACD"))(50),
        main = "Network heatmap plot, all genes")
dev.off()
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs = orderMEs(MEs_col)
### 模块与表型数据关联
design=model.matrix(~0+ datTraits$group)
design=as.matrix(design)
colnames(design)=c("Sev","Control")
moduleColors <- labels2colors(net$colors)
nSamples = nrow(datExpr)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
#pdf("Fig2A_Module-trait_relationships.pdf",width = 6,height = 8)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               xLabelsAngle = 0,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               #showCols = 1
)
dev.off()
table(moduleColors)

### 计算模块与基因的相关性矩阵
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(datExpr, MEs, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

#性状与每个基因表达量相关性在各个模块的均值作为该性状在该模块的显著性，显著性最大模块与该性状最相关
GS1 <- as.numeric(WGCNA::cor(design[,1],datExpr,use="p",method="pearson"))
# 显著性绝对值：
GeneSignificance <- abs(GS1)
# 获得该性状在每个模块中的显著性：
ModuleSignificance <- tapply(GeneSignificance,moduleColors,mean,na.rm=T)
ModuleSignificance <- as.matrix(t(ModuleSignificance))
sd <- tapply(GeneSignificance,moduleColors,sd,na.rm=T)
SE.mean <- sd/sqrt(table(moduleColors))
#pdf("Fig2B_module_signif.pdf",width = 8,height = 6)
barplot_SS <- barplot(ModuleSignificance[1,],names.arg=F,
                      ylab = "Gene Significance",ylim = c(0,1),
                      col = colnames(ModuleSignificance),)
arrows(barplot_SS,ModuleSignificance[1,]+SE.mean,barplot_SS,ModuleSignificance[1,]-SE.mean,length =0.05, angle = 90, code = 3)
dev.off()

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(datExpr, design, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(datExpr, design, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "yellow"
pheno = "Sev"
modNames = substring(colnames(MEs), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module
#pdf("Fig2C_.pdf",width = 6,height = 6)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

##保存模块基因
moduleGenes=colnames(datExpr)[moduleGenes==TRUE]
save(moduleGenes,file='moduleGenes.Rdata')
write.csv(moduleGenes,"moduleGenes.csv")