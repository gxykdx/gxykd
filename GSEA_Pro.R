rm(list = ls())
library(patchwork)
library(org.Rn.eg.db)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
set.seed(666666)
load("E:/麻醉/03_差异分析/deg_Pro.Rdata")
deg2$SYMBOL<-rownames(deg2)
df<-deg2[,c(1,7)]
head(df)#查看前面几行
dim(df)#数据总共几行几列
df_id<-bitr(df$SYMBOL, #转换的列是df数据框中的SYMBOL列
            fromType = "SYMBOL",#需要转换ID类型
            toType = "ENTREZID",#转换成的ID类型
            OrgDb = "org.Hs.eg.db")#对应的物种，小鼠的是org.Mm.eg.db
df_all<-merge(df,df_id,by="SYMBOL",all=F)
head(df_all)
dim(df_all)
df_all_sort<-df_all[order(df_all$logFC, decreasing = T),]#降序排序
gene.fc = df_all_sort$logFC#把foldchange按照从大到小提取出来
head(gene.fc)
names(gene.fc) <- df_all_sort$ENTREZID#给上面提取的foldchange对应上ENTREZID
head(gene.fc)
df_all_sort <- df_all[order(df_all$logFC, decreasing = T),]#先按照logFC降序排序
gene_fc = df_all_sort$logFC #把foldchange按照从大到小提取出来
head(gene_fc)
names(gene_fc) <- df_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(gene_fc)
KEGG <- gseKEGG(gene_fc, organism = "hsa")
BP <- gseGO(
  gene_fc, #gene_fc
  ont = "BP",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)
head(KEGG)
KEGG1<- KEGG[order(KEGG$enrichmentScore, decreasing = T),]
paths <- c("hsa04657","hsa04668","hsa05144","hsa04061","hsa04064","hsa04917","hsa04625","hsa04060","hsa04659")
p1<-gseaplot2(KEGG,paths, pvalue_table = F)

BP1<- BP[order(BP$enrichmentScore, decreasing = T),]
paths <- c("GO:0051918","GO:0051917","GO:0051412","GO:0010273","GO:1990169","GO:0048245","GO:0042730","GO:0051447","GO:0061687")
p2<-gseaplot2(BP,paths, pvalue_table = F)

CC <- gseGO(
  gene_fc, #gene_fc
  ont = "CC",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

CC1<- CC[order(CC$enrichmentScore, decreasing = T),]
paths <- c("GO:0005667","GO:0005882","GO:0072562","GO:0019814","GO:0042571")
p3<-gseaplot2(CC,paths, pvalue_table = F)

MF <- gseGO(
  gene_fc, #gene_fc
  ont = "MF",# "BP"、"MF"和"CC"或"ALL"
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
)

MF1<- MF[order(MF$enrichmentScore, decreasing = T),]
paths <- c("GO:0017017","GO:0035259","GO:0033549","GO:0048019","GO:0045236","GO:0008009","GO:0035014","GO:0042379","GO:0048020")
p4<-gseaplot2(MF,paths, pvalue_table = F)
p2+p3+p4+p1

write.csv(BP1,"BP_Pro.csv")
write.csv(CC1,"CC_Pro.csv")
write.csv(MF1,"MF_Pro.csv")
write.csv(KEGG1,"KEGG_Pro.csv")
