rm(list = ls())
library(tidyverse)
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("GOplot")
library(patchwork)
getwd()
sig_RNA<-read.table("gene.txt",header = T)
head(sig_RNA)

gene.df <- bitr(sig_RNA$SYMBOL,
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = org.Hs.eg.db)


id_gsym<-gene.df
id_gsym$ENTREZID <- as.character(id_gsym$ENTREZID)

library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(ggplot2)
library(org.Hs.eg.db)
library(GOSemSim)
library(reshape2)
library(ggplot2)

Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
rt<-id_gsym

#用godata()函数来构建相应物种的Molecular Function本体的GO DATA
mf <- godata('org.Hs.eg.db', ont="MF", computeIC = FALSE)
#用godata()函数来构建相应物种的Cellular Component本体的GO DATA
cc <- godata('org.Hs.eg.db', ont="CC", computeIC = FALSE)
#用godata()函数来构建相应物种的Biological Process本体的GO DATA
bp <- godata('org.Hs.eg.db', ont="BP", computeIC = FALSE)

#用mgeneSim来计算MF本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simmf <- mgeneSim(id_gsym$ENTREZID, semData = mf, measure = "Wang", drop = NULL, combine = "BMA")
#用mgeneSim来计算CC本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simcc <- mgeneSim(id_gsym$ENTREZID, semData = cc, measure = "Wang", drop = NULL, combine = "BMA")
#用mgeneSim来计算BP本体，基因之间的语义相似度，结果为一个行列相同的矩阵
simbp <- mgeneSim(id_gsym$ENTREZID, semData = bp, measure = "Wang", drop = NULL, combine = "BMA")

#计算基因在MF本体和CC本体下的几何平均值，一个打分值同时包括基因的分子功能和细胞定位两个信息
fsim <- sqrt(simmf * simcc)
#或者计算基因在MF、CC、BP本体下的几何平均值
#fsim <- (simmf * simcc * simbp)^(1/3)

#将基因的名字由ENTREZID改为gene symbol，方便看懂。
colnames(fsim) = id_gsym$SYMBOL
rownames(fsim) = id_gsym$SYMBOL

#将基因自己和自己的相似度设为NA，方便接下来去掉。
for (i in 1:ncol(fsim)){
  fsim[i,i] <- NA
}

y <- melt(fsim) #把宽格式数据转化成长格式，其实就是把正方形矩阵转成三列
y <- y[!is.na(y$value),] #删掉带NA的行

# 把每两个基因之间的相似度保存到文件，只需要保存第一列基因名和第三列数值
write.csv(y[,c(1,3)], "very_easy_input.csv", row.names = F)
## 开始画图
y <- read.csv("very_easy_input.csv")
head(y)

#计算每个基因跟其他基因相似度的平均值
y.mean <- aggregate(.~Var1,y,mean) 
m <- y.mean$value
names(m) <- y.mean$Var1
#按平均值给基因名排序，便于画图
y$Var1 <- factor(y$Var1, levels=names(sort(m)))

f <- function(y) {
  r <- quantile(y, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  r[3] <- mean(y)
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

p1 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  scale_fill_brewer(palette="Set3") + #配色
  guides(fill=FALSE) + #不显示图例
  
  stat_summary(fun.data= f, geom='boxplot') + 
  geom_hline(aes(yintercept=0.75), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细 
p1

p2 <- ggplot(y, aes(Var1, value, fill = factor(Var1))) + 
  guides(fill=FALSE) + #不显示图例
  
  stat_summary(fun.data= f, geom='boxplot') + 
  #geom_hline(aes(yintercept=0.75), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细 
p2

#devtools::install_github("GuangchuangYu/gglayer")
require(gglayer)
library(gglayer)

y <- read.csv("very_easy_input.csv")
head(y)
unique(y$Var1)

#计算每个分组的平均值
y.mean <- aggregate(.~Var1,y,mean) 
m <- y.mean$value
names(m) <- y.mean$Var1
#按平均值给分组排序，便于画图
y$Var1 <- factor(y$Var1, levels=names(sort(m))) 

source("R_rainclouds.R")

p3 <- ggplot(y, aes(Var1, value, fill = Var1)) +
  #scale_fill_brewer(palette="Set2") + #配色
  guides(fill=FALSE) +
  geom_flat_violin(position=position_nudge(x=.1)) +
  
  #分散不重叠的点图
  #geom_jitter(aes(color=Var1), width=.15) + guides(color=FALSE) +
  #堆叠的点图
  geom_dotplot(binaxis="y", stackdir="down", dotsize=.35) +
  
  geom_boxplot(width=.1, position=position_nudge(x=.1)) +
  #geom_hline(aes(yintercept=0.75), linetype="dashed") + #画一条虚线
  
  coord_flip() + # x、y坐标轴互换
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(family = "Arial", size = 16, face = "bold"),
        axis.text.y = element_text(family = "Arial", size = 16, face = "bold")) + 
  theme_bw() + 
  theme(panel.border=element_rect(size=1)) #边框粗细 

p3
