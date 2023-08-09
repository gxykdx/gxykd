rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
#dir.create("08_CIBERSORT")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


library(dplyr)
library(tidyr)
library(tidyverse)

load("E:/麻醉/08_CIBERSORT/CIBEROSTR.Rdata")
cibersort_raw <- as.data.frame(my_cibersort)
colnames(cibersort_raw)<-gsub("\\_","\\ ",colnames(cibersort_raw))## 将空格替换为_;为了下游分析时更为方便
# 通过管道符一步步先将CIBERSORT_Results读入R语言中，并将其第一列列名“Mixture”修改为“Sample”。
#并赋值给cibersort_raw。
load("E:/麻醉/03_差异分析/差异分析_input.Rdata")
load("E:/麻醉/03_差异分析/差异分析_input_Pro.Rdata")
cli$X<-rownames(cli)
cli<-cli_Pro
cli$X<-rownames(cli)
cibersort_raw<-cibersort_raw[rownames(cli),]
cibersort_raw$X<-rownames(cibersort_raw)
cibersort_raw<-merge(cibersort_raw,cli,by="X")
cibersort_raw$Group<-cibersort_raw$group

### 调整数据
library(dplyr)
library(tidyr)
dd1 <- cibersort_raw %>% 
  pivot_longer(cols=2:23,
               names_to= "celltype",
               values_to = "NES")

library(ggplot2)
library(ggpubr)
### 箱线图
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = Group),outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=Group), label = "p.signif")

### 小提琴
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_violin(aes(fill = Group),position = position_dodge(1),scale = "width")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=Group), label = "p.signif")

### 混合叠加
ggplot(data =dd1, aes(x = celltype, y = NES))+
  geom_boxplot(aes(fill = Group),position = position_dodge(1),width=.3,outlier.shape = NA)+
  geom_violin(aes(colour = Group),position = position_dodge(1),scale = "width",fill=NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1, colour = "black"))+
  stat_compare_means(aes(group=Group), label = "p.signif")
