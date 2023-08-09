rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
dir.create("08_CIBERSORT")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

library(e1071)
library(parallel)
library(pacman)

getwd()
dir.create("")
source("CIBERSORT.R")

sig_matrix <- "LM22.txt"

## 指定 表达矩阵
df<-read.csv(Sys.glob("E:/麻醉/01_完成注释/*.csv"),check.names = F,row.names = 1)

df[is.na(df)] <- 0
df$gene<-rownames(df)
df<-df[,c(41,1:40)]

write.table(df,"E:/麻醉/08_CIBERSORT/exprMatt.txt",sep = "\t",row.names = F)

mixture_file = 'E:/麻醉/08_CIBERSORT/exprMatt.txt'

my_cibersort <- CIBERSORT(sig_matrix, ## LM22基因集 对应22种immnue cells
                          mixture_file, ## 表达矩阵(是否log2不重要;无需normalized;若RNA-seq:推荐TPM)
                          perm=1000, ##置换检验1000次
                          QN=TRUE) ##芯片=T;RNA-seq=F; QN:Quntile Normalization

my_cibersort<-my_cibersort[,1:22]
colnames(my_cibersort)<-gsub("\\ ","\\_",colnames(my_cibersort))## 将空格替换为_;为了下游分析时更为方便

save(my_cibersort,file = "E:/麻醉/08_CIBERSORT/CIBEROSTR.Rdata")
write.csv(my_cibersort,file = "E:/麻醉/08_CIBERSORT/CIBERSORT.csv",row.names = T,quote = F)
