rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


library("BiocManager")
library("devtools") 
library("dplyr")
library("tidyr")
library("GEOquery")
library("Biobase")
library(devtools)
library(GEOmirror)
library(idmap1)
library(idmap2)
library(idmap3)
library(dplyr)
library(tidyr)
dir.create("01_下载GEO数据/")
setwd("01_下载GEO数据/")

### 想要下载的数据集
GSE_ID <- c('GSE4386')
gset<-lapply(GSE_ID,function(GSE_ID){
  geoChina(gse=GSE_ID)
})

getwd()
setwd("../")

### Load刚下载的数据
GSE_ID <- c('GSE4386')
GSE_file<-paste0("01_下载GEO数据/",GSE_ID,"_eSet.Rdata") 

load(GSE_file)
class(gset)
length(gset)

exprSet <- gset[[1]]
str(exprSet, max.level = 2)

### 未经注释的芯片表达矩阵【assayData】
assayData <- exprs(exprSet)
dim(assayData)
assayData[1:5, 1:6]

### 临床表型【phenoData】
phenoData <- pData(exprSet)
dim(phenoData)
head(phenoData[,1:5])

### 想要的列数的临床信息【meta】
col<-c("title","source_name_ch1")
meta<-phenoData[, col]
table(meta[,2])
### 平台GPL信息
gpl <- exprSet@annotation

### 法1: 使用idamp1/2/3下载soft文件注释 【推荐】
# featureData =get_soft_IDs(gpl) # 使用idamp1
featureData=getIDs(gpl) # 使用idamp2
# featureData=get_pipe_IDs(gpl) # 使用idamp3


head(featureData)[,1:3]
head(assayData)[,1:3]
colnames(featureData)
featureData<-featureData[,c("probe_id","symbol")]
colnames(featureData)<-c("ID","symbol")
featureData <- featureData[featureData$symbol != '', ]##gene

index<-intersect(rownames(assayData),featureData$ID)
assayData<-assayData[index,]
rownames(featureData)<-featureData$ID
featureData<-featureData[index,]
identical(rownames(assayData),featureData$ID)

class(assayData)
newAssayData<-as.data.frame(assayData)
identical(rownames(featureData),rownames(newAssayData))
featureData$max <- apply(newAssayData, 1, max) 

featureData[1:15,1:3]
featureData <- featureData[order(featureData$symbol, ##gene
                                 featureData$max,
                                 decreasing = T), ]

dim( featureData )
featureData <- featureData[!duplicated(featureData$symbol), ]##表达矩阵仅保留最大表达量的探针
dim( featureData )


### 到这里表达矩阵的处理就结束了，代码比较繁杂
### 也可以选择，直接舍弃“///”列，或者使用R包去注释数据
newAssayData<-newAssayData[featureData$ID,]
identical(rownames(newAssayData),rownames(featureData))

rownames(newAssayData) <-featureData$symbol
newAssayData[1:5, 1:6]

dir.create("01_完成注释")
write.csv(newAssayData,file = paste0("01_完成注释/",GSE_ID,"_raw_Assay.csv"),row.names = T,quote = F)
save(phenoData,file = paste0("01_完成注释/",GSE_ID,"_cli.Rdata"))
