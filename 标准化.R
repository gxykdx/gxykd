### 环境设置
rm(list = ls())
options(stringsAsFactors = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") 
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

getwd()
dir.create("02_normalized")
library(sva)
library(limma)

getwd()
GSE_ID <- c('GSE4386')
GSE_file<-paste0("01_完成注释/",GSE_ID,"_raw_Assay.csv")
rt<-read.csv(file = GSE_file,row.names = 1,check.names = F)

### 判断是否需要log
expset<-rt
qx <- as.numeric(quantile(expset, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { expset[which(expset <= 0)] <- NaN
expset<- log2(expset)
print("log2 transform finished")}else{print("log2 transform not needed")}



### 判断是否需要normalization
par(cex = 0.7)
n.sample=ncol(expset)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
p1<-boxplot(expset, col = cols,main="expression value",las=2)

getwd()

expset[1:5,1:5]
expset<- normalizeBetweenArrays(as.matrix(expset,method="scale"))
p2<-boxplot(expset,col = cols,main="expression value",las=2)
write.csv(expset,file = paste0("02_normalized/",GSE_ID,"_normalized.csv"),row.names = T,quote = F)
save(expset,file = 'normalized_exprset4386.Rdata')
