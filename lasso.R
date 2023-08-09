library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)
library(patchwork)
library(limma)
rm(list = ls())
getwd()
setwd("E:/麻醉/WGCNA_Pro_top5000")
GSE="Pro"
#rt=read.table(paste0(GSE,".txt"),sep="\t",header=T,check.names=F)
load("E:/麻醉/03_差异分析/差异分析_input_Pro.Rdata")
load("E:/麻醉/WGCNA_Pro_top5000/moduleGenes.Rdata")
exp<-exp[moduleGenes,]
rt<-exp
rt$geneNames<-rownames(rt)
rt<-rt[,c(31,1:30)]
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)
data<-rt
data=t(data)
#data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample_Pro.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)


#控制组放置最前
afcon=as.matrix(table(sample[,1]))[1,1]
afcon=as.vector(afcon)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])

set.seed(666666)
cvfit = cv.glmnet(x, y,family = "binomial", nlambda=100, alpha=1,nfolds = 10) #这里alpha=1为LASSO回归，如果等于0就是岭回归，10乘交叉验证
#参数 family 规定了回归模型的类型：
#family="gaussian" 适用于一维连续因变量（univariate）
#family="mgaussian" 适用于多维连续因变量（multivariate）
#family="poisson" 适用于非负次数因变量（count）
#family="binomial" 适用于二元离散因变量（binary）
#family="multinomial" 适用于多元离散因变量（category）
#我们这里结局指标是2分类变量，所以使用binomial

fit <- glmnet(x,y,family = "binomial")
cvfit$lambda.min

#提取信息及预测风险
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.table(geneCoef, file="geneCoef.xls", sep="\t", quote=F, row.names=F)
write.table(file="lassoset.txt",lassoGene,sep="\t",quote=F,col.names=F,row.names=F) #文件名

#######################简单作图########################################
pdf("lasso.pdf",height = 5,width = 7)
layout(matrix(c(1,1,2,2), 2, 2, byrow = F))   #两行两列，图一占前俩格，图二占后两格，按列排
#pdf("lambda.pdf")
plot(fit,xvar = 'lambda')
#dev.off()
#pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
#dev.off()
dev.off()
