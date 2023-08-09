rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
library(Biobase)
library(limma)

GSE_ID <- c('GSE4386')
cliRdata = paste0("01_完成注释/",GSE_ID,"_cli.Rdata")
expData = paste0("02_normalized/",GSE_ID,"_normalized.csv")


### 准备分组信息
exp<-read.csv(file =expData ,row.names = 1,check.names = F)
load(cliRdata)
cli_Sev<-read.csv("cli_Sev.CSV",row.names = 1,header = T)
cli<-cli_Sev
table(cli$group)
cli_back<-cli
table(cli$group)
cli$group<-ifelse(cli$group=="Control","Control","Sev")
cli<-cli[order(cli$group),]
head(cli)[,1:2]
exp<-exp[,rownames(cli)]
dir.create("03_差异分析")
identical(rownames(cli),colnames(exp))
save(cli,exp,file="03_差异分析/差异分析_input.Rdata")

### 环境设置
rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
load("03_差异分析/差异分析_input.Rdata")

do_limma_array <- function(exprSet,group_list){
  suppressMessages(library(limma))
  design <- model.matrix(~0+factor(group_list))
  colnames(design)=levels(factor(group_list))
  rownames(design)=colnames(exprSet)
  design
  fit <- lmFit(exprSet, design)
  group_list
  cont.matrix=makeContrasts(contrasts=c('me-other'),levels = design)
  fit2=contrasts.fit(fit,cont.matrix)
  fit2=eBayes(fit2)
  tempOutput = topTable(fit2, coef='me-other', n=Inf)
  DEG_limma = na.omit(tempOutput)
  head(DEG_limma) 
  return(DEG_limma)
}

group_list=ifelse(cli$group=="Control",'other','me')

deg1=do_limma_array(exp,group_list)
head(deg1)[,1:5]

dir.create("04_火山图")
save(deg1,file = "04_火山图/volcano.Rdata")
save(deg1,cli,file = "03_差异分析/deg.Rdata")
