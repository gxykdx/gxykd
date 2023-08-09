rm(list = ls())
library(reshape2)
library(ggpubr)
outFile="figure_top5000pro.pdf"
load("E:/麻醉/03_差异分析/差异分析_input.Rdata")
load("E:/麻醉/03_差异分析/差异分析_input_Pro.Rdata")
anoikis<-read.table("hubgenetop5000.txt")$V1
rt<-exp[anoikis,]
rt<-as.data.frame(t(rt))
rt$sample<-rownames(rt)
cli<-cli_Pro
cli$sample<-rownames(cli)
rt_1<-merge(rt,cli,by="sample")
rownames(rt_1)<-rt_1$sample
data<-rt_1[,c(6,2:4)]
rt<-data
colnames(rt)[1]="Type"
#把数据转换成ggplot2数据文件
x=colnames(rt)[1]
data=melt(rt,id.vars=c("Type"))
colnames(data)=c("Type","Gene","Expression")
#绘制小提琴图
p=ggviolin(data, x="Gene", y="Expression", color = "Type", 
           ylab="Gene expression",
           xlab=x,
           legend.title=x,
           add.params = list(fill="white"),
           palette = c("#4DBBD5","#E64B35"),
           width=1, add = "boxplot")
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Type),
                        method="t.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
                        label = "p.signif")

#输出
pdf(file=outFile, width=7, height=4)
print(p1)
dev.off()

