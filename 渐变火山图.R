rm(list = ls())
# 加载包：
library(ggplot2)

# 读取数据：
#data <- read.csv("data02.csv",row.names = 1)
load("E:/麻醉/03_差异分析/deg.Rdata")
load("E:/麻醉/03_差异分析/deg_Pro.Rdata")
deg1<-deg2
data<-deg1
# 新增一列用于存储label信息，将需要显示的label列出即可：
data$label <- c(rownames(data)[1:10],rep(NA,(nrow(data)-10)))

ggplot(data,aes(logFC, -log10(adj.P.Val)))+
  # 横向水平参考线：
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  # 纵向垂直参考线：
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  # 散点图:
  geom_point(aes(size=-log10(adj.P.Val), color= -log10(adj.P.Val)))+
  # 指定颜色渐变模式：
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  # 指定散点大小渐变模式：
  scale_size_continuous(range = c(1,3))+
  # 主题调整：
  theme_bw()+
  theme(panel.grid = element_blank())
