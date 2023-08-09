rm(list=ls())
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
##dir.create("08_CIBERSORT")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))


library(dplyr)
library(tidyr)
library(tidyverse)

load("E:/麻醉/03_差异分析/差异分析_input.Rdata")
AS_exrp<-as.data.frame(t(exp))
load("E:/麻醉/08_CIBERSORT/CIBEROSTR.Rdata")
my_cibersort<-as.data.frame(my_cibersort)
cli<-cli_Pro
my_cibersort<-my_cibersort[rownames(cli),]
cibersort_raw<-my_cibersort

library(dplyr)
## 获取表达量数据
expr_data <- AS_exrp

immu_data <- cibersort_raw
#immu_data<-immu_data[,-1]

gene <- "XXXX"
y <- as.numeric(expr_data[,gene])

cor_data <- do.call(rbind,lapply(colnames(immu_data),function(x){
  dd <- cor.test(as.numeric(immu_data[,x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}))


### 画图展示全貌
library(dplyr)
library(ggplot2)
cor_data %>% 
  filter(p.value <0.05) %>% 
  ggplot(aes(cor,forcats::fct_reorder(cell,cor)))+
  geom_segment(aes(xend=0,yend=cell))+
  geom_point(aes(col=p.value,size=abs(cor)))+
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3"))+
  #scale_color_viridis_c(begin = 0.5, end = 1)+
  scale_size_continuous(range =c(2,8))+
  theme_bw()+
  ylab(NULL)

### 筛选p值有意义的细胞
imucells <- cor_data %>% 
  filter(p.value <0.05) %>% 
  arrange(desc(cor)) %>% 
  pull(cell) %>% 
  as.vector()

imucells

library(ggplot2)
corr_eqn <- function(x,y,digits=2) {
  test <- cor.test(x,y,method="spearman")
  paste(paste0("n = ",length(x)),
        paste0("r = ",round(test$estimate,digits),"(pearson)"),
        paste0("p.value= ",round(test$p.value,digits)),
        sep = ", ")
}

imucell_1 <- "T_cells_CD8"

plot_df_1 <- data.frame(
  gene = expr_data[,gene],
  imucell_1 = immu_data[,imucell_1]
)

plot_df_1<-as.data.frame(lapply(plot_df_1,as.numeric))


## 作图
plot_df_1 %>% 
  ggplot(aes(gene,imucell_1))+
  geom_point(col="#984ea3")+
  geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
  geom_rug(col="#7fc97f")+
  theme_minimal()+
  xlab(gene)+
  ylab(paste0(imucell_1," (NES)"))+
  labs(title = paste0(corr_eqn(plot_df_1$gene,plot_df_1$imucell_1)))+
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(1, 1, 1, 1, "cm"))


