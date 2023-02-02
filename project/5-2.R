# 做了几件事：
# 给画了热图，然后对于每个细胞，计算细胞IC与四个分子IC之间的多元统计关系

library(stringr)
library(ggpubr) 
library(ggplot2)



data_path = '/y/Bondi/jupyter/hl/5-2周老师PI3K多因子分析/PI3K相关化合物分子及细胞数据汇总.csv'


data <- read.csv(data_path, encoding="UTF-8",skip=2)
rownames(data) = data$X
data = na.omit(data[,c(2:9)])
rn  = rownames(data)

data = as.data.frame(lapply(data,as.character), stringsAsFactors=FALSE)

for(i in 1:dim(data)[1]){
    for(j in 1:dim(data)[2]){
        if(substr(data[i,j],1,1)=='>'){
            tmp=as.integer( str_trim(substr(data[i,j],2,nchar(as.character(data[i,j])))) )
            data[i,j]=tmp }
        if(stringr::str_detect(data[i,j],'±')){
            data[i,j]=unlist(strsplit(data[i,j],'±'))[1]
        }
    }

}
data = as.data.frame(lapply(data,as.numeric))
rownames(data) = rn

cell_df = data[,c(5,6,7,8)]
mol_df = data[,c(1,2,3,4)]

# heatmap
library(pheatmap)

p_res = pheatmap(log(data),clustering_method = "ward.D",clustering_distance_rows = "correlation")

# regression
## 做这些细胞系的IC数据和分子数据的多元分析  Kasumi.1,NCI.H929,RS411,HCC1954
cell_line = 'HCC1954'

lm_res = lm(get(cell_line)~`PI3Kα`+`PI3Kδ`+`PI3Kβ`+`PI3Kγ`,data) # 可以改

lm_res_df = list()
lm_res_df = as.data.frame(summary(lm_res)$coefficients )
lm_res_df$Var = rownames(lm_res_df)

lm_res_df

dataset = list()
dataset$Varnames = rownames(lm_res_df)
dataset$Mean = lm_res_df$Estimate
dataset$CI = lm_res_df$'Std. Error'
dataset$log10P = log10(lm_res_df$'Pr(>|t|)')
dataset = as.data.frame(dataset)[-1,]

p <- ggplot(dataset, aes(Mean, Varnames, col=log10P)) # 不同形状shape= Factor
p <- p + geom_point(size=5) +
  geom_errorbarh(aes(xmax = Mean + CI, xmin = Mean - CI), height = 0.1) +
  #scale_x_continuous(limits= c(0.1, 2.6), breaks= seq(0, 2.5, 0.5)) +
  geom_vline(aes(xintercept = 0),linetype="dotted") +
  xlab('Coefficient±Std. Error') + ylab(' ') +
  theme(axis.text.x =element_text(size=15), axis.text.y=element_text(size=15),axis.title.x =element_text(size=15))

plot(p)

for(i in 1:dim(dataset)[1]){
if (dataset[i,'log10P'] < log10(0.05)){
    #print(ggplot(data=data, aes(x=dataset[1,'Varnames'], y=HCC1954)) + geom_point())

    p=ggplot(data=log10(data.frame(x=data[,as.character(dataset[i,'Varnames'])],y=data[,cell_line])), 
        aes(x=x, y=y)) + geom_point() + 
        labs(title=dataset[i,'Varnames'], x=as.character(dataset[i,'Varnames']), y=cell_line)+
        geom_smooth(method = lm) +    
        stat_cor(method = "pearson")
    plot(p)
    }
    }








































































































































