library(tidyverse)
library(ggplot2)






corr_df = read_csv('/mnt/d/my_project/hl/AML大组学/数据/cor_pro_rna.csv')

##############
# hist plot
##############
p <- ggplot(corr_df,aes(x=cor)) + 
  geom_histogram(aes(y=..density.., fill=..x..),color='black',size=1) +
  geom_density(color='black', size=1.5) +
  scale_fill_gradient2("Count", low = "gold", high = "navy") +
  theme_bw() +
  labs(x="Pearson's Correlation Coefficients", y='Probability DEnsity') +
  geom_vline(aes(xintercept=median(cor)), color='black',size=1, linetype="dashed")

ggsave('histogram.pdf', dpi=50, width=7, height=7)
plot(p)

##############
# rank plot
##############
corr_df_sorted = corr_df[base::order(corr_df$cor,decreasing = TRUE),]
corr_df_sorted$X1 = factor(corr_df_sorted$X1,level=corr_df_sorted$X1)

# Barplot
ggplot(corr_df_sorted) + 
    geom_bar(aes(x=X1, y=cor, color=cor), stat = "identity", width=0.2) +
    scale_colour_gradient2(low = "Gold", mid = "white",high = "Navy",midpoint = 0) +
    labs(x='Protein No.',y='Correlation') +
    theme_few() +
    theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

ggsave('rank.pdf', dpi=50, width=7, height=3)



##############
# C图
## KEGG
##############
library("clusterProfiler")
library(ggthemes)

gmtfile <- '/mnt/d/my_project/data/pathway/c2.cp.kegg.v7.2.symbols.gmt'
c5 <- read.gmt(gmtfile)

df_for_plot = merge(corr_df,c5,by.x='X1',by.y='gene')
df_for_plot = df_for_plot%>%filter(term%in%c('KEGG_CITRATE_CYCLE_TCA_CYCLE','KEGG_PENTOSE_PHOSPHATE_PATHWAY',
                                    'KEGG_GLYCOLYSIS_GLUCONEOGENESIS','KEGG_PROTEIN_EXPORT'))


chosen_pw = list(pos_pw=c('KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION','KEGG_FC_GAMMA_R_MEDIATED_PHAGOCYTOSIS','KEGG_CHEMOKINE_SIGNALING_PATHWAY' ),
                          back_pw=c('KEGG_SPLICEOSOME','KEGG_OXIDATIVE_PHOSPHORYLATION','KEGG_RIBOSOME'))



df_for_plot = merge(corr_df,c5,by.x='X1',by.y='gene')
df_for_plot = df_for_plot%>%filter(term%in%c(chosen_pw$pos_pw, chosen_pw$back_pw))



cols <- c('TRUE' = "navy", 'FALSE' = "Gold")
ggplot(df_for_plot ,aes(x=cor, y=term, fill=term)) +
    #geom_violin() +
    geom_point(aes(color=cor>0), shape='|',size=10,position = position_dodge(0)) +
    theme_bw() +
    scale_colour_manual(values = cols) +
    theme(legend.position = 'none') +
    ylab('') + xlab("Pearson's Correlation Coefficients")  +
    scale_x_continuous(breaks=seq(-0.5, 1, 0.25))

ggsave('enrichment.pdf', dpi=50, width=11, height=3)












