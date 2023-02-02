library(tidyverse)
library(ggalluvial)




mut_df_anno_fun_ctype_df = read_tsv('mut_df_anno_fun_ctype_df.tsv')


mut_df_anno_fun_ctype_df = mut_df_anno_fun_ctype_df[c('DepMap_ID','Variant_Classification','TransactivationClass','StructureFunctionClass','Variant_annotation')]

# DepMap_ID	Variant_Classification	TransactivationClass	StructureFunctionClass	Variant_annotation
# <chr>	<chr>	<chr>	<chr>	<chr>
# ACH-000003	Nonsense_Mutation	NA	NA	damaging
# ACH-000003	Nonsense_Mutation	NA	NA	damaging
# ACH-000003	Nonsense_Mutation	NA	NA	damaging
# ACH-000003	Nonsense_Mutation	NA	NA	damaging
# ACH-000004	Missense_Mutation	non-functional	non-functional	other non-conserving

# draw sanky graph
p = ggplot(data = mut_df_anno_fun_ctype_df,
       aes(axis1 = Variant_Classification, axis2 = TransactivationClass,axis3=StructureFunctionClass)) +
  scale_x_discrete(limits = c("Variant_Classification", "TransactivationClass",'StructureFunctionClass'), expand = c(.2, .05)) +
  #xlab("Demographic") +
  geom_alluvium(aes(fill = Variant_Classification)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()
ggsave('tp53_sanky.pdf',dpi=50,width=10, height=8)




































