

data3 <- reshape2::melt(phosph_df, id.vars = c("sample"))

df_for_cal_phosph_site = merge(group_df,data3)%>%filter(!is.na(value))

kw_res_df_phosph_site = df_for_cal_phosph_site%>%filter() %>% 
  nest(-variable) %>% 
  mutate(kw=map(data,~kruskal.test(.x$value ~ .x$group))) %>%
  mutate(tidied = map(kw, tidy)) %>% 
  unnest(tidied, .drop = T)
kw_res_df_phosph_site[,c(-2,-3)]




