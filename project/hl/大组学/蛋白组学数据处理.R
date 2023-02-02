library(tidyverse)
library(plotly)
options (warn = -1)




accumulate_plot <- function(){


    df_PHOgroupsite = readxl::read_excel('./log2-PHOgroupsite-sum-ref0-medin-ref.xlsx')

    df_PHOgroupsite = df_PHOgroupsite %>% column_to_rownames('SAMPLE-name')

    data = df_PHOgroupsite
    # ================
    long_data = data %>% gather(key = "variables", value = "values")
    long_data$values=as.numeric(long_data$values)



    long_data$gene = rep(rownames(data), dim(data)[2])

    data = data.frame(sapply(data, as.numeric), row.names=rownames(data))

    data_for_count = as.data.frame(is.na(data)==FALSE)
    # data_for_count

    for(i in 2:dim(data_for_count)[2]){
        data_for_count[[i]] = data_for_count[[i]]|data_for_count[[i-1]]
    }

    count_ls = sapply(data_for_count,sum)
    count_ls = count_ls[base::order(count_ls)]

    df_for_plot = data.frame(gene=names(count_ls),number=c(1:length(count_ls)),count=count_ls)


    p = ggplot() + geom_point(data=df_for_plot,aes(x=number, y=count)) + geom_line(data=df_for_plot,aes(x=number, y=count))


    return(p)


}







