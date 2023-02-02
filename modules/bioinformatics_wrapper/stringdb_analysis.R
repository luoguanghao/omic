
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

library(STRINGdb)
library(igraph)
library(ggraph)

library(visNetwork)

# 创建STRINGdb对象
string_db <- STRINGdb$new( version="11", species=9606,  score_threshold=400, input_directory="")

get_stringdb_network <- function(gene, string_db){
    # 将Gene Symbol转换为Entrez ID
    gene <- gene %>% bitr(fromType = "SYMBOL", 
                          toType = "ENTREZID", 
                          OrgDb = "org.Hs.eg.db", 
                          drop = T)

    data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                    removeUnmappedRows = TRUE)
    
    # string_db$plot_network( data_mapped$STRING_id )

    data_links <- data_mapped$STRING_id %>% string_db$get_interactions()

    # 转换stringID为Symbol，只取前两列和最后一列
    links <- data_links %>%
      mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
      mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
      dplyr::select(from, to , last_col()) %>% 
      dplyr::rename(weight = combined_score)
    # 节点数据
    nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()
    # 创建网络图
    # 根据links和nodes创建
    net <- igraph::graph_from_data_frame(d=links,vertices=nodes,directed = F)
    # 添加一些参数信息用于后续绘图
    # V和E是igraph包的函数，分别用于修改网络图的节点（nodes）和连线(links)
    igraph::V(net)$deg <- igraph::degree(net) # 每个节点连接的节点数
    igraph::V(net)$size <- igraph::degree(net)/5 #
    igraph::E(net)$width <- igraph::E(net)$weight/10

    p_ggraph = ggraph(net,layout = "linear", circular = TRUE)+
      geom_edge_arc(aes(edge_width=width), color = "lightblue", show.legend = F)+
      geom_node_point(aes(size=size), color="orange", alpha=0.7)+
      geom_node_text(aes(filter=deg>5, label=name), size = 1, repel = F)+
      scale_edge_width(range = c(0.2,1))+
      scale_size_continuous(range = c(1,10) )+
      guides(size=F)+
      theme_graph()
    
    return(list(net=net, data_mapped=data_mapped, p_ggraph=p_ggraph))
}


plot_visNetwork <- function(net,result_dir,pj_name){
    netdf <- as_data_frame(net,what = "both")
    nodedf <- netdf$vertices
    edagedf <- netdf$edges

    nodesize <- degree(net)

    ## 设置网络的节点数据
    Newnodes <- data.frame(id=nodedf$name, # 节点的id
                           #label = nodedf$label, # 节点的标签
                           #group = paste("Group",nodedf$Faction), # 节点的分组
                           title = nodedf$name,# tooltip of the node
                           size = 10+1.1**(nodesize)
                           )
    ## 设置网络的边数据
    Newedages <- data.frame(from = edagedf$from,#边的起点
                            to = edagedf$to#, # 边的终点
                            )

    #nodes <- data.frame(id = 1:5, group = c(rep("A", 2), rep("B", 3)))
    #edges <- data.frame(from = c(2,5,3,3), to = c(1,2,4,2))

    p_visNetwork = visNetwork(Newnodes, Newedages, width = "100%", background = "white") #%>% 
    
    visSave(p_visNetwork, file = sprintf("%s/%s_p_visNetwork.html",result_dir,pj_name)) 
    
    return(p_visNetwork)
}


