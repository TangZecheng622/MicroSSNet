#'  网络可视化，根据物种信息上色
#'  @param g 输入需要可视化的网络
#'  @param tag 网络名
#'  @param tax 物种分类信息
#'  @param df 微生物丰度文件，提供则按照丰度中位数作为网络点的大小
#'  @export



plotnetwork_tax <- function(g ,tag=tag,tax=tax,df = sel_species_df1){

  g <- delete_vertices(g, which(degree(g)==0) )
  V(g)$label <- NA
  ###设置边颜色，正相关红色，负相关蓝色，以及边的宽度
  df_edge=E(g)$correlation
  df_edge_color=ifelse(df_edge>0,"#FFE4E1",ifelse(df_edge<0,"#79CDCD","#EEE9E9"))
  E(g)$color=as.character(df_edge_color)
  E(g)$width = abs(df_edge)*1

  ###设置点大小###
  # V(g)$size = rowMeans(df[V(g),])
  V(g)$size = apply(df[V(g), ], 1, median)
  V(g)$size = log1p(V(g)$size)*2
  # min_size <- min(V(g)$size)
  # max_size <- max(V(g)$size)
  # V(g)$size <- (V(g)$size - min_size) / (max_size - min_size)*5

  ###设置点颜色###
  tax_df <- tax[V(g)$name,]
  for(i in colnames(tax_df)){

    tax_count <- table(tax_df[[i]])
    tax_count <- sort(tax_count,decreasing = TRUE)
    tax_df[[i]] <- factor(tax_df[[i]],levels = names(tax_count))
    g_col <- tax_df[[i]]
    ###生成颜色
    # set.seed(1)
    color <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423",
                      "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")
    levels(g_col)[1:18] = color
    levels(g_col)[19:length(g_col)] = "#C1C1C1"
    V(g)$color <- as.character(g_col)


    ###绘图
    # V(g)$membership <- modularity_igraph(g,"cluster_walktrap")[[1]]
    file_name <- paste(tag,"_",i,"_network_plot.pdf")
    sub_net_layout <- layout_with_fr(g, niter=999,grid = 'nogrid')

    ## 可视化并输出
    pdf(file = file_name, width = 10,height = 10)
    par(font.main=2)
    plot(g,layout=sub_net_layout, edge.color = E(g)$color,vertex.frame.color=NA)
    if(length(levels(tax_df[[i]]))>18){
      legend("bottomleft",levels(tax_df[[i]])[1:18],pch=21,col="black",pt.bg=color,cex = 0.8,box.lty=0, bg = "transparent",ncol = 2)
    }else{
      legend("bottomleft",levels(tax_df[[i]]),pch=21,col="black",pt.bg=color,cex = 0.8,box.lty=0, bg = "transparent",ncol = 2)
    }
    title(main = paste0('Nodes=',length(V(g)$name),', ','Edges=',nrow(data.frame(as_edgelist(g)))),cex.main=2)

    dev.off()
  }

}


#
# setwd("C:\\Users\\86135\\Desktop\\shi/")
#
#
# tax <- read.csv(file ="tax.csv",row.names = 1)
# otu <- read.csv(file = "SGB-shi.csv",row.names = 1)
# otu <- sweep(otu,2,colSums(otu),"/")
# rownames(otu) <- sub(";.*", "", rownames(otu))
#
# cor_pval <- corMicro(otu = otu, method = "spearman",R = 10,ncpus = 10)
# cor <- cor_pval[[1]]
# pval <- cor_pval[[2]]
# cor[abs(cor)<0.3|pval>0.05] = 0
# diag(cor) = 0
# g <-  graph_from_adjacency_matrix(cor, weighted = TRUE, mode = 'undirected')#构建邻接矩阵
#
# E(g)$correlation <- E(g)$weight #赋予正负相关性
# E(g)$weight <- abs(E(g)$weight) #不考虑正负的相关性
#
# V(g)$label <- NA
# V(g1)$membership <- modularity_igraph(g1,"cluster_fast_greedy")[[1]]
# g1 <- delete_vertices(g, which(degree(g)==0) )
#
# plotnetwork_tax(g=g,tag = "test",tax = tax,df = otu)
# plotnetwork_tax(g=g1,tag = "test2",tax = tax,df = otu)
