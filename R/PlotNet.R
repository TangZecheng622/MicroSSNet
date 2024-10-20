#'  网络可视化
#'  @param g 输入需要可视化的网络
#'  @param clu_method 网络聚类方法.
#'  @param tag 网络名
#'  @param df 微生物丰度文件，提供则按照平均丰度作为网络点的大小
#'  @param node.cluster 每个模块点数量阈值，大于这个阈值才上色
#'  @export

plotnetwork <- function(g,clu_method="cluster_fast_greedy",tag=tag,df = NULL,node.cluster = 10){


  # cols <- c("#DEB99B" ,"#5ECC6D", "#5DAFD9", "#7ED1E4", "#EA9527", "#F16E1D" ,"#6E4821", "#A4B423","#F7B6D2FF","#aec7e8ff",
  #           "#C094DF" ,"#DC95D8" ,"#326530", "#50C0C9", "#67C021" ,"#DC69AF", "#8C384F", "#30455C", "#F96C72","#5ED2BF")
  base_colors <- pal_d3("category20")(20)
  color_to_remove <- c("#7F7F7FFF","#C7C7C7FF")
  base_colors <- base_colors[!(base_colors %in% color_to_remove)]
  col_g <- "#C7C7C7FF"
  V(g)$membership <- modularity_igraph(g,clu_method)[[1]]
  # modu_sort <- V(g)$membership %>% table() %>% sort(decreasing = T)
  modu_sort <- V(g)$membership %>% table()
  # top_num <- 18
  modu_sort <- modu_sort[modu_sort>=node.cluster]

  new_colors <- colorRamp(base_colors)(seq(0, 1, length.out = length(modu_sort)))
  cols <- rgb(new_colors, maxColorValue = 255)


  modu_name <- names(modu_sort)

  modu_cols <- cols[1:length(modu_name)]
  names(modu_cols) <- modu_name
  V(g)$label <- NA
  ###color
  V(g)$color <- V(g)$membership

  V(g)$color[!(V(g)$color %in% modu_name)] <- col_g
  V(g)$color[(V(g)$color %in% modu_name)] <- modu_cols[match(V(g)$color[(V(g)$color %in% modu_name)],modu_name)]
  V(g)$frame.color <- V(g)$color
  ###size
  if(!is.null(df)){
    V(g)$size = rowMeans(df[V(g),])
    V(g)$size = log1p(V(g)$size)
    min_size <- min(V(g)$size)
    max_size <- max(V(g)$size)
    V(g)$size <- (V(g)$size - min_size) / (max_size - min_size)*3
  }else{
    V(g)$size = 1
  }


  E(g)$color <- col_g
  for ( i in modu_name){
    col_edge <- cols[which(modu_name==i)]
    otu_same_modu <-V(g)$name[which(V(g)$membership==i)]
    E(g)$color[(data.frame(as_edgelist(g))$X1 %in% otu_same_modu)&(data.frame(as_edgelist(g))$X2 %in% otu_same_modu)] <- col_edge
  }
  g1 <- delete_vertices(g, which(degree(g)==0) ) #去除度为0
  sub_net_layout1 <- layout_with_fr(g1, niter=999,grid = 'nogrid')
  ## 可视化并输出
  pdf(file = paste(tag,"_network_plot.pdf"), width = 8,height = 8)
  par(font.main=4)
  plot(g1,layout=sub_net_layout1, edge.color = E(g)$color)
  title(main = paste0('Nodes=',length(V(g1)$name),', ','Edges=',nrow(data.frame(as_edgelist(g1)))),cex.main=2)
  dev.off()

  curves <- autocurve.edges2(g)
  file_name <- paste(tag,"_network_plot2.pdf")
  sub_net_layout <- layout_with_fr(g, niter=999,grid = 'nogrid')

  ## 可视化并输出
  pdf(file = file_name, width = 8,height = 8)
  par(font.main=4)
  plot(g,layout=sub_net_layout, edge.color = E(g)$color,#vertex.size=2,
       vertex.frame.width=0.1,
       #vertex.label=nodes_CQ1$Name,
       vertex.label.cex=0.1,
       edge.curved=curves,
       #vertex.size=V(phage_network)$size,
       #vertex.size2	=V(phage_network)$size,
       vertex.frame.color="grey")

  title(main = paste0('Nodes=',length(V(g)$name),', ','Edges=',nrow(data.frame(as_edgelist(g)))),cex.main=2)

  dev.off()
}
autocurve.edges2 <-function (graph, start = 0.5)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}


