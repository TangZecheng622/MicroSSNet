#'  Zipi值分布可视化
#'  @param igraph 输入需要可视化的网络
#'  @param clu_method 网络聚类方法
#'  @param output 输出路径
#'  @param tag 网络名字
#'  @export
ZiPiplot = function(igraph = igraph,clu_method = "cluster_fast_greedy",output = "./",tag=NULL){

  comm_membership <- modularity_igraph(igraph,clu_method)[[1]]

  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
  }

  # 计算每个节点的度
  node_degree <- degree(igraph)

  # 计算每个模块的平均度和标准差
  module_mean <- tapply(node_degree, comm_membership, mean)
  module_sd <- tapply(node_degree, comm_membership, sd)

  # 创建一个空向量，储存每个节点与在同一模块的其他节点的度
  node_degree_module <- numeric()


  # 创建一个矩阵来存储每个节点与每个模块的连接数
  module_connections <- matrix(0, nrow=length(V(igraph)), ncol=max(comm_membership))

  # 计算每个节点与每个模块的连接数
  for (v in V(igraph)) {
    neighbors_of_v <- igraph::neighbors(igraph, v)
    modules_of_neighbors <- comm_membership[neighbors_of_v]
    node_degree_module[names(V(igraph))[v]] <- sum(modules_of_neighbors == comm_membership[v])
    module_connections[v, ] <- table(factor(modules_of_neighbors, levels=1:max(comm_membership)))
  }

  # 计算Zi，并处理标准差为0的情况
  Zi <- ifelse(module_sd[comm_membership] != 0 & !is.na(module_sd[comm_membership]),
               (node_degree_module - module_mean[comm_membership]) / module_sd[comm_membership],
               0)

  # 计算Pi
  Pi <- 1 - rowSums((module_connections / node_degree)^2)

  # 将结果整合到一个数据框
  zi_pi_metrics <- data.frame(
    name = V(igraph)$name,
    Zi = Zi,
    Pi = Pi
  )

  # 节点分类
  zi_pi_metrics$type <- ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi < 0.62, "Module hubs",
                               ifelse(zi_pi_metrics$Zi < 2.5 & zi_pi_metrics$Pi > 0.62, "Connectors",
                                      ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi > 0.62, "Network hubs",
                                             "Peripherals")))

  # 导出CSV
  write.csv(zi_pi_metrics, paste(output,tag,"_zi_pi_metrics.csv"), row.names = FALSE)

  p <- ggplot(zi_pi_metrics, aes(x = Pi, y = Zi, color = type)) +
    geom_point(alpha = 0.7, size = 3) +
    theme_minimal() +
    labs(
      x = "Among-module conectivities (Pi)",
      y = "Within-module conectivities (Zi)",
      color = "Node Type"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    geom_vline(xintercept = 0.62, linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = 2.5, linetype = "dashed", alpha = 0.5)

  ggsave(filename= paste(output,tag,"_zi_pi.pdf"),plot=p,width = 8,height = 8,dpi = 800)

  return(zi_pi_metrics)

}
