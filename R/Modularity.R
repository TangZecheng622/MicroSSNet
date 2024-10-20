#'  网络模块化
#'  @param g 输入需要模块化的网络
#'  @param clu_method 网络聚类方法
#'  @export

modularity_igraph = function(g,clu_method = "cluster_fast_greedy"){
  if (clu_method == "cluster_walktrap" ) {
    fc <- igraph::cluster_walktrap(g,weights =  igraph::E(g)$weight)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }

  if (clu_method == "cluster_edge_betweenness" ) {
    fc <- igraph::cluster_edge_betweenness(g,weights =  igraph::E(g)$weight)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (clu_method == "cluster_fast_greedy" ) {
    fc <- igraph::cluster_fast_greedy(g,weights =  igraph::E(g)$weight)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  if (clu_method == "cluster_spinglass" ) {
    fc <- igraph::cluster_spinglass(g,weights =  igraph::E(g)$weight)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
  }
  membership_r <- membership(fc)
  modularity_r <- igraph::modularity(g,membership_r)
  return(list(membership_r,modularity_r))
}

