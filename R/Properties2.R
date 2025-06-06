#'  网络属性计算
#' @param igraph `igraph`对象，待分析的网络图
#' @param clu_method 字符串，网络聚类方法，默认为 "cluster_walktrap"
#' @param zipi 逻辑值，是否计算Zipi值，默认为TRUE
#' @param RM 逻辑值，是否计算相对模块化，默认为TRUE
#' @param tag 字符串，网络名称，默认为NULL
#' @param output 字符串，输出文件路径，默认为 "./zipi/"
#' @return 返回包含局部和全局指标的列表
#' @importFrom igraph delete_vertices degree transitivity betweenness closeness eigen_centrality page_rank coreness V E edge_density edge_connectivity average.path.length diameter centr_degree centr_betw centr_clo centr_eigen count_components average_local_efficiency erdos.renyi.game modularity
#' @importFrom dplyr select everything
#' @importFrom stats na.omit
#' @importFrom utils head
#' @export
node_properties2<-function(igraph,clu_method="cluster_walktrap" ,zipi = FALSE,RM = TRUE,tag="NA",output = "./zipi/"){

  if (!inherits(igraph, "igraph")) stop("Error: 'igraph' must be an 'igraph' object.")
  if (!is.character(clu_method)) stop("Error: 'clu_method' must be a character string.")
  if (!is.logical(zipi)) stop("Error: 'zipi' must be a logical value.")
  if (!is.logical(RM)) stop("Error: 'RM' must be a logical value.")
  if (!is.character(output)) stop("Error: 'output' must be a character string specifying the directory.")
  #############################################局部指标(Local indicators)################################################
  igraph = delete_vertices(igraph, which(degree(igraph)==0) )

  ## Nodal degree度
  igraph.degree<-igraph::degree(igraph)

  ## Clustering coefficient(local)局部聚类系数
  igraph.locclu.coeff <- igraph::transitivity(igraph,type = "local")

  ##Betweenness Centrality(介数中心性)
  igraph.bc <- igraph::betweenness(igraph)

  ##Closeness Centrality (接近中心性)
  igraph.cc <- igraph::closeness(igraph)

  ##Eigenvector Centrality(特征向量中心性)
  igraph.ec <- igraph::eigen_centrality(igraph)$vector

  ##PageRank(页面排名)
  igraph.pagerank <- igraph::page_rank(igraph)$vector

  ##k-core Decomposition
  igraph.kcore <- igraph::coreness(igraph)

  ##degree distribution度分布
  igraph.degree.distribution <- igraph::degree_distribution(igraph)


  ################################################全局指标(Global indicators)##############################################


  ## N节点数
  igraph.nodes <- length(V(igraph))

  ## E边数量
  igraph.edges <- length(E(igraph))

  ## Average local clustering coefficient (平均局部聚类系数)
  igraph.aveclu.coeff <- igraph::transitivity(igraph,type = "average")

  ## Connectance (连通性)
  connectance <- edge_density(igraph,loops=FALSE)

  ## Edge Connectivity (边连通性)
  edge.connectivity <- edge_connectivity(igraph)

  ##Average Degree平均度
  igraph.average.degree <- mean(igraph.degree)

  ##Average Shortest Path Length平均最短路径长度
  igraph.average.path <- igraph::average.path.length(igraph)

  ##Diameter直径
  igraph.diameter <- igraph::diameter(igraph,directed = FALSE, unconnected = TRUE, weights = NULL)

  ##Degree centralization(度中心化)
  igraph.cen.degree<- igraph::centr_degree(igraph)$centralization

  ##Betweenness Centralization(介数中心化)
  igraph.cen.bet<-igraph::centr_betw(igraph)$centralization

  ##Closeness centralization (接近中心化)
  igraph.cen.clo<-igraph::centr_clo(igraph)$centralization

  ##Eigenvector centralization(特征向量中心化)
  igraph.cen.eigen <- igraph::centr_eigen(igraph)$centralization

  ##Assortativity (同配性)
  #  igraph.Assortativity <- igraph::assortativity(igraph)

  ##聚集cluster个数
  no.clusters <- count_components(igraph)

  #平均局部效率
  # igraph.aveloc.eff <- igraph::average_local_efficiency(igraph)

  #计算相对模块化
  calculate_relative_modularity <- function(igraph, clu_method) {
    mod1 <- modularity_igraph(igraph, clu_method = clu_method)[[2]]
    rand.g <- igraph::erdos.renyi.game(igraph::V(igraph), igraph::E(igraph), type = "gnm")
    mod2 <- modularity_igraph(rand.g, clu_method = clu_method)[[2]]
    (mod1 - mod2) / mod2
  }

  ############################################RM#############################################
  if(RM == TRUE){

    RM = calculate_relative_modularity(igraph, clu_method)

  }else{
    RM = 0
  }


  ############################################Summary###########################################
  # igraph.local.indicators <- cbind(igraph.degree,igraph.locclu.coeff,igraph.bc,igraph.cc,igraph.ec,igraph.pagerank,igraph.kcore)
  # colnames(igraph.local.indicators) <- c("Degree","Local_Clustering_Coefficient","Betweeness_Centrality","Closeness_Centrality","Eigenvector_Centrality","PageRank","K_core")
  igraph.local.indicators <- data.frame(
    Degree = igraph.degree,
    Local_Clustering_Coefficient = igraph.locclu.coeff,
    Betweeness_Centrality = igraph.bc,
    Closeness_Centrality = igraph.cc,
    Eigenvector_Centrality = igraph.ec,
    PageRank = igraph.pagerank,
    K_core = igraph.kcore
  )

  # igraph.global.indicators <- cbind(igraph.nodes,igraph.edges,igraph.aveclu.coeff,connectance,edge.connectivity,igraph.average.degree,igraph.average.path,
  #                                   igraph.diameter,igraph.cen.degree,igraph.cen.bet,igraph.cen.clo,igraph.cen.eigen,no.clusters,igraph.aveloc.eff,RM)
  # colnames(igraph.global.indicators)<-c("Nodes","Edges","Average_Local_Clustering_Coefficient","Connectance","Edge_Connectivity","Average_Degree","Average_Shortest_Path_Length","Diameter",
  # "Degree_Cenrealization","Betweenness_Centralization","Closeness_Centralization","Eigenvector_Centralization ","No_Cluster","Average_Local_Efficiency","Relative_Modularity")
  igraph.global.indicators <- data.frame(
    Name = tag,
    Nodes = igraph.nodes,
    Edges = igraph.edges,
    Average_Local_Clustering_Coefficient = igraph.aveclu.coeff,
    Connectance = connectance,
    Edge_Connectivity = edge.connectivity,
    Average_Degree = igraph.average.degree,
    Average_Shortest_Path_Length = igraph.average.path,
    Diameter = igraph.diameter,
    Degree_Centralization = igraph.cen.degree,
    Betweenness_Centralization = igraph.cen.bet,
    Closeness_Centralization = igraph.cen.clo,
    Eigenvector_Centralization = igraph.cen.eigen,
    No_Cluster = no.clusters,
    # Average_Local_Efficiency = igraph.aveloc.eff,
    Relative_Modularity = RM
  )

  ############################################ZiPi###########################################
  if (zipi==TRUE){
    zipi_m <- ZiPiplot(igraph,clu_method=clu_method,tag=tag,output_dir =  output)
    igraph.local.indicators1 <- cbind(igraph.local.indicators,zipi_m)
    igraph.local.indicators1 <- igraph.local.indicators1 %>% dplyr::select(name,everything())

  }else{
    igraph.local.indicators1 <- igraph.local.indicators
  }

  return(list(igraph.local.indicators1,igraph.global.indicators))
}
