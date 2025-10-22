#' Calculate Network Properties
#'
#' This function calculates local and global properties of a given network.
#'
#' @param graph An `igraph` object representing the network.
#' @param clu_method Character, network clustering method. Default is "cluster_walktrap".
#' @param zipi Logical, whether to calculate ZiPi values. Default is FALSE.
#' @param RM Logical, whether to calculate relative modularity. Default is TRUE.
#' @param tag Character, network name. Default is "NA".
#' @param output Character, output file path. Default is "./zipi/".
#' @return A list containing local and global indicators.
#' @importFrom igraph delete_vertices degree transitivity betweenness closeness eigen_centrality page_rank coreness V E edge_density edge_connectivity average.path.length diameter centr_degree centr_betw centr_clo centr_eigen count_components average_local_efficiency erdos.renyi.game modularity
#' @importFrom dplyr select everything
#' @importFrom stats na.omit
#' @importFrom utils head
#' @export
node_properties <- function(graph, clu_method = "cluster_walktrap", zipi = FALSE, RM = TRUE, tag = "NA", output = "./zipi/",
                            remove_isolated = TRUE) {
  if (!inherits(graph, "igraph")) stop("Error: 'graph' must be an 'igraph' object.")
  if (!is.character(clu_method)) stop("Error: 'clu_method' must be a character string.")
  if (!is.logical(zipi)) stop("Error: 'zipi' must be a logical value.")
  if (!is.logical(RM)) stop("Error: 'RM' must be a logical value.")
  if (!is.character(output)) stop("Error: 'output' must be a character string specifying the directory.")

  # Remove vertices with zero degree
  graph <- igraph::delete_vertices(graph, which(igraph::degree(graph) == 0))


  #Make sure the node has a name attribute
  if (is.null(igraph::V(graph)$name)) {
    igraph::V(graph)$name <- as.character(seq_along(igraph::V(graph)))
  }
  # Local indicators
  graph.degree <- igraph::degree(graph)
  graph.locclu.coeff <- igraph::transitivity(graph, type = "local")
  graph.locclu.coeff[is.na(graph.locclu.coeff)] <- 0
  graph.bc <- igraph::betweenness(graph)
  graph.cc <- igraph::closeness(graph)
  graph.ec <- igraph::eigen_centrality(graph)$vector
  graph.pagerank <- igraph::page_rank(graph)$vector
  graph.kcore <- igraph::coreness(graph)

  # Global indicators
  graph.nodes <- length(igraph::V(graph))
  graph.edges <- length(igraph::E(graph))
  graph.aveclu.coeff <- igraph::transitivity(graph, type = "average")
  connectance <- igraph::edge_density(graph, loops = FALSE)
  edge.connectivity <- igraph::edge_connectivity(graph)
  graph.average.degree <- mean(graph.degree)
  graph.average.path <- igraph::average.path.length(graph,weights = NULL)
  graph.diameter <- igraph::diameter(graph, directed = FALSE, unconnected = TRUE, weights = NULL)
  graph.cen.degree <- igraph::centr_degree(graph)$centralization
  graph.cen.bet <- igraph::centr_betw(graph)$centralization
  graph.cen.clo <- igraph::centr_clo(graph)$centralization
  graph.cen.eigen <- igraph::centr_eigen(graph)$centralization
  no.clusters <- igraph::count_components(graph)
  graph.aveloc.eff <- igraph::average_local_efficiency(graph)

  # Calculate Relative Modularity
  #计算相对模块化
  calculate_relative_modularity <- function(graph, clu_method) {
    mod1 <- modularity_igraph(graph, clu_method = clu_method)[[2]]
    rand.g <- igraph::erdos.renyi.game(length(igraph::V(graph)), length(igraph::E(graph)), type = "gnm")
    mod2 <- modularity_igraph(rand.g, clu_method = clu_method)[[2]]
    (mod1 - mod2) / mod2
  }

  if (RM == TRUE) {
    RM <- calculate_relative_modularity(graph, clu_method)
  } else {
    RM <- NA
  }

  # Summarize indicators
  graph.local.indicators <- data.frame(
    name = igraph::V(graph)$name,
    Degree = graph.degree,
    Local_Clustering_Coefficient = graph.locclu.coeff,
    Betweenness_Centrality = graph.bc,
    Closeness_Centrality = graph.cc,
    Eigenvector_Centrality = graph.ec,
    PageRank = graph.pagerank,
    K_core = graph.kcore
  )

  graph.global.indicators <- data.frame(
    Name = tag,
    Nodes = graph.nodes,
    Edges = graph.edges,
    Average_Local_Clustering_Coefficient = graph.aveclu.coeff,
    Connectance = connectance,
    Edge_Connectivity = edge.connectivity,
    Average_Degree = graph.average.degree,
    Average_Shortest_Path_Length = graph.average.path,
    Diameter = graph.diameter,
    Degree_Centralization = graph.cen.degree,
    Betweenness_Centralization = graph.cen.bet,
    Closeness_Centralization = graph.cen.clo,
    Eigenvector_Centralization = graph.cen.eigen,
    No_Cluster = no.clusters,
    Average_Local_Efficiency = graph.aveloc.eff,
    Relative_Modularity = RM
  )

  # Calculate ZiPi if required
  if (zipi == TRUE) {
    zipi_m <- ZiPiplot(graph, clu_method = clu_method, tag = tag, output_dir = output)
    graph.local.indicators <- dplyr::left_join(graph.local.indicators,zipi_m,by = "name")

    graph.local.indicators <- dplyr::select(graph.local.indicators, name, everything())
  }

  return(list(graph.local.indicators, graph.global.indicators))
}
