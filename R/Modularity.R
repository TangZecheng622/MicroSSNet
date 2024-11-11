#' Network Modularity Calculation
#'
#' This function calculates the modularity of a network and returns the membership of nodes and the modularity value.
#'
#' @param g An `igraph` object representing the network.
#' @param clu_method Character string specifying the network clustering method. Supported methods are "cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass".
#' @return A list containing the membership vector and the modularity value.
#' @importFrom igraph E cluster_walktrap cluster_edge_betweenness cluster_fast_greedy cluster_spinglass membership modularity

#' @export
modularity_igraph <- function(g, clu_method = "cluster_fast_greedy") {
  # Parameter checks
  if (!inherits(g, "igraph")) stop("Error: 'g' must be an 'igraph' object.")
  supported_methods <- c("cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass")
  if (!clu_method %in% supported_methods) {
    stop("Error: 'clu_method' must be one of ", paste(supported_methods, collapse = ", "), ".")
  }
  # Check if edge weights exist
  if (is.null(igraph::E(g)$weight)) {
    igraph::E(g)$weight <- 1  # Set default weight to 1
  }

  fc <- switch(clu_method,
               "cluster_walktrap" = igraph::cluster_walktrap(g, weights = igraph::E(g)$weight),
               "cluster_edge_betweenness" = igraph::cluster_edge_betweenness(g, weights = igraph::E(g)$weight),
               "cluster_fast_greedy" = igraph::cluster_fast_greedy(g, weights = igraph::E(g)$weight),
               "cluster_spinglass" = igraph::cluster_spinglass(g, weights = igraph::E(g)$weight))
  membership_r <- igraph::membership(fc)
  modularity_r <- igraph::modularity(g, membership_r)
  return(list(membership = membership_r, modularity = modularity_r))
}
