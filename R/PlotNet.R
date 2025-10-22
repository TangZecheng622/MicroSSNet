#' Network Visualization
#'
#' This function visualizes the network with nodes colored by modules and sized by average abundance.
#'
#' @param g An `igraph` object representing the network.
#' @param clu_method Character string specifying the network clustering method.
#' @param tag Character string specifying the network name.
#' @param df Optional data frame containing species abundance data to size nodes.
#' @param node_cluster Integer, the minimum number of nodes in a module to assign a unique color.
#' @param output_dir Character string specifying the output directory.
#' @importFrom ggsci pal_d3
#' @importFrom scales rescale
#' @importFrom igraph as_edgelist delete_vertices degree V E layout_with_fr
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics par plot title
#' @export
plotnetwork <- function(g, clu_method = "cluster_fast_greedy", tag = "network", df = NULL, node_cluster = 10, output_dir = "./") {
  if (!inherits(g, "igraph")) stop("Error: 'g' must be an 'igraph' object.")
  if (!is.character(clu_method)) stop("Error: 'clu_method' must be a character string.")
  if (!is.character(tag)) stop("Error: 'tag' must be a character string.")
  if (!is.null(df) && !is.data.frame(df)) stop("Error: 'df' must be a data frame if provided.")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Color settings
  base_colors <- ggsci::pal_d3("category20")(20)
  color_to_remove <- c("#7F7F7FFF", "#C7C7C7FF")
  base_colors <- base_colors[!(base_colors %in% color_to_remove)]
  col_g <- "#C7C7C7FF"

  igraph::V(g)$membership <- modularity_igraph(g, clu_method)$membership
  modu_sort <- table(igraph::V(g)$membership)
  modu_sort <- modu_sort[modu_sort >= node_cluster]
  num_colors <- length(modu_sort)
  cols <- base_colors[seq_len(num_colors)]

  modu_name <- names(modu_sort)
  modu_cols <- cols
  names(modu_cols) <- modu_name

  igraph::V(g)$label <- NA
  # Assign colors to nodes
  igraph::V(g)$color <- ifelse(igraph::V(g)$membership %in% modu_name,
                       modu_cols[as.character(igraph::V(g)$membership)],
                       col_g)
  igraph::V(g)$frame.color <- igraph::V(g)$color
  node_sizes <- rowMeans(df[igraph::V(g)$name, , drop = FALSE], na.rm = TRUE)
  # Assign sizes to nodes
  if (!any(node_sizes< 0) ) {
    igraph::V(g)$size <- scales::rescale(log1p(node_sizes), to = c(0.5, 2))  # 标准化至 1-5 之间
  } else {
    igraph::V(g)$size <- scales::rescale(log1p(exp(node_sizes)), to = c(0.5, 2))  # 标准化至 1-5 之间
  }

  # Assign colors to edges
  igraph::E(g)$color <- col_g
  edge_list <- as.data.frame(igraph::as_edgelist(g))
  for (i in modu_name) {
    col_edge <- modu_cols[i]
    otu_same_modu <- igraph::V(g)$name[igraph::V(g)$membership == i]
    idx <- which(edge_list$V1 %in% otu_same_modu & edge_list$V2 %in% otu_same_modu)
    igraph::E(g)$color[idx] <- col_edge
  }

  # Remove isolated nodes
  g1 <- igraph::delete_vertices(g, which(igraph::degree(g) == 0))
  curves <- autocurve.edges2(g1)
  sub_net_layout1 <- igraph::layout_with_fr(g1, niter = 999, grid = 'nogrid')

  # Visualization and output
  output_file_pdf <- file.path(output_dir, paste0(tag, "_Network_Plot.pdf"))
  grDevices::pdf(file = output_file_pdf, width = 8, height = 8)
  par(font.main = 4)
  plot(g1,
       layout = sub_net_layout1,
       edge.curved = curves,
       edge.color = igraph::E(g1)$color,
       vertex.size = igraph::V(g1)$size,
       vertex.label = igraph::V(g1)$label,
       vertex.color = igraph::V(g1)$color,
       vertex.frame.color = igraph::V(g1)$frame.color)
  title(main = paste0('Nodes = ', length(igraph::V(g1)$name), ', Edges = ', length(igraph::E(g1))), cex.main = 2)
  grDevices::dev.off()
}


#' Adjust Edge Curvature for Visualization
#'
#' This function calculates the curvature of edges in a graph for better visualization.
#'
#' @param graph An `igraph` object.
#' @param start Numeric value indicating the starting curvature.
#' @return A numeric vector of curvature values for edges.
#' @importFrom igraph count_multiple is.mutual get.edgelist
#' @export
autocurve.edges2 <- function(graph, start = 0.5) {
  cm <- igraph::count_multiple(graph)
  mut <- igraph::is.mutual(graph)
  el <- apply(igraph::get.edgelist(graph, names = FALSE), 1, paste, collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut_obs <- mut[ord[p]]
    idx <- p:(p + m - 1)
    r <- if (m == 1 & !mut_obs) 0 else seq(-start, start, length.out = m)
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}
