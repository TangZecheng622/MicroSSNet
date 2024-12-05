#' Zi-Pi Plot Visualization
#'
#' This function visualizes the Zi-Pi distribution of nodes in a network and classifies them into roles.
#'
#' @param igraph An `igraph` object representing the network.
#' @param clu_method Character string specifying the network clustering method. Supported methods are "cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass".
#' @param output_dir Character string specifying the output directory. Defaults to current working directory.
#' @param tag Character string specifying the network name. Defaults to "group".
#' @return A data frame containing Zi and Pi values and node types.
#' @importFrom igraph V E neighbors degree
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs theme element_text geom_vline geom_hline ggsave
#' @export
ZiPiplot <- function(igraph, clu_method = "cluster_fast_greedy", output_dir = "./", tag = "group") {
  if (!inherits(igraph, "igraph")) stop("Error: 'igraph' must be an 'igraph' object.")
  if (igraph::gorder(igraph) == 0 || igraph::gsize(igraph) == 0) stop("Error: 'igraph' must have at least one node and one edge.")
  supported_methods <- c("cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass")
  if (!clu_method %in% supported_methods) {
    stop("Error: 'clu_method' must be one of ", paste(supported_methods, collapse = ", "), ".")
  }
  if (!is.character(output_dir)) stop("Error: 'output_dir' must be a character string.")
  if (!is.character(tag)) stop("Error: 'tag' must be a character string.")

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_file_csv <- file.path(output_dir, paste0(tag, "_zi_pi_metrics.csv"))
  output_file_pdf <- file.path(output_dir, paste0(tag, "_zi_pi.pdf"))

  comm_membership <- modularity_igraph(igraph, clu_method)$membership

  # Calculate node degrees
  node_degree <- igraph::degree(igraph)

  # Calculate mean and standard deviation of degrees within modules
  module_mean <- tapply(node_degree, comm_membership, mean)
  module_sd <- tapply(node_degree, comm_membership, sd)

  # Initialize vectors
  node_degree_module <- numeric(length(V(igraph)))
  module_connections <- matrix(0, nrow = length(V(igraph)), ncol = max(comm_membership))

  # Calculate node-degree within modules and module connections
  for (i in seq_along(V(igraph))) {
    v <- V(igraph)[i]
    neighbors_of_v <- igraph::neighbors(igraph, v)
    modules_of_neighbors <- comm_membership[neighbors_of_v]
    node_degree_module[i] <- sum(modules_of_neighbors == comm_membership[v])
    module_connections[i, ] <- tabulate(modules_of_neighbors, nbins = max(comm_membership))
  }

  # Calculate Zi, handling cases where module SD is zero
  Zi <- ifelse(module_sd[comm_membership] != 0 & !is.na(module_sd[comm_membership]),
               (node_degree_module - module_mean[comm_membership]) / module_sd[comm_membership],
               0)

  # Calculate Pi
  Pi <- 1 - rowSums((module_connections / node_degree)^2)

  # Combine results into a data frame
  zi_pi_metrics <- data.frame(
    name = V(igraph)$name,
    Zi = Zi,
    Pi = Pi
  )

  # Classify node types
  zi_pi_metrics$type <- ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi < 0.62, "Module hubs",
                               ifelse(zi_pi_metrics$Zi < 2.5 & zi_pi_metrics$Pi > 0.62, "Connectors",
                                      ifelse(zi_pi_metrics$Zi > 2.5 & zi_pi_metrics$Pi > 0.62, "Network hubs",
                                             "Peripherals")))

  # Export CSV
  write.csv(zi_pi_metrics, output_file_csv, row.names = FALSE)

  # Plotting
  p <- ggplot2::ggplot(zi_pi_metrics, ggplot2::aes(x = Pi, y = Zi, color = type)) +
    ggplot2::geom_point(alpha = 1 , size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Among-module connectivities (Pi)",
      y = "Within-module connectivities (Zi)",
      color = "Node Type"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    ggplot2::geom_vline(xintercept = 0.62, linetype = "dashed", alpha = 0.5) +
    ggplot2::geom_hline(yintercept = 2.5, linetype = "dashed", alpha = 0.5)

  ggplot2::ggsave(filename = output_file_pdf, plot = p, width = 8, height = 8, dpi = 800)

  # g <- ggplot2::ggplotGrob(p)
  # width <- grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE)
  # height <- grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE)
  #
  # ggplot2::ggsave("plot.png", plot = p, width = width, height = height, units = "in")

  return(zi_pi_metrics)
}
