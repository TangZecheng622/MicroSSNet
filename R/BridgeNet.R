#' Construct and Visualize Bridge Network between Two Domains
#'
#' This function generates a bridge network visualization between two domains based on OTU data and correlation matrix.
#'
#' @param mix_otu_df Data frame, containing OTU data from both domains.
#' @param bircor Matrix, correlation matrix between the two domains.
#' @param table1name Character string specifying the name prefix for species in `table1`.
#' @param table2name Character string specifying the name prefix for species in `table2`.
#' @param output_path Character string specifying the directory path for saving output files.
#' @param sel_group Character string specifying the current group name.
#' @importFrom networktools bridge
#' @importFrom mgm mgm
#' @importFrom qgraph qgraph
#' @importFrom grDevices pdf dev.off
#' @importFrom stats setNames
#' @export
bridge_network <- function(
    mix_otu_df,
    bircor,
    table1name = "domain1",
    table2name = "domain2",
    output_path = getwd(),
    sel_group = "Group"
) {
  # Check if output_path exists
  if (!dir.exists(output_path)) stop("Error: 'output_path' does not exist.")

  # Load required packages
  if (!requireNamespace("networktools", quietly = TRUE)) {
    stop("Package 'networktools' is required but not installed.")
  }
  if (!requireNamespace("mgm", quietly = TRUE)) {
    stop("Package 'mgm' is required but not installed.")
  }
  if (!requireNamespace("qgraph", quietly = TRUE)) {
    stop("Package 'qgraph' is required but not installed.")
  }

  # 1. Data Processing
  bircor_filtered <- bircor[rowSums(abs(bircor)) > 0, colSums(abs(bircor)) > 0]
  otu_data <- mix_otu_df[, -ncol(mix_otu_df), drop = FALSE]
  otu_data <- otu_data[colnames(bircor_filtered), , drop = FALSE]

  # Convert OTU data to numeric matrix, remove NAs
  fit_matrix <- prepare_numeric_matrix(otu_data)

  # 2. Define Groups and Colors
  otu_positions <- grep(paste0("^", table1name), colnames(fit_matrix))
  sup_positions <- grep(paste0("^", table2name), colnames(fit_matrix))
  groups <- stats::setNames(list(otu_positions, sup_positions), c(table1name, table2name))
  group_colors <- stats::setNames(c("#B1CE46", "#9DC3E7"), c(table1name, table2name))

  # 3. Generate Main Network Plot
  pdf(file.path(output_path, paste0(sel_group, "_birnetwork_plot.pdf")), width = 7, height = 5)
  plot_main_network(bircor_filtered, groups, group_colors)
  dev.off()

  # 4. Bridge Network Construction
  bridge_nodes <- identify_bridge_nodes(bircor_filtered, table1name, table2name)
  sub_matrix <- mix_otu_df[bridge_nodes, , drop = FALSE]
  sub_network <- bircor[bridge_nodes, bridge_nodes]

  # Refined fit matrix for prediction
  refined_fit_matrix <- prepare_numeric_matrix(sub_matrix)
  refined_groups <- list(
    otu = grep(paste0("^", table1name), colnames(refined_fit_matrix)),
    sup = grep(paste0("^", table2name), colnames(refined_fit_matrix))
  )

  # 5. Predictive Analysis
  prediction <- perform_prediction(refined_fit_matrix)

  # 6. Bridge Centrality Analysis
  bridge_centrality <- networktools::bridge(
    sub_network,
    communities = refined_groups
  )

  # Plot bridge centrality
  pdf(file.path(output_path, paste0(sel_group, "_bridge_centrality_plot.pdf")), width = 6, height = 12)
  plot_bridge_centrality(bridge_centrality)
  dev.off()

  # 7. Bridge Network Visualization
  top_bridge_nodes <- get_top_bridge_nodes(bridge_centrality)
  groups_new <- define_new_groups(refined_groups$otu, refined_groups$sup, top_bridge_nodes, table1name, table2name)
  group_colors_extended <- c("#B1CE46", "#9DC3E7", "#f27970")
  names(group_colors_extended) <- c(table1name, table2name, "Bridge Nodes")

  pdf(file.path(output_path, paste0(sel_group, "_bridge_network_plot.pdf")), width = 7, height = 5)
  plot_bridge_network(sub_network, groups_new, group_colors_extended, prediction)
  dev.off()
}

# ---- Helper Functions ----

#' Prepare Numeric Matrix
#'
#' This function converts an OTU data frame to a numeric matrix, transposes it,
#' and removes any rows containing NA values.
#'
#' @param data Data frame containing OTU data, where each column represents a sample
#'        and each row represents an OTU feature.
#' @return A numeric matrix where rows correspond to OTU features, and columns correspond to samples.
#'         Any rows with NA values are removed.
#' @examples
#' # Example OTU data frame
#' otu_df <- data.frame(
#'     OTU1 = c(1, 2, NA),
#'     OTU2 = c(3, 4, 5),
#'     OTU3 = c(6, 7, 8)
#' )
#' prepare_numeric_matrix(otu_df)
#' @export
prepare_numeric_matrix <- function(data) {
  matrix_data <- as.matrix(t(data))
  matrix_data <- apply(matrix_data, 2, as.numeric)
  matrix_data <- na.omit(matrix_data)
  return(matrix_data)
}

#' Main Network Plot
#'
#' @param network Matrix, filtered correlation matrix.
#' @param groups List, positions of groups in the network.
#' @param colors Named vector of colors for each group.
#' @importFrom qgraph qgraph
plot_main_network <- function(network, groups, colors) {
  qgraph::qgraph(
    network,
    vsize = 1,
    layout = "spring",
    color = colors,
    groups = groups,
    labels = "",
    label.cex = 1,
    label.color = 'grey',
    negDashed = FALSE,
    legend = TRUE,
    legend.cex = 0.45,
    legend.mode = 'style3',
    border.color = NA
  )
}

#' Identify Bridge Nodes
#'
#' Extracts nodes directly related between two domains based on the correlation matrix.
#'
#' @param cor_matrix Matrix, correlation matrix between the two domains.
#' @param table1name Character, name prefix for species in the first domain.
#' @param table2name Character, name prefix for species in the second domain.
#' @return A vector of node names that act as bridges between the two domains.
#' @export
identify_bridge_nodes <- function(cor_matrix, table1name, table2name) {
  table1_positions <- grep(paste0("^", table1name), colnames(cor_matrix))
  table2_positions <- grep(paste0("^", table2name), colnames(cor_matrix))
  direct_related <- cor_matrix[table1_positions, table2_positions, drop = FALSE]
  direct_related <- direct_related[rowSums(abs(direct_related)) > 0, colSums(abs(direct_related)) > 0, drop = FALSE]
  return(unique(c(rownames(direct_related), colnames(direct_related))))
}

#' Perform Prediction Analysis
#'
#' Uses mgm to perform predictive analysis on a matrix.
#' @param data_matrix Matrix, the numeric data matrix used for prediction.
#' @importFrom mgm mgm
#' @importFrom stats predict
#' @return A vector of prediction error values.
#' @export
perform_prediction <- function(data_matrix) {
  p <- ncol(data_matrix)
  fit_obj <- mgm::mgm(data = data_matrix, type = rep("g", p), level = rep(1, p), lambdaSel = 'EBIC', ruleReg = 'OR', pbar = TRUE)
  pred_obj <- stats::predict(object = fit_obj, data = data_matrix, errorCon = 'R2')
  return(pred_obj$error[, 2])
}

#' Plot Bridge Centrality
#'
#' Visualizes bridge centrality using specified metrics.
#' @importFrom networktools bridge
#' @param bridge_centrality List, bridge centrality analysis results.
#' @export
plot_bridge_centrality <- function(bridge_centrality) {
  plot(bridge_centrality, include = c("Bridge Strength"), order = "value", zscore = TRUE)
}

#' Get Top Bridge Nodes
#'
#' Extracts nodes with top bridge strength.
#' @param bridge_centrality List, bridge centrality analysis result containing bridge strength values.
#' @importFrom stats quantile
get_top_bridge_nodes <- function(bridge_centrality) {
  bridge_strength <- bridge_centrality$`Bridge Strength`
  threshold <- stats::quantile(bridge_strength, probs = 0.90, na.rm = TRUE)
  top_bridges <- names(bridge_strength[bridge_strength > threshold])
  return(which(names(bridge_strength) %in% top_bridges))
}

#' Define New Groups for Visualization
#'
#' Creates new group definitions with bridge nodes.
#' @param otu_positions Integer vector, positions of OTU nodes in the network.
#' @param sup_positions Integer vector, positions of supplementary nodes in the network.
#' @param bridge_nodes Integer vector, positions of bridge nodes connecting the two domains.
#' @param table1name Character, name of the first domain group.
#' @param table2name Character, name of the second domain group.
#' @return A named list of node positions for each group.
#' @export
define_new_groups <- function(otu_positions, sup_positions, bridge_nodes, table1name, table2name) {
  filtered_otu <- otu_positions[!otu_positions %in% bridge_nodes]
  filtered_sup <- sup_positions[!sup_positions %in% bridge_nodes]
  return(stats::setNames(list(filtered_otu, filtered_sup, bridge_nodes), c(table1name, table2name, "Bridge Nodes")))
}

#' Bridge Network Plot
#'
#' @param network Matrix, sub-network for bridge analysis.
#' @param groups List, positions of nodes in each group.
#' @param colors Named vector of colors for each group.
#' @param prediction Vector, prediction error values for pie charts.
#' @importFrom qgraph qgraph
plot_bridge_network <- function(network, groups, colors, prediction) {
  qgraph::qgraph(
    network,
    vsize = 3,
    layout = "spring",
    color = colors,
    groups = groups,
    labels = colnames(network),
    label.cex = 1,
    label.color = 'gray',
    negDashed = FALSE,
    legend = TRUE,
    legend.cex = 0.5,
    legend.mode = 'style3',
    repulsion = 0.9,
    pie = prediction,
    border.color = NA
  )
}
