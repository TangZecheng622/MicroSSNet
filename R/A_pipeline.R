#' Network Analysis Pipeline Function
#'
#' This function performs network construction and analysis based on microbial abundance data.
#'
#' @param table1 Data frame containing species abundance data (rows as species, columns as samples).
#' @param tax Optional data frame containing taxonomic information of species.
#' @param cor_table_list Optional list containing pre-calculated correlation matrices.
#' @param group_df Optional data frame containing sample grouping information.
#' @param vscol1 Character string specifying the column name in `group_df` that indicates group labels.
#' @param pvl.threshold Numeric value between 0 and 1 for prevalence threshold in filtering species.
#' @param top Numeric value. If `top >= 1`, selects the top `top` species based on average abundance.
#' If `0 < top < 1`, selects species with average abundance greater than `top`.
#' @param r.threshold Numeric value for the correlation coefficient threshold.
#' @param p.threshold Numeric value for the p-value threshold in correlation analysis.
#' @param method Character string specifying the correlation method ("pearson", "spearman", "sparcc", "SpiecEasi").
#' @param R Integer, number of bootstrap replicates for SparCC method.
#' @param ncpus Integer, number of CPU cores to use for parallel processing.
#' @param clu_method Character string specifying the network clustering method. Supported methods are "cluster_walktrap", "cluster_edge_betweenness", "cluster_fast_greedy", "cluster_spinglass".
#' @param p.adj Character string specifying the p-value adjustment method. Options are "holm", "hochberg", "hommel", "bonferroni", "BH" (Benjamini & Hochberg), "BY" (Benjamini & Yekutieli), "fdr", "none".
#' @param step Integer, number of steps in robustness analysis.
#' @param node.cluster Integer, minimum number of nodes in a module to assign a unique color in visualization.
#' @param rm.p.list Numeric vector of node removal proportions in robustness analysis.
#' @param zipi Logical, whether to calculate ZiPi values.
#' @param compare_net Logical, whether to compare the network with a random network.
#' @param calculate_vul Logical, whether to calculate vulnerability in robustness analysis.
#' @param calculate_rob Logical, whether to calculate robustness in robustness analysis.
#' @param calculate_cpx Logical, whether to calculate complexity in robustness analysis.
#' @param output_dir Character string specifying the output directory. Defaults to current working directory.
#'
#' @return The function outputs various files and plots in the specified output directory.
#' @export
#'
#' @examples
#' # Example usage:
#' otu_data <- data.frame(
#' OTU_ID = paste0("OTU_", 1:10),
#' CK1 = c(10, 20, 30, 30, 5, 6, 7, 8, 12, 10),
#' CK2 = c(12, 24, 28, 28, 6, 7, 8, 9, 11, 9),
#' CK3 = c(11, 22, 29, 27, 7, 5, 7, 8, 13, 11),
#' CK4 = c(10, 21, 30, 30, 6, 6, 8, 8, 12, 10),
#' CK5 = c(11, 23, 31, 31, 5, 6, 7, 7, 13, 12),
#' CK6 = c(12, 25, 28, 29, 6, 5, 6, 9, 11, 11),
#' PD1 = c(10, 20, 30, 30, 5, 6, 7, 8, 12, 10),
#' PD2 = c(12, 24, 28, 28, 6, 7, 8, 9, 11, 9),
#' PD3 = c(11, 22, 29, 27, 7, 5, 7, 8, 13, 11),
#' PD4 = c(10, 21, 30, 30, 6, 6, 8, 8, 12, 10),
#' PD5 = c(11, 23, 31, 31, 5, 6, 7, 7, 13, 12),
#' PD6 = c(12, 25, 28, 29, 6, 5, 6, 9, 11, 11)
#' )
#' rownames(otu_data) <- otu_data$OTU_ID
#' otu_data <- otu_data[, -1]
#'
#' aggregation_netpipeline (
#'   table1 = otu_data,
#'   tax = NULL,
#'   group_df = NULL,
#'   vscol1 = "group",
#'   pvl.threshold = 0,
#'   top = NULL,
#'   r.threshold = 0,
#'   p.threshold = 1,
#'   method = "pearson",
#'   ncpus = 1,
#'   output_dir = "./network_analysis_test"
#' )


aggregation_netpipeline <- function(
    table1,
    tax = NULL,
    cor_table_list = NULL,
    group_df = NULL,
    vscol1 = "Group",
    pvl.threshold = 0.5,
    top = NULL,
    r.threshold = 0.8,
    p.threshold = 0.05,
    method = "pearson",
    R = 10,
    ncpus = 4,
    clu_method = "cluster_fast_greedy",
    p.adj = "fdr",
    step = 10,
    node.cluster = 5,
    rm.p.list = seq(0, 0.5, by = 0.1),
    zipi = TRUE,
    compare_net = TRUE,
    calculate_vul = TRUE,
    calculate_rob = TRUE,
    calculate_cpx = TRUE,
    output_dir = getwd()
) {
  # Load required packages

  # Parameter validation
  if (!is.data.frame(table1)) stop("Error: 'table1' must be a data frame.")
  if (!is.null(tax) && !is.data.frame(tax)) stop("Error: 'tax' must be a data frame or NULL.")
  if (!is.null(cor_table_list) && !is.list(cor_table_list)) stop("Error: 'cor_table_list' must be a list or NULL.")
  if (!is.null(group_df) && !is.data.frame(group_df)) stop("Error: 'group_df' must be a data frame or NULL.")
  if (!is.character(vscol1)) stop("Error: 'vscol1' must be a character string.")
  if (!is.numeric(pvl.threshold) || pvl.threshold < 0 || pvl.threshold > 1) stop("Error: 'pvl.threshold' must be between 0 and 1.")
  if (!is.null(top) && (!is.numeric(top) || top <= 0)) stop("Error: 'top' must be a positive number or NULL.")
  if (!is.numeric(r.threshold)) stop("Error: 'r.threshold' must be numeric.")
  if (!is.numeric(p.threshold)) stop("Error: 'p.threshold' must be numeric.")
  if (!method %in% c("pearson", "spearman", "sparcc", "SpiecEasi")) stop("Error: Invalid 'method' selected.")
  if (!is.numeric(R) || R <= 0) stop("Error: 'R' must be a positive integer.")
  if (!is.numeric(ncpus) || ncpus <= 0) stop("Error: 'ncpus' must be a positive integer.")
  if (!is.character(clu_method)) stop("Error: 'clu_method' must be a character string.")
  if (!is.character(p.adj)) stop("Error: 'p.adj' must be a character string.")
  if (!is.numeric(step) || step <= 0) stop("Error: 'step' must be a positive integer.")
  if (!is.numeric(node.cluster) || node.cluster <= 0) stop("Error: 'node.cluster' must be a positive integer.")
  if (!is.numeric(rm.p.list) || any(rm.p.list < 0) || any(rm.p.list > 1)) stop("Error: 'rm.p.list' must be numeric values between 0 and 1.")
  if (!is.logical(zipi)) stop("Error: 'zipi' must be a logical value.")
  if (!is.logical(compare_net)) stop("Error: 'compare_net' must be a logical value.")
  if (!is.logical(calculate_vul)) stop("Error: 'calculate_vul' must be a logical value.")
  if (!is.logical(calculate_rob)) stop("Error: 'calculate_rob' must be a logical value.")
  if (!is.logical(calculate_cpx)) stop("Error: 'calculate_cpx' must be a logical value.")
  if (!is.character(output_dir)) stop("Error: 'output_dir' must be a character string.")

  # Create output directory if it doesn't exist
  create_output_dir(output_dir)

  # Initialize lists to store robustness results
  rrob_sum_list <- list()
  rrob_detail_list <- list()

  # Process group and table data
  check_result <- process_group_and_table(table = table1, group_df = group_df, vscol = vscol1)
  table1 <- check_result$table
  group_df <- check_result$group_df
  group_list <- check_result$group_list

  for (sel_group in group_list) {
    # sel_group <- group_list[[1]]
    message("Processing group: ", sel_group)
    # Create output directory if it doesn't exist
    group_output_dir <- file.path(output_dir, sel_group)
    create_output_dir(group_output_dir)
    # Extract abundance table for the selected group
    sel_samples <- group_df$sample[group_df$group == sel_group]
    sel_species_df <- table1[, sel_samples, drop = FALSE]

    # Filter species based on prevalence and average abundance
    sel_species_df <- sel_species_df[rowSums(sel_species_df > 0) > 0, , drop = FALSE]
    sel_species_df <- filter_OTU2(sel_species_df, Pre = pvl.threshold)
    sel_species_df <- filter_OTU(sel_species_df, Top = top)

    # Calculate correlation matrix
    if (!is.null(cor_table_list) && !is.null(cor_table_list[[sel_group]])) {
      cor <- cor_table_list[[sel_group]]

    } else {
      cor_result <- corMicro(
        table1 = sel_species_df,
        r.threshold = r.threshold,
        p.threshold = p.threshold,
        method = method,
        R = R,
        ncpus = ncpus,
        p.adj = p.adj,
        sel_group = sel_group,
        output_dir = group_output_dir
      )
      cor <- cor_result[[1]]
    }

    # Construct network graph
    g <- igraph::graph_from_adjacency_matrix(as.matrix(cor), weighted = TRUE, mode = 'undirected', diag = FALSE)
    igraph::E(g)$correlation <- igraph::E(g)$weight
    igraph::E(g)$weight <- abs(igraph::E(g)$weight)

    # Save network graph
    igraph::write_graph(g, file.path(group_output_dir, paste0(sel_group, '.Net.graphml')), format = 'graphml')
    igraph::write_graph(g, file.path(group_output_dir, paste0(sel_group, '.Net.gml')), format = 'gml')

    # Calculate network properties
    zipi_output_dir <- file.path(group_output_dir,"ZiPi")
    create_output_dir(zipi_output_dir)

    pro <- node_properties(g, clu_method = clu_method, zipi = zipi, tag = sel_group, output =  zipi_output_dir)
    local_pro <- pro[[1]]
    global_pro <- pro[[2]]

    # Save network properties
    write.csv(local_pro, file.path(group_output_dir, paste0(sel_group, "_Local_Properties.csv")), row.names = FALSE)
    write.csv(global_pro, file.path(group_output_dir, paste0(sel_group, "_Global_Properties.csv")), row.names = FALSE)

    # Network visualization
    if (!is.null(tax)) {
      # If taxonomic information is provided, use plotnetwork_tax function (implementation needed)
      message("Taxonomic information provided, but 'plotnetwork_tax' function is not implemented in this context.")
    } else {
      plotnetwork(g = g, clu_method = clu_method, tag = sel_group, df = sel_species_df, node_cluster = node.cluster, output_dir = group_output_dir)
    }

    # Compare with random network
    if (compare_net) {
      compare_rmc <- grobal_pro_compare(graph = g, step = step, netName = sel_group, ncpus = ncpus)
      compare_table <- compare_rmc[[1]]
      compare_p <- compare_rmc[[2]]
      write.csv(compare_table, file.path(group_output_dir, paste0(sel_group, "_Compare_RMC.csv")), row.names = TRUE)
      ggplot2::ggsave(filename = file.path(group_output_dir, paste0(sel_group, "_Degree_Distribution.pdf")),
                      plot = compare_p, width = 10, height = 6, dpi = 300)
    }

    # Robustness analysis
    if (any(c(calculate_vul, calculate_rob, calculate_cpx))) {
      random_rob <- Feature.Random.removal(
        cor = cor,
        otu = sel_species_df,
        tag = sel_group,
        rm.p.list = rm.p.list,
        nperm = step,
        ncpus = ncpus,
        calculate_vul = calculate_vul,
        calculate_rob = calculate_rob,
        calculate_cpx = calculate_cpx
      )

      rrob_sum_list[[paste0(sel_group, "_rrob_sum")]] <- random_rob$summary
      rrob_detail_list[[paste0(sel_group, "_rrob_detail")]] <- random_rob$detail
    }
  }

  # Compare robustness and vulnerability across groups

  compare_rob_vul(rrob_sum_list, rrob_detail_list, group_list, sel_group = paste(group_list, collapse = "_"), output_dir = output_dir)

}
