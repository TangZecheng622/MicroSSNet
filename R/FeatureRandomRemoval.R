#' Robustness Analysis by Random Node Removal
#'
#' This function simulates random node removal to assess network robustness, vulnerability, and complexity.
#'
#' @param cor Numeric matrix representing the correlation matrix of the network.
#' @param otu Data frame of species abundance (rows as species, columns as samples).
#' @param tag Character string for labeling output files.
#' @param nperm Integer, number of permutations for simulation. Default is 100.
#' @param rm.p.list Numeric vector of node removal proportions. Default is seq(0.05, 1, by = 0.05).
#' @param ncpus Integer, number of CPU cores for parallel processing.
#' @param calculate_vul Logical, whether to calculate network vulnerability.
#' @param calculate_rob Logical, whether to calculate network robustness.
#' @param calculate_cpx Logical, whether to calculate network complexity.
#' @return A list containing summary and detailed data frames of the simulation results.
#' @export
Feature.Random.removal <- function(
    cor,
    otu,
    tag = "default",
    nperm = 100,
    rm.p.list = seq(0.05, 1, by = 0.05),
    ncpus = 2,
    calculate_vul = TRUE,
    calculate_rob = TRUE,
    calculate_cpx = TRUE,
    seed = 1
) {
  # Parameter checks
  if (!is.matrix(cor) && !is.data.frame(cor)) stop("Error: 'cor' must be a matrix or data frame.")
  if (!is.data.frame(otu)) stop("Error: 'otu' must be a data frame.")
  if (!is.character(tag)) stop("Error: 'tag' must be a character string.")
  if (!is.numeric(nperm) || nperm <= 0) stop("Error: 'nperm' must be a positive integer.")
  if (!is.numeric(rm.p.list) || any(rm.p.list < 0) || any(rm.p.list > 1)) stop("Error: 'rm.p.list' values must be between 0 and 1.")
  if (!is.numeric(ncpus) || ncpus <= 0) stop("Error: 'ncpus' must be a positive integer.")

  # Prepare correlation matrix and abundance weights
  cor1 <- cor
  diag(cor1) <- 0
  sp.ra <- rowMeans(otu)
  cor1_raw <- cor1[rowSums(abs(cor1)) > 0, colSums(abs(cor1)) > 0]
  sp.ra2 <- sp.ra[rowSums(abs(cor1)) > 0]
  options(future.globals.maxSize = 4 * 1024^3)
  # Run simulation
  results <- rmsimu3(
    netRaw = cor1_raw,
    rm.p.list = rm.p.list,
    sp.ra = sp.ra2,
    abundance.weighted = TRUE,
    nperm = nperm,
    ncpus = ncpus,
    calculate_vul = calculate_vul,
    calculate_rob = calculate_rob,
    calculate_cpx = calculate_cpx,
    tag = tag,
    seed = seed
  )

  return(results)
}

#' Simulate Random Node Removal
#'
#' This helper function simulates random node removal to assess network robustness, vulnerability, and complexity.
#'
#' @param netRaw Numeric matrix representing the network correlation matrix.
#' @param rm.p.list Numeric vector of node removal proportions.
#' @param sp.ra Numeric vector of species abundances.
#' @param abundance.weighted Logical, whether to weight interactions by abundance.
#' @param nperm Integer, number of permutations.
#' @param ncpus Integer, number of CPU cores for parallel processing.
#' @param calculate_vul Logical, whether to calculate vulnerability.
#' @param calculate_rob Logical, whether to calculate robustness.
#' @param calculate_cpx Logical, whether to calculate complexity.
#' @param tag Character string for labeling output.
#' @return A list containing summary and detailed data frames of the simulation results.

rmsimu3 <- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm = 10, ncpus = 2,
                    calculate_vul = TRUE, calculate_rob = TRUE, calculate_cpx = TRUE, tag = "default",seed = 1) {
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install 'future.apply'")

  Vulnerability <- function(cor) {
    # Parameter checks
    if (!is.matrix(cor) && !is.data.frame(cor)) stop("Error: 'cor' must be a matrix or data frame.")
    if (nrow(cor) != ncol(cor)) stop("Error: 'cor' must be a square matrix.")

    cor1 <- as.matrix(cor)
    diag(cor1) <- 0
    cor1[abs(cor1) > 0] <- 1
    g <- igraph::graph_from_adjacency_matrix(cor1, mode = "undirected", weighted = NULL, diag = FALSE)

    # Remove isolated nodes
    iso_node_id <- which(igraph::degree(g) == 0)
    if (length(iso_node_id) > 0) {
      g <- igraph::delete_vertices(g, iso_node_id)
    }

    # Calculate network efficiency
    network_efficiency <- function(graph) {
      if (!igraph::is_igraph(graph)) stop("Error: Please provide a valid igraph object.")
      dist_inv <- 1 / igraph::distances(graph)
      diag(dist_inv) <- NA
      mean(dist_inv, na.rm = TRUE)
    }

    net_eff <- network_efficiency(g)

    # Calculate vulnerability for each node
    info_centrality_vertex <- function(graph, net_eff) {
      if (!igraph::is_igraph(graph)) stop("Error: Please provide a valid igraph object.")
      n <- igraph::vcount(graph)
      vulnerabilities <- numeric(n)
      for (i in seq_len(n)) {
        g_removed <- igraph::delete_vertices(graph, i)
        eff_removed <- network_efficiency(g_removed)
        vulnerabilities[i] <- (net_eff - eff_removed) / net_eff
      }
      return(vulnerabilities)
    }

    node_vul <- info_centrality_vertex(g, net_eff)

    # Return the maximum node vulnerability
    max_vul <- max(node_vul, na.rm = TRUE)
    return(max_vul)
  }

  boot_feature <- function(i, netRaw, rm.percent, sp.ra, abundance.weighted,
                           calculate_vul = TRUE, calculate_rob = TRUE, calculate_cpx = FALSE) {
    # Randomly remove a proportion of nodes
    id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
    net.Modified <- netRaw
    net.Modified[id.rm, ] <- 0
    net.Modified[, id.rm] <- 0

    # Apply abundance weighting if enabled
    net.strength <- if (abundance.weighted) net.Modified * sp.ra else net.Modified
    sp.meanInteraction <- rowMeans(net.strength)
    id.rm2 <- which(sp.meanInteraction <= 0)

    remain.percent <- if (calculate_rob) (nrow(netRaw) - length(id.rm2)) / nrow(netRaw) else NA
    vul <- if (calculate_vul) Vulnerability(net.strength) else NA
    cpx <- if (calculate_cpx) {
      gi <- igraph::graph_from_adjacency_matrix(net.strength, weighted = TRUE, mode = 'undirected', diag = FALSE)
      length(igraph::E(gi)) / length(igraph::V(gi))
    } else NA

    return(list(remain = remain.percent, vul = vul, complexity = cpx))
  }
  future::plan(future::multisession, workers = ncpus)

  dat1_results <- data.frame()
  dat2_results <- data.frame()

  for (rm.percent in rm.p.list) {
    # 包装器，显式 set.seed 让每次迭代都随机但不重复
    boot_feature_wrapper <- function(i) {
      set.seed(seed + i)
      boot_feature(i, netRaw, rm.percent, sp.ra,
                   abundance.weighted, calculate_vul,
                   calculate_rob, calculate_cpx)
    }

    results_list <- future.apply::future_lapply(1:nperm, boot_feature_wrapper, future.seed = NULL)

    Proportion.remain <- sapply(results_list, `[[`, "remain")
    vulnerability     <- sapply(results_list, `[[`, "vul")
    complexity        <- sapply(results_list, `[[`, "complexity")

    dat1_results <- rbind(dat1_results, data.frame(
      group = tag,
      Proportion.removed = rm.percent,
      Robustness.mean = mean(Proportion.remain, na.rm = TRUE),
      Robustness.sd = sd(Proportion.remain, na.rm = TRUE),
      Vulnerability.mean = mean(vulnerability, na.rm = TRUE),
      Vulnerability.sd = sd(vulnerability, na.rm = TRUE),
      Complexity.mean = mean(complexity, na.rm = TRUE),
      Complexity.sd = sd(complexity, na.rm = TRUE)
    ))

    dat2_results <- rbind(dat2_results, data.frame(
      group = rep(tag, nperm),
      iteration = 1:nperm,
      Proportion.removed = rep(rm.percent, nperm),
      Robustness = Proportion.remain,
      Vulnerability = vulnerability,
      Complexity = complexity
    ))
  }

  return(list(summary = dat1_results, detail = dat2_results))
}
# rmsimu3 <- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE, nperm = 10, ncpus = 2,
#                     calculate_vul = TRUE, calculate_rob = TRUE, calculate_cpx = TRUE, tag = "default") {
#   # Parameter checks
#   if (!is.matrix(netRaw)) stop("Error: 'netRaw' must be a matrix.")
#   if (!is.numeric(rm.p.list) || any(rm.p.list < 0) || any(rm.p.list > 1)) stop("Error: 'rm.p.list' values must be between 0 and 1.")
#   if (!is.numeric(sp.ra)) stop("Error: 'sp.ra' must be a numeric vector.")
#   if (!is.numeric(nperm) || nperm <= 0) stop("Error: 'nperm' must be a positive integer.")
#   if (!is.numeric(ncpus) || ncpus <= 0) stop("Error: 'ncpus' must be a positive integer.")
#
#   # Initialize cluster
#   cl <- parallel::makeCluster(ncpus)
#   safe_stop_cluster <- function(cl) {
#     tryCatch({
#       parallel::stopCluster(cl)
#     }, error = function(e) {
#       message("Warning: stopCluster failed. It may already be closed.")
#     })
#   }
#   on.exit(safe_stop_cluster(cl))
#
#
#   dat1_results <- data.frame()
#   dat2_results <- data.frame()
#
#   for (rm.percent in rm.p.list) {
#     # Export necessary variables
#     vars_to_export <- c("netRaw", "boot_feature", "rm.percent",
#                         "abundance.weighted", "sp.ra",
#                         "calculate_vul", "calculate_rob", "calculate_cpx")
#
#     if (calculate_vul) {
#       vars_to_export <- c(vars_to_export, "Vulnerability")
#     }
#
#     parallel::clusterEvalQ(cl, { library(igraph) })
#     parallel::clusterExport(cl, varlist = vars_to_export, envir = environment())
#
#     # Run simulations in parallel
#     results_list <- parallel::parLapply(cl, 1:nperm, function(i) {
#       local_netRaw <- netRaw
#       local_rm.percent <- rm.percent
#       local_sp.ra <- sp.ra
#       local_abundance.weighted <- abundance.weighted
#       local_calculate_vul <- calculate_vul
#       local_calculate_rob <- calculate_rob
#       local_calculate_cpx <- calculate_cpx
#
#       tryCatch({
#         boot_feature(i, local_netRaw, local_rm.percent, local_sp.ra,
#                      local_abundance.weighted, local_calculate_vul,
#                      local_calculate_rob, local_calculate_cpx)
#       }, error = function(e) {
#         message(sprintf("Iteration %d failed: %s", i, e$message))
#         return(list(remain = NA, vul = NA, complexity = NA))
#      })
#     })
#
#     # Extract results
#     Proportion.remain <- sapply(results_list, `[[`, "remain")
#     vulnerability <- sapply(results_list, `[[`, "vul")
#     complexity <- sapply(results_list, `[[`, "complexity")
#
#     # Compile summary results
#     dat1_results <- rbind(dat1_results, data.frame(
#       group = tag,
#       Proportion.removed = rm.percent,
#       Robustness.mean = mean(Proportion.remain, na.rm = TRUE),
#       Robustness.sd = sd(Proportion.remain, na.rm = TRUE),
#       Vulnerability.mean = mean(vulnerability, na.rm = TRUE),
#       Vulnerability.sd = sd(vulnerability, na.rm = TRUE),
#       Complexity.mean = mean(complexity, na.rm = TRUE),
#       Complexity.sd = sd(complexity, na.rm = TRUE)
#     ))
#
#     # Compile detailed results
#     dat2_results <- rbind(dat2_results, data.frame(
#       group = rep(tag, nperm),
#       iteration = 1:nperm,
#       Proportion.removed = rep(rm.percent, nperm),
#       Robustness = Proportion.remain,
#       Vulnerability = vulnerability,
#       Complexity = complexity
#     ))
#   }
#
#   return(list(summary = dat1_results, detail = dat2_results))
# }

#' Bootstrapping Function for Random Node Removal
#'
#' This helper function performs a single iteration of random node removal simulation.
#'
#' @param i Integer, iteration index.
#' @param netRaw Numeric matrix representing the network correlation matrix.
#' @param rm.percent Numeric, proportion of nodes to remove.
#' @param sp.ra Numeric vector of species abundances.
#' @param abundance.weighted Logical, whether to weight interactions by abundance.
#' @param calculate_vul Logical, whether to calculate vulnerability.
#' @param calculate_rob Logical, whether to calculate robustness.
#' @param calculate_cpx Logical, whether to calculate complexity.
#' @return A list containing the remaining proportion, vulnerability, and complexity.
boot_feature <- function(i, netRaw, rm.percent, sp.ra, abundance.weighted,
                         calculate_vul = TRUE, calculate_rob = TRUE, calculate_cpx = FALSE) {
  # Randomly remove a proportion of nodes
  id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
  net.Modified <- netRaw
  net.Modified[id.rm, ] <- 0
  net.Modified[, id.rm] <- 0

  # Apply abundance weighting if enabled
  net.strength <- if (abundance.weighted) net.Modified * sp.ra else net.Modified
  sp.meanInteraction <- rowMeans(net.strength)
  id.rm2 <- which(sp.meanInteraction <= 0)

  remain.percent <- if (calculate_rob) (nrow(netRaw) - length(id.rm2)) / nrow(netRaw) else NA
  vul <- if (calculate_vul) Vulnerability(net.strength) else NA
  cpx <- if (calculate_cpx) {
    gi <- igraph::graph_from_adjacency_matrix(net.strength, weighted = TRUE, mode = 'undirected', diag = FALSE)
    length(igraph::E(gi)) / length(igraph::V(gi))
  } else NA

  return(list(remain = remain.percent, vul = vul, complexity = cpx))
}
