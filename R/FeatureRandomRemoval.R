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
    calculate_cpx = TRUE
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
    tag = tag
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
                    calculate_vul = TRUE, calculate_rob = TRUE, calculate_cpx = TRUE, tag = "default") {
  # Parameter checks
  if (!is.matrix(netRaw)) stop("Error: 'netRaw' must be a matrix.")
  if (!is.numeric(rm.p.list) || any(rm.p.list < 0) || any(rm.p.list > 1)) stop("Error: 'rm.p.list' values must be between 0 and 1.")
  if (!is.numeric(sp.ra)) stop("Error: 'sp.ra' must be a numeric vector.")
  if (!is.numeric(nperm) || nperm <= 0) stop("Error: 'nperm' must be a positive integer.")
  if (!is.numeric(ncpus) || ncpus <= 0) stop("Error: 'ncpus' must be a positive integer.")

  # Initialize cluster
  cl <- parallel::makeCluster(ncpus)
  on.exit(parallel::stopCluster(cl))

  dat1_results <- data.frame()
  dat2_results <- data.frame()

  for (rm.percent in rm.p.list) {
    # Export necessary variables
    parallel::clusterExport(cl, c("netRaw", "boot_feature", "Vulnerability", "rm.percent",
                                  "abundance.weighted", "sp.ra", "calculate_vul", "calculate_rob", "calculate_cpx"),
                            envir = environment())
    # Run simulations in parallel
    results_list <- parallel::parLapply(cl, 1:nperm, function(i) {
      boot_feature(i, netRaw, rm.percent, sp.ra, abundance.weighted, calculate_vul, calculate_rob, calculate_cpx)
    })

    # Extract results
    Proportion.remain <- sapply(results_list, `[[`, "remain")
    vulnerability <- sapply(results_list, `[[`, "vul")
    complexity <- sapply(results_list, `[[`, "complexity")

    # Compile summary results
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

    # Compile detailed results
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
