#' Calculate Sum of Weights for Network Nodes
#'
#' Calculates the sum of edge weights (absolute or non-absolute) for each node in a network across multiple samples.
#'
#' @param network Data frame. Input network with columns for nodes and weight values.
#' @param offset Logical. If TRUE, removes the first column if it’s an index.
#' @param log Logical. If TRUE, applies a log transformation to the result.
#' @param abs Logical. If TRUE, calculates the sum of absolute weights.
#' @return Data frame. Sum of weights for each node across samples.
#' @export
Calculate_sum_of_weights <- function(network, offset = FALSE, log = FALSE, abs = TRUE) {
  # Remove the index column if specified
  if (offset) {
    network <- network[, -1]
  }
  
  # Extract unique nodes and initialize a data frame for results
  genes <- unique(c(as.character(network[[1]]), as.character(network[[2]])))
  SOW <- matrix(0, nrow = length(genes), ncol = ncol(network) - 2)
  rownames(SOW) <- genes
  colnames(SOW) <- colnames(network)[-c(1, 2)]
  
  # Loop over genes to calculate sum of edge weights
  for (gene in genes) {
    gene_indices <- which(network[[1]] == gene | network[[2]] == gene)
    
    # Sum weights for this gene across all samples
    if (abs) {
      SOW[gene, ] <- colSums(abs(network[gene_indices, -c(1, 2)]))
    } else {
      SOW[gene, ] <- colSums(network[gene_indices, -c(1, 2)])
    }
  }
  
  # Log transformation if specified
  if (log) {
    SOW <- log1p(SOW)
  }
  
  return(as.data.frame(SOW))
}
