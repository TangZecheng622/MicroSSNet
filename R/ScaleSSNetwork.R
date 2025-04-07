#' Scale Network Weights
#'
#' This function scales the weights of a network (correlation matrix) to be within the range .
#' It optionally saves the scaled network to a file.
#'
#' @param network Data frame. The input network with columns for edges and weight values.
#' @param offset Logical. If TRUE, removes the first column (for cases with unnecessary index column).
#' @param vose Logical. If TRUE, saves the scaled network to a specified output file.
#' @param group Character. Used as the identifier in the output file name when `vose = TRUE`.
#' @return Data frame. Scaled network with weights normalized to the range.
#' @importFrom data.table fwrite
#' @export
scale_network <- function(network, offset = FALSE, vose = FALSE, group = "Total") {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Handle potential index column removal
  if (offset) {
    network <- network[, -1]
  }

  # Separate edge columns and weight columns
  edges <- network[, 1:2]
  weights <- network[, 3:ncol(network)]
  weights <- as.data.frame(lapply(weights, as.numeric))


  # Scale weights to [-1, 1] using the maximum absolute value
  max_weight <- max(abs(weights), na.rm = TRUE)
  if (max_weight == 0) {
    stop("Error: Maximum weight is zero, scaling cannot be performed.")
  }

  weights <- weights / max_weight

  # Combine edges and scaled weights
  network_scaled <- cbind(edges, weights)

  # Optionally save the scaled network to a file
  if (vose) {
    file_name <- paste0("./SSN/ssPCC/Scaled_", group, "_-1_1scaled.tsv")
    dir_path <- dirname(file_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    data.table::fwrite(network_scaled, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(network_scaled)
}
