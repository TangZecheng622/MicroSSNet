#' Select Top Edges from Network
#'
#' Selects the top `n` edges from a network based on edge weights for each sample.
#' Edges not in the top `n` are set to zero.
#'
#' @param n Integer. Number of top edges to select.
#' @param table Data frame. Input network with columns for edges and weight values.
#' @param group Character. Identifier used in the output file name when `vose = TRUE`.
#' @param remove_zero_rows Logical. If TRUE, removes rows where all edge weights are zero across samples.
#' @param offset Logical. If TRUE, removes the first column if it's an index.
#' @param vose Logical. If TRUE, saves the filtered network to a file.
#' @return Data frame. Network with only the top `n` edges retained for each sample.
#' @importFrom data.table fwrite
#' @export
select_top_edges <- function(n, table, ssn_model = ssn_model,ssn_dir = ssn_dir, remove_zero_rows = FALSE, offset = FALSE, vose = TRUE) {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Handle potential index column removal
  if (offset) {
    table <- table[, -1]
  }

  # Select top edges for each sample
  for (i in 3:ncol(table)) {
    idx <- order(abs(table[[i]]), decreasing = TRUE)
    threshold_idx <- idx[(n + 1):length(idx)]
    table[threshold_idx, i] <- 0
  }

  # Remove rows with all zeros if specified
  if (remove_zero_rows) {
    table <- table[rowSums(abs(table[, -c(1, 2)])) > 0, ]
  }

  # Optionally save the filtered network to a file
  if (vose) {
    file_name <- paste0(ssn_dir, ssn_model, "_top_", n, "_edges.tsv")
    dir_path <- dirname(file_name)
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
    }
    data.table::fwrite(table, file = file_name, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  return(table)
}
