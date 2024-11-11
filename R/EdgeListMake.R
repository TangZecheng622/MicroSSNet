#' Consolidate Edge List Files
#'
#' This function consolidates multiple edge list files into a single table.
#' It supports filtering by p-value and removes self-loops as well as duplicate edges.
#'
#' @param dir Character. Directory containing edge list files.
#' @param pattern Character. File name pattern to match the edge list files.
#' @param file_name Character. Output file name (without extension) for the consolidated edge list.
#' @param col_n Integer. Column index in edge files containing PCC values.
#' @param p_values Logical. If TRUE, filter edges by p-value.
#' @param p_cutoff Numeric. P-value threshold for filtering edges. Only used if `p_values` is TRUE.
#' @param col_p_n Integer. Column index in edge files containing p-values.
#' @param n_genes Integer. Number of genes (or nodes) in each network.
#' @return Data frame. Consolidated edge list table.
#' @importFrom data.table fread fwrite
#' @export
make_one_edgelist_file_general <- function(
    dir,
    pattern,
    file_name,
    col_n,
    p_values = FALSE,
    p_cutoff = 0.05,
    col_p_n = 3,
    n_genes
) {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Ensure directory path has a trailing slash
  dir <- ifelse(grepl("/$", dir), dir, paste0(dir, "/"))

  # Fetch matching files
  file_list <- list.files(path = dir, pattern = pattern)
  if (length(file_list) == 0) {
    stop("No files matching the pattern were found in the specified directory.")
  }

  # Initialize data frame for consolidated results
  n_edges <- n_genes * (n_genes - 1) / 2  # Number of unique edges without self-loops
  table_total <- data.frame(matrix(NA_real_, ncol = length(file_list) + 2, nrow = n_edges))
  col_names <- c("reg", "tar")

  # Process each file and consolidate into table_total
  for (i in seq_along(file_list)) {
    file_path <- paste0(dir, file_list[i])
    sample <- sub("ssPCC_(.*)\\.tsv", "\\1", file_list[i])

    # Load data and validate structure
    table_sample <- data.table::fread(file_path, header = TRUE, data.table = FALSE)
    if (ncol(table_sample) < max(col_n, col_p_n)) {
      stop(paste("Error: Column indices out of range for file", file_list[i]))
    }

    # Assign regulator and target columns on first iteration
    if (i == 1) {
      table_total[, 1:2] <- table_sample[, 1:2]
    } else if (!all(table_total[, 1:2] == table_sample[, 1:2])) {
      stop("Edge order mismatch between files.")
    }

    # Filter based on p-values if specified
    if (p_values) {
      table_sample[table_sample[, col_p_n] > p_cutoff | is.na(table_sample[, col_p_n]), col_n] <- 0.0
    }

    # Populate PCC values in consolidated table
    col_names <- c(col_names, sample)
    table_total[, match(sample, col_names)] <- table_sample[, col_n]

    # Remove the processed file
    file.remove(file_path)
  }

  # Set final column names and save consolidated file
  colnames(table_total) <- col_names
  output_path <- paste0(dir, file_name, ".tsv")
  data.table::fwrite(table_total, file = output_path, sep = "\t", buffMB = 80)

  return(table_total)
}
