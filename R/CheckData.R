#' Process Group and Table Data
#'
#' A helper function that aligns a data table with a grouping data frame (`group_df`).
#' The function verifies and harmonizes the sample columns and grouping columns in
#' both `table` and `group_df`. If `group_df` is not provided, it auto-generates groups
#' based on the column names of the table.
#'
#' @param table Data frame or matrix containing the data to be grouped.
#' @param group_df Data frame containing sample and group information (optional).
#' @param vscol Character string specifying the group column in `group_df`.
#' @return A list with three elements:
#'   \item{table}{Modified data table aligned with `group_df` samples.}
#'   \item{group_df}{Processed grouping data frame with consistent sample and group columns.}
#'   \item{group_list}{Unique groups extracted from `group_df`.}
#' @export
process_group_and_table <- function(table, group_df = NULL, vscol) {
  # group_list <- NULL
  # Check if group_df is provided
  if (!is.null(group_df)) {
    # Detect column names in group_df
    strings <- colnames(group_df)
    sample_match <- grepl("sample", strings, ignore.case = TRUE)
    group_match <- grepl(vscol, strings, ignore.case = TRUE)

    # Check for 'sample' column
    if (any(sample_match)) {
      sample_col <- strings[sample_match][1]  # First matching column as sample column
    } else {
      stop("Error: 'sample' column not found in group_df")
    }

    # Check for 'group' column matching vscol
    if (any(group_match)) {
      group_col <- strings[group_match][1]  # First matching column as group column
    } else {
      stop(paste("Error:", vscol, "column not found in group_df"))
    }

    # Set standardized columns for 'group' and 'sample' in group_df
    group_df$group <- group_df[[group_col]]
    group_df$sample <- group_df[[sample_col]]
    group_list <- unique(group_df$group)
  } else {
    # If no group_df, create sample and group data based on column names of table
    samples <- colnames(table)
    group_df <- data.frame(
      sample = samples,
      group = as.factor(gsub("[[:digit:]]|[[:punct:]]", "", samples))
    )
    group_list <- unique(group_df$group)
  }

  # Verify alignment between 'table' columns and 'group_df' samples
  if (!setequal(colnames(table), group_df$sample)) {
    message("Warning: table columns and group_df samples do not match. Adjusting data to match common samples.")

    # Find common samples between table and group_df
    common_samples <- intersect(colnames(table), group_df$sample)

    # Check if there are common samples
    if (length(common_samples) == 0) {
      stop("Error: No common samples found between 'table' and 'group_df'.")
    }

    # Filter both table and group_df to include only common samples
    table <- table[, common_samples, drop = FALSE]
    group_df <- group_df[group_df$sample %in% common_samples, ]
    group_list <- unique(group_df$group)
  }

  # Return modified table, group_df, and group_list as a list
  list(table = table, group_df = group_df, group_list = group_list)
}
