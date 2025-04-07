#' Filter OTUs Based on Prevalence
#'
#' This function filters OTUs based on their prevalence across samples.
#'
#' @param otu_table Data frame where rows are OTUs and columns are samples.
#' @param Pre Numeric value between 0 and 1. Filters OTUs that are present in at least `Pre * 100%` of samples.
#' @return A filtered OTU table.
#' @export
filter_OTU2 <- function(otu_table, Pre = NULL) {
  if (!is.data.frame(otu_table)) {
    stop("Error: 'otu_table' must be a data frame.")
  }
  if (!is.null(Pre) && (Pre < 0 || Pre > 1)) {
    stop("Error: 'Pre' must be between 0 and 1.")
  }

  if (!is.null(Pre)) {
    prevalence <- rowSums(otu_table != 0)
    otu_table_filtered <- otu_table[prevalence >= round(ncol(otu_table) * Pre), ]
    return(otu_table_filtered)
  } else {
    return(otu_table)
  }
}
