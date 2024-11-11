#' Filter OTUs Based on Average Abundance
#'
#' This function filters OTUs (Operational Taxonomic Units) from the abundance table based on average abundance.
#'
#' @param otu_table Data frame where rows are OTUs and columns are samples.
#' @param Top Numeric value. If `Top >= 1`, selects the top `Top` OTUs based on average abundance.
#' If `0 < Top < 1`, selects OTUs with average abundance greater than `Top`.
#' @return A filtered OTU table.
#' @importFrom dplyr arrange desc
#' @importFrom utils head
#' @export
filter_OTU <- function(otu_table, Top = NULL) {
  if (!is.data.frame(otu_table)) {
    stop("Error: 'otu_table' must be a data frame.")
  }

  # Check for samples with zero total counts
  sample_sums <- colSums(otu_table)
  if (any(sample_sums == 0)) {
    stop("Error: One or more samples have zero total counts.")
  }

  if (!is.null(Top) && Top != 0) {
    relative_abundance <- sweep(otu_table, 2, sample_sums, FUN = "/")
    relative_abundance$mean <- rowMeans(relative_abundance, na.rm = TRUE)
    relative_abundance$ID <- rownames(relative_abundance)

    relative_abundance <- dplyr::arrange(relative_abundance, dplyr::desc(mean))

    if (Top > 0 && Top < 1) {
      subtab <- relative_abundance[relative_abundance$mean > Top, ]
    } else if (Top >= 1) {
      subtab <- utils::head(relative_abundance, Top)
    } else {
      stop("Error: 'Top' must be a positive number.")
    }

    otu_table_filtered <- otu_table[subtab$ID, , drop = FALSE]
    return(otu_table_filtered)
  } else {
    return(otu_table)
  }
}
