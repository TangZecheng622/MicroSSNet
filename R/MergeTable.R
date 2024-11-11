#' Merge and Filter OTU Tables from Two Domains
#'
#' This function processes and merges two OTU tables with filtering options.
#'
#' @param table1 Data frame of the first OTU table.
#' @param table2 Data frame of the second OTU table.
#' @param table1name Character string specifying the prefix name for species in `table1`. Default is "domain1".
#' @param table2name Character string specifying the prefix name for species in `table2`. Default is "domain2".
#' @param top Numeric value. If `top >= 1`, selects the top `top` species based on average abundance.
#' If `0 < top < 1`, selects species with average abundance greater than `top`.
#' @param pvl.threshold Numeric value between 0 and 1 for prevalence threshold in filtering species.
#' @return Data frame of merged OTU tables with filtered data.
#' @export
merge_bio <- function(
    table1,
    table2,
    table1name = "domain1",
    table2name = "domain2",
    top = NULL,
    pvl.threshold = 0.5
) {
  # Helper function to process and filter OTU data
  process_otu <- function(data, label, top, threshold) {
    if (!is.null(data)) {
      data_filtered <- data
      if (!is.null(pvl.threshold)) {
        data_filtered <- filter_OTU2(data_filtered, Pre = threshold)
      }
      if (!is.null(top)) {
        data_filtered <- filter_OTU(data_filtered, Top = top)
      }
      
      otu_table <- as.data.frame(data_filtered)
      rownames(otu_table) <- paste(label, rownames(otu_table), sep = "_")
      otu_table$filed <- rep(label, nrow(otu_table))
      
      return(otu_table)
    } else {
      return(NULL)
    }
  }
  
  # Parameter validation
  if (!is.data.frame(table1) || !is.data.frame(table2)) {
    stop("Error: Both 'table1' and 'table2' must be data frames.")
  }
  if (!is.character(table1name) || !is.character(table2name)) {
    stop("Error: 'table1name' and 'table2name' must be character strings.")
  }
  if (!is.null(top) && (!is.numeric(top) || top <= 0)) {
    stop("Error: 'top' must be a positive number or NULL.")
  }
  if (!is.numeric(pvl.threshold) || pvl.threshold < 0 || pvl.threshold > 1) {
    stop("Error: 'pvl.threshold' must be a numeric value between 0 and 1.")
  }
  
  # Process each OTU table
  otu_data1 <- process_otu(table1, table1name, top, pvl.threshold)
  otu_data2 <- process_otu(table2, table2name, top, pvl.threshold)
  
  # Combine processed OTU tables
  merged_otu_table <- do.call(rbind, Filter(Negate(is.null), list(otu_data1, otu_data2)))
  
  return(merged_otu_table)
}
