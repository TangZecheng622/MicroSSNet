#' Determine Network Characteristics
#'
#' This function calculates local and global properties for each sample in a network.
#' It extracts subgraphs with non-zero weights, computes node-level properties, and saves results.
#'
#' @param network Data frame representing the network, where the first two columns are edges and subsequent columns are sample weights.
#' @return Saves local and global network properties as .tsv files and returns a data frame of global properties.
#' @importFrom igraph graph.data.frame E
#' @export



determineCharacteristics <- function(network) {


  output_dir <- './SSN/Properties/'
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  global_list <- list()
  n <- ncol(network) - 2

  for (i in seq_len(n)) {

    sample <- colnames(network)[i + 2]


    non_null_ind <- which(network[[i + 2]] != 0)


    sample_netw <- network[non_null_ind, c(1, 2, i + 2), drop = FALSE]
    colnames(sample_netw) <- c("From", "To", "Weight")


    if (nrow(sample_netw) > 0) {
      single_sample_graph <- igraph::graph.data.frame(sample_netw, directed = FALSE)
      igraph::E(single_sample_graph)$correlation <- sample_netw$Weight
      igraph::E(single_sample_graph)$weight <- abs(sample_netw$Weight)


      ssn_p <- node_properties2(igraph = single_sample_graph, zipi = TRUE, tag = sample, RM = FALSE, output = file.path(output_dir, "ZiPi/"))


      local_pro <- ssn_p[[1]]
      write.table(local_pro, file = file.path(output_dir, paste0(sample, '_SSN_Properties.tsv')), sep = "\t", quote = FALSE, row.names = FALSE)


      global_pro <- ssn_p[[2]]
      global_list[[i]] <- global_pro
    } else {
      message("sample: ", sample, " skipping")
    }
  }


  combined_df <- do.call(rbind, global_list)
  write.table(combined_df, file = file.path(output_dir, 'SSN_Global_properties.tsv'), sep = "\t", quote = FALSE, row.names = FALSE)

  return(combined_df)
}
