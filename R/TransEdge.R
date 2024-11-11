#' Get Full Edge List from Matrix
#'
#' This function converts a matrix into an edge list data frame.
#'
#' @param mat Numeric matrix.
#' @return A data frame containing edge pairs and weights.
#' @export
get_full_edge_list <- function(mat) {
  all_pairs <- which(!is.na(mat), arr.ind = TRUE)

  node1 <- rownames(mat)[all_pairs[, 1]]
  node2 <- colnames(mat)[all_pairs[, 2]]

  weights <- mat[all_pairs]

  edge_list <- data.frame(Node1 = node1, Node2 = node2, Weight = weights, stringsAsFactors = FALSE)

  return(edge_list)
}
