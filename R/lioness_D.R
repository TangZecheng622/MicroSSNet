lioness_D <- function(sel_otu_table, group_df, vscol1 = "group",vscol2 = vscol1) {
  # Load required package
  if (!requireNamespace("lionessR", quietly = TRUE)) {
    stop("Package 'lionessR' is required but not installed.")
  }

  # Create output directory
  output_dir <- file.path("SSN", "LIONESS-D")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Check and process group information
  group_df[[vscol1]] <- as.factor(group_df[[vscol1]])
  group_df[[vscol2]] <- as.factor(group_df[[vscol2]])
  cor_list <- list()

  for(i in unique(group_df[[vscol1]])){
    # i <- unique(group_df[[group]])[[1]]
    # print(paste("Group:",i))
    sel_i <- group_df$sample[group_df[[vscol1]] == i]
    sel_i <- sel_i[sel_i %in% colnames(sel_otu_table)]
    sel_otu_table_i <- sel_otu_table[, sel_i, drop = FALSE]
    for (j in unique(group_df[[vscol2]])) {
      # j <- unique(group_df[[vscol2]])[[2]]
      sel_j <- group_df$sample[group_df[[vscol2]] == j]
      sel_j <- sel_j[sel_j %in% colnames(sel_otu_table_i)]
      if(length(sel_j) == 0){
        next
      }else{
        # sel_j <- sel_j[sel_j %in% colnames(sel_otu_table_i)]
        sel_otu_table_j <- sel_otu_table_i[, sel_j, drop = FALSE]
        cor_ij <- lionessR::lioness(sel_otu_table_j)
        n <- paste(i,"_",j)
        cor_list[[n]] <- cor_ij
        cor_ij <- NULL
      }

    }

  }
  edge_keys <- cor_list[[1]][, 1:2]
  sample_matrices <- lapply(cor_list, function(df) df[, -(1:2)])
  names(sample_matrices) <- NULL  # 移除 list 名
  all_samples <- do.call(cbind, sample_matrices)
  combined_matrix <- cbind(edge_keys, all_samples)
  return(combined_matrix)
}


