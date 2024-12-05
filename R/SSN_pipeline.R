#' Single Sample Network Analysis Pipeline Function
#'
#' This function performs network construction and analysis for single-sample networks based on microbial abundance data.
#'
#' @param table1 Data frame containing species abundance data (rows as species, columns as samples).
#' @param group_df Optional data frame containing sample grouping information. Must include columns specified in `vscol1` and `vscol2`.
#' @param vscol1 Character string specifying the column name in `group_df` that indicates the primary group labels.
#' @param vscol2 Character string specifying the column name in `group_df` that indicates the secondary group labels.
#' @param ssn_method Character string specifying the method for single-sample network construction. Options are "SSN" or "Lioness". Default is "SSN".
#' @param control Character string specifying the control group.
#'        Can take one of the following values:
#'        - "all" : Use all samples for network construction (no control group).
#'        - "pergroup" : Perform single-sample PCC for each group separately.
#'        - A specific group name (e.g., "control", "CK") : Specify a group name within `group_df` to use as the control group for network construction.
#'        Default is "control".#' @param top Numeric value specifying the number of top edges to select. If `top > 0`, selects top `top` edges based on weight. Default is 0 (no filtering).
#' @param pvl_threshold Numeric value between 0 and 1 for prevalence threshold in filtering species. Default is 0.5.
#' @param r_threshold Numeric value for the correlation coefficient threshold. Edges with absolute correlation below this value will be removed. Default is 0.3.
#' @param log Logical value indicating whether to table1 to the log(table1 + 10e-16) . Default is TRUE.
#' @param scale Logical value indicating whether to scale network weights to the range . Default is TRUE.
#' @param save Logical value indicating whether to save the outputs to files. Default is TRUE.
#' @return The function outputs various files and plots in the `./ssn/` directory.
#' @importFrom dplyr full_join
#' @importFrom purrr map reduce
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom stats cor
#' @importFrom methods is
#' @export
#'
#' @examples
#' # Example usage:
#' otu_data <- data.frame(
#' OTU_ID = paste0("OTU_", 1:10),
#' CK1 = c(10, 20, 30, 30, 5, 6, 7, 8, 12, 10),
#' CK2 = c(12, 24, 28, 28, 6, 7, 8, 9, 11, 9),
#' CK3 = c(11, 22, 29, 27, 7, 5, 7, 8, 13, 11),
#' CK4 = c(10, 21, 30, 30, 6, 6, 8, 8, 12, 10),
#' CK5 = c(11, 23, 31, 31, 5, 6, 7, 7, 13, 12),
#' CK6 = c(12, 25, 28, 29, 6, 5, 6, 9, 11, 11),
#' PD1 = c(10, 20, 30, 30, 5, 6, 7, 8, 12, 10),
#' PD2 = c(12, 24, 28, 28, 6, 7, 8, 9, 11, 9),
#' PD3 = c(11, 22, 29, 27, 7, 5, 7, 8, 13, 11),
#' PD4 = c(10, 21, 30, 30, 6, 6, 8, 8, 12, 10),
#' PD5 = c(11, 23, 31, 31, 5, 6, 7, 7, 13, 12),
#' PD6 = c(12, 25, 28, 29, 6, 5, 6, 9, 11, 11),
#' PDD1 = c(15, 25, 18, 16, 12, 14, 16, 10, 15, 18),
#' PDD2 = c(18, 30, 20, 19, 15, 17, 19, 12, 16, 20),
#' PDD3 = c(17, 28, 19, 18, 14, 16, 18, 11, 17, 19),
#' PDD4 = c(16, 27, 20, 17, 13, 15, 17, 13, 18, 17),
#' PDD5 =  c(19, 32, 22, 21, 16, 18, 20, 14, 17, 21),
#' PDD6 =  c(18, 31, 21, 20, 15, 17, 19, 12, 19, 20)
#' )
#' rownames(otu_data) <- otu_data$OTU_ID
#' otu_data <- otu_data[, -1]
#'
#' ssn_pipeline(
#'   table1 = otu_data,
#'   group_df = NULL,
#'   vscol1 = "group",
#'   vscol2 = "group",
#'   ssn_method = "SSN",
#'   control = "CK",
#'   top = NULL,
#'   pvl_threshold = 0,
#'   r_threshold = 0,
#'   scale = TRUE,
#'   save = TRUE
#' )


ssn_pipeline <- function(
    table1,
    group_df = NULL,
    vscol1,
    vscol2 = vscol1,
    ssn_method = "ssPCC",
    control = "pergroup",
    log = TRUE,
    top = NULL,
    pvl_threshold = 0.5,
    r_threshold = 0.3,
    scale = TRUE,
    save = TRUE,
    pca = TRUE,
    pcoa = TRUE,
    limma = TRUE,
    property = TRUE
) {

  # Parameter validation
  if (!is.data.frame(table1)) stop("Error: 'table1' must be a data frame.")
  if (!is.null(group_df) && !is.data.frame(group_df)) stop("Error: 'group_df' must be a data frame or NULL.")
  if (!is.character(vscol1)) stop("Error: 'vscol1' must be a character string.")
  if (!is.character(vscol2)) stop("Error: 'vscol2' must be a character string.")
  if (!ssn_method %in% c("ssPCC", "Lioness")) stop("Error: 'ssn_method' must be either 'ssPCC' or 'Lioness'.")
  if (!is.character(control) || !(control %in% c("all", "pergroup"))) {
    stop("Error: 'control' must be one of the following: 'all', 'pergroup', or a character string specifying a valid control group name (e.g., 'control').")
  }
  if (!is.null(top) && (!is.numeric(top) || top <= 0)) stop("Error: 'top' must be a positive number or NULL.")
  if (!is.numeric(pvl_threshold) || pvl_threshold < 0 || pvl_threshold > 1) stop("Error: 'pvl_threshold' must be between 0 and 1.")
  if (!is.numeric(r_threshold)) stop("Error: 'r_threshold' must be numeric.")
  if (!is.logical(scale)) stop("Error: 'scale' must be a logical value.")
  if (!is.logical(log)) stop("Error: 'log' must be a logical value.")
  if (!is.logical(save)) stop("Error: 'save' must be a logical value.")
  if (!is.logical(pca)) stop("Error: 'pca' must be a logical value.")
  if (!is.logical(pcoa)) stop("Error: 'pcoa' must be a logical value.")
  if (!is.logical(limma)) stop("Error: 'limma' must be a logical value.")
  if (!is.logical(property)) stop("Error: 'property' must be a logical value.")

  # Process group and table data
  check_result <- process_group_and_table(table = table1, group_df = group_df, vscol = vscol1)
  table1 <- check_result$table
  group_df <- check_result$group_df
  group_list <- check_result$group_list
  group_list2 <- unique(group_df[[vscol2]])
  table_i_list <- list()
  for (i in group_list) {
    i <- as.character(i)
    # i = group_list[[1]]

    # 根据group_df过滤table1中的列
    selected_table <- table1[, group_df[,1][group_df[[vscol1]] == i], drop = FALSE]

    for (j in group_list2) {
      # j <- group_list2[[2]]
      j <- as.character(j)
      cols_to_select <- group_df[, 1][group_df[[vscol2]] == j]
      valid_cols <- cols_to_select[cols_to_select %in% colnames(selected_table)]
      if(length(valid_cols) == 0){
        next
      }else{
        table_i <- selected_table[, valid_cols]
        table_i <- filter_OTU2(table_i, Pre = pvl_threshold)
        n <- paste(i,"_",j)
        table_i_list[[n]] <- table_i
      }
    }

    merged_table <- table_i_list %>%
      purrr::map(~ .x %>% tibble::rownames_to_column(var = "row_name")) %>% # 将行名转换为一列
      purrr::reduce(dplyr::full_join, by = "row_name") %>%                  # 按"row_name"列合并
      replace(is.na(.), 0) %>%                                              # 填充缺失值为0
      tibble::column_to_rownames(var = "row_name")                           # 将"row_name"列转换回行名
  }
  table1 <- table1[rownames(merged_table),]

  if(log == TRUE){
    table1 <- log(t(table1)+10e-16)
  }

  # Create directory for output
  if (!dir.exists("SSN")) dir.create("SSN")

  # Single Sample Network Construction
  if (ssn_method == "ssPCC") {

    if (control == "all") {
      sspcc_cal(sel_otu_table = table1)
    }else if (control == "pergroup"){
      sspcc_cal3(sel_otu_table = table1,group_df,group = vscol1)
    }else{
      sspcc_cal2(sel_otu_table = table1, group_df = group_df, ck = control, group = vscol1)
      group_df[[vscol2]][group_df[[vscol1]] == control] <- "control"
    }
    n_genes <- nrow(table1)
    cor <- make_one_edgelist_file_general(
      dir = "./SSN/ssPCC/",
      pattern = "ssPCC.*\\.tsv",
      file_name = "SSN-Total",
      col_n = 4,
      p_values = TRUE,
      col_p_n = 3,
      n_genes = n_genes
    )
  }

  if (ssn_method == "Lioness") {
    netFun <- function(x) {
      cor(x, method = "pearson")
    }
    cor <- lionessR::lioness(table1, netFun)
  }

  cor[is.na(cor)] <- 0
  cor <- cor[rowSums(abs(cor[, -c(1, 2)]) > 0) > 0, ]

  # Scale network weights
  if (scale) {
    cor <- scale_network(network = cor, vose = save)
  }

  # Apply correlation threshold
  if (!is.null(r_threshold)) {
    cor <- as.data.frame(cor)
    cor[, -c(1, 2)][abs(cor[, -c(1, 2)]) < r_threshold] <- 0
  }

  # Filter top edges
  if (top > 0 && !is.null(top)) {
    cor_top <- select_top_edges(n = top, table = cor, group = "TotalSSN", remove_zero_rows = TRUE, vose = save)
  } else {
    cor_top <- cor[rowSums(abs(cor[, -c(1, 2)])) > 0, ]
  }

  # Calculate Sum of Weights (SOW) for nodes
  SOW <- Calculate_sum_of_weights(network = cor_top)

  # PCA and PCoA analysis

  if(pca == TRUE){
    PCA_draw(SOW, group_df, vscol2, save = save, showType = "centroid", mode = "PCA", offset = FALSE)
    PCA_draw(cor_top, group_df, vscol2, save = save, showType = "centroid", mode = "PCA", offset = TRUE)
  }
  if(pcoa == TRUE){
    PCA_draw(SOW, group_df, vscol2, save = save, showType = "centroid", mode = "PCOA", offset = FALSE)
    PCA_draw(cor_top, group_df, vscol2, save = save, showType = "centroid", mode = "PCOA", offset = TRUE)
  }


  # Differential analysis using limma
  if(limma == TRUE){
    Differential_nodes(SOW, group_df, offset = FALSE, log_transform = FALSE, vscol = vscol2, save = save, plot = TRUE)
    Differential_nodes(cor_top, group_df, offset = TRUE, log_transform = FALSE, vscol = vscol2, save = save, plot = FALSE)
  }


  # Network properties analysis
  if(property == TRUE){
    determineCharacteristics(cor_top)
  }

}
