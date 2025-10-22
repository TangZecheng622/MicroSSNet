#' Single-Sample Network Construction using PCC
#'
#' Constructs single-sample networks by calculating pairwise PCC (Pearson Correlation Coefficient)
#' differences between each sample and the entire dataset.
#'
#' @param sel_otu_table Data frame. The abundance table of selected OTUs, where rows are OTUs and columns are samples.
#' @return None. Saves each sample's single-sample PCC to individual TSV files in the `./ssn/sspcc/` directory.
#' @importFrom data.table fwrite
#' @importFrom stats pnorm
#' @export
sspcc_cal <- function(sel_otu_table,base_correlation = NULL,
                      r_emergent_threshold = 0.6) {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Create output directory
  output_dir <- file.path("SSN", "ssPCC")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Process input data

  sel_samples <- colnames(sel_otu_table)
  total_n <- length(sel_samples)
  if (total_n < 10) stop("Need at least 10 samples to compute PCC.")
  sel_otu_table_lg <- sel_otu_table


  if (!is.null(base_correlation)){
    r_new_all = base_correlation
  }else{
    cat("Base correlation is NULL, computing base correlation...\n")
    r_new_all <- corMicro(sel_otu_table_lg, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = "ssn", vose = FALSE)[[1]]
  }

  # Filter base correlation matrix
  r_new_all[is.na(r_new_all)] <- 0
  r_new_all <- pmin(pmax(r_new_all, -(1-1e-10)), (1-1e-10))
  upper_tri <- upper.tri(r_new_all, diag = FALSE)

  # Construct single-sample PCC for each sample
  for (num in seq_len(ncol(sel_otu_table))) {
    sampleID <- colnames(sel_otu_table)[num]
    cat("Processing Sample:", sampleID, "\n")

    sample_table <- sel_otu_table_lg[, -num, drop = FALSE]

    r_bg_mat <- corMicro(sample_table, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = sampleID, vose = FALSE)[[1]]
    r_bg_mat[is.na(r_bg_mat)] <- 0

    # 限制 PCC 值范围

    r_bg_mat <- pmin(pmax(r_bg_mat, -(1-1e-10)), (1-1e-10))

    # Compute differential PCC and significance
    r_bg_mat <- r_bg_mat[rownames(r_new_all), colnames(r_new_all)]

    delta_PCC <- r_new_all - r_bg_mat

    Zvals <- delta_PCC / (1 - (r_bg_mat^2)) * (total_num - 2)

    node_1 <- rownames(delta_PCC)[row(delta_PCC)[upper_tri]]
    node_2 <- colnames(delta_PCC)[col(delta_PCC)[upper_tri]]
    r_bg   <- r_bg_mat[upper_tri]
    r_new  <- r_new_all[upper_tri]
    dlt    <- delta_PCC[upper_tri]
    Zv     <- Zvals[upper_tri]

    # p 与 FDR（仅上三角）
    p_val <- 2 * stats::pnorm(-abs(Zv))
    q_val <- p.adjust(p_val, method = "BH")

    # 方向与标签
    change    <- ifelse(abs(r_new) > abs(r_bg), "gained", "lost")
    sign_flip <- sign(r_bg) * sign(r_new) == -1
    emergent  <- abs(r_bg) < r_emergent_threshold & abs(r_new) >= r_emergent_threshold
    emergent[is.na(emergent)] <- FALSE

    # 汇总与保存
    res <- data.frame(
      node_1    = node_1,
      node_2    = node_2,
      r_bg      = as.numeric(r_bg),         # 参考（去掉该样本）
      r_new     = as.numeric(r_new),        # 加入该样本（全体）
      Delta_PCC = as.numeric(dlt),
      Z         = as.numeric(Zv),
      p         = as.numeric(p_val),
      q_BH      = as.numeric(q_val),
      change    = change,       # gained / lost（以 |r| 比较）
      sign_flip = sign_flip,    # TRUE / FALSE
      emergent  = emergent,     # TRUE / FALSE
      stringsAsFactors = FALSE
    )
    colnames(res)[5] <- paste(sampleID,"PCC_value",  sep = "_")

    # Save results
    file_name <- file.path(output_dir, paste0("ssPCC_", sampleID, ".tsv"))
    data.table::fwrite(res, file = file_name, sep = "\t")
  }
}

#' Comparison-based Single-Sample Network Construction
#'
#' Constructs single-sample networks by calculating PCC differences between each sample in the specified group
#' and a control group.
#'
#' @param sel_otu_table Data frame. The abundance table of selected OTUs.
#' @param group_df Data frame. Contains group and sample information.
#' @param ck Character. The name of the control group in `group_df`.
#' @param group Character. The column name in `group_df` specifying the group variable.
#' @return None. Saves each sample's single-sample PCC to individual TSV files in the `./ssn/sspcc/` directory.
#' @importFrom data.table fwrite
#' @importFrom stats pnorm
#' @export
sspcc_cal2 <- function(sel_otu_table, group_df, ck = "CK", group = "group",r_emergent_threshold = 0.6) {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Create output directory
  output_dir <- file.path("SSN", "ssPCC")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Check and process group information
  group_df[[group]] <- as.factor(group_df[[group]])
  sel_ck <- group_df$sample[group_df[[group]] == ck]
  sel_ck <- sel_ck[sel_ck %in% colnames(sel_otu_table)]
  ck_otu_table <- sel_otu_table[, sel_ck, drop = FALSE]
  ck_otu_table_lg <- ck_otu_table

  # Compute base PCC matrix for control group
  ck_matrix <- corMicro(ck_otu_table_lg, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = ck, vose = FALSE)[[1]]
  sel_otu_table_lg <- sel_otu_table[, !colnames(sel_otu_table) %in% sel_ck, drop = FALSE]
  # sel_otu_table_lg <- sel_otu_table
  total_num <- ncol(ck_otu_table)
  if (total_num < 10) stop("Too few control samples (n<10) for PCC.")
  sel_samples <- colnames(sel_otu_table)
  # sspcc_cal(ck_otu_table,base_correlation = ck_matrix)
  ck_matrix <- pmin(pmax(ck_matrix, -(1-1e-10)), (1-1e-10))
  upper_tri <- upper.tri(ck_matrix, diag = FALSE)

  for (num in seq_len(ncol(sel_otu_table_lg))) {
    # num=7
    sampleID <- colnames(sel_otu_table_lg)[num]
    cat("Processing Sample:", sampleID, "\n")

    # Construct sample matrix including control group
    sample_table <- cbind(ck_otu_table, sel_otu_table_lg[, num, drop = FALSE])
    sample_correlation <- corMicro(sample_table, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = sampleID, vose = FALSE)[[1]]

    # Compute differential PCC and significance
    sample_correlation <- sample_correlation[rownames(ck_matrix), rownames(ck_matrix)]
    sample_correlation[is.na(sample_correlation)] <- 0
    # 限制 PCC 值范围
    sample_correlation <- pmin(pmax(sample_correlation, -(1-1e-10)), (1-1e-10))

    delta_PCC <- sample_correlation - ck_matrix

    Zvals <- delta_PCC / ((1 - (ck_matrix^2)) / (total_num - 1))

    node_1 <- rownames(delta_PCC)[row(delta_PCC)[upper_tri]]
    node_2 <- colnames(delta_PCC)[col(delta_PCC)[upper_tri]]
    r_bg   <- ck_matrix[upper_tri]           # 基线 r
    r_new  <- sample_correlation[upper_tri]  # 加入后的 r'
    dlt    <- delta_PCC[upper_tri]           # ΔPCC
    Zv     <- Zvals[upper_tri]               # Z

    p_val  <- 2 * stats::pnorm(-abs(Zv))
    q_val  <- p.adjust(p_val, method = "BH")

    # 加强/减弱/翻转标注（以 |r'| vs |r| 为准）
    gained_lost <- ifelse(abs(r_new) > abs(r_bg), "gained", "lost")
    sign_flip   <- sign(r_bg) * sign(r_new) == -1
    # emergent（参考接近0时的新生关联）
    emergent    <- abs(r_bg) < r_emergent_threshold & abs(r_new) >= r_emergent_threshold
    emergent[is.na(emergent)] <- FALSE

    # 汇总结果
    res <- data.frame(
      node_1 = node_1,
      node_2 = node_2,
      r_bg      = as.numeric(r_bg),         # 参考（去掉该样本）
      r_new     = as.numeric(r_new),        # 加入该样本（全体）
      Delta_PCC = as.numeric(dlt),
      Z         = as.numeric(Zv),
      p         = as.numeric(p_val),
      q_BH      = as.numeric(q_val),
      change  = gained_lost,   # gained / lost
      sign_flip = sign_flip,   # TRUE / FALSE
      emergent  = emergent,    # TRUE / FALSE
      stringsAsFactors = FALSE
    )

    colnames(res)[5] <- paste(sampleID,"PCC_value",  sep = "_")

    # Save results
    file_name <- file.path(output_dir, paste0("ssPCC_", sampleID, ".tsv"))
    data.table::fwrite(res, file = file_name, sep = "\t")
  }
}

sspcc_cal3 <- function(sel_otu_table, group_df, group = "group",base_correlation = NULL) {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Create output directory
  output_dir <- file.path("SSN", "ssPCC")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Check and process group information
  group_df[[group]] <- as.factor(group_df[[group]])
  for(i in unique(group_df[[group]])){
    # i <- unique(group_df[[group]])[[1]]
    print(paste("Group:",i))
    sel_i <- group_df$sample[group_df[[group]] == i]
    sel_i <- sel_i[sel_i %in% colnames(sel_otu_table)]
    sel_otu_table_lg <- sel_otu_table[, sel_i, drop = FALSE]


    cat("Base correlation is NULL, computing base correlation...\n")
    base_correlation <- corMicro(sel_otu_table_lg, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = "ssn", vose = FALSE)[[1]]
    base_correlation[is.na(base_correlation)] <- 0
    sel_samples <- colnames(sel_otu_table_lg)
    total_num <- length(sel_samples)
    # Construct single-sample PCC for each sample
    for (num in seq_len(ncol(sel_otu_table_lg))) {
      # num <- 1
      sampleID <- colnames(sel_otu_table_lg)[num]
      cat("Processing Sample:", sampleID, "\n")

      sample_table <- sel_otu_table_lg[, -num, drop = FALSE]
      sample_correlation <- corMicro(sample_table, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = sampleID, vose = FALSE)[[1]]
      sample_correlation[is.na(sample_correlation)] <- 0

      # 限制 PCC 值范围
      base_correlation <- pmin(pmax(base_correlation, -(1-1e-10)), (1-1e-10))
      sample_correlation <- pmin(pmax(sample_correlation, -(1-1e-10)), (1-1e-10))
      # Compute differential PCC and significance

      sample_correlation <- sample_correlation[rownames(base_correlation),colnames(base_correlation)]
      delta_PCC <- base_correlation - sample_correlation
      Zvals <- delta_PCC / (1 - (sample_correlation^2)) * (total_num - 2)
      Zvals_vect <- c(Zvals)
      p_value <- 2 * stats::pnorm(-abs(Zvals_vect))
      adjusted_p_value <- p.adjust(p_value, method = "BH")

      # Extract upper triangular matrix
      upper_tri <- upper.tri(delta_PCC, diag = FALSE)
      pstat_s <- data.frame(
        node_1 = rownames(delta_PCC)[row(delta_PCC)[upper_tri]],
        node_2 = colnames(delta_PCC)[col(delta_PCC)[upper_tri]],
        p_value = adjusted_p_value[upper_tri],
        PCC_value = delta_PCC[upper_tri]
      )
      colnames(pstat_s)[4] <- paste("PCC_value", sampleID, sep = "_")

      # Save results
      file_name <- file.path(output_dir, paste0("ssPCC_", sampleID, ".tsv"))
      data.table::fwrite(pstat_s, file = file_name, sep = "\t")
    }
    base_correlation <- NULL

  }

}
