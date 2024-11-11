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
sspcc_cal <- function(sel_otu_table) {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Create output directory
  output_dir <- file.path("ssn", "sspcc")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Process input data
  sel_samples <- colnames(sel_otu_table)
  sel_otu_table_lg <- log1p(sel_otu_table)
  base_correlation <- corMicro(sel_otu_table_lg, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = "ssn", vose = FALSE)[[1]]

  # Filter base correlation matrix
  base_correlation[is.na(base_correlation)] <- 0
  total_num <- length(sel_samples)

  # Construct single-sample PCC for each sample
  for (num in seq_len(ncol(sel_otu_table))) {
    sampleID <- colnames(sel_otu_table)[num]
    cat("Processing Sample:", sampleID, "\n")

    sample_table <- sel_otu_table_lg[, -num, drop = FALSE]
    sample_correlation <- corMicro(sample_table, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = sampleID, vose = FALSE)[[1]]

    # Compute differential PCC and significance
    sample_correlation[is.na(sample_correlation)] <- 0
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
sspcc_cal2 <- function(sel_otu_table, group_df, ck = "CK", group = "group") {
  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required but not installed.")
  }

  # Create output directory
  output_dir <- file.path("ssn", "sspcc")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Check and process group information
  group_df[[group]] <- as.factor(group_df[[group]])
  sel_ck <- group_df$sample[group_df[[group]] == ck]
  sel_ck <- sel_ck[sel_ck %in% colnames(sel_otu_table)]
  ck_otu_table <- sel_otu_table[, sel_ck, drop = FALSE]
  ck_otu_table_lg <- log1p(ck_otu_table)

  # Compute base PCC matrix for control group
  ck_matrix <- corMicro(ck_otu_table_lg, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = ck, vose = FALSE)[[1]]
  sel_otu_table <- sel_otu_table[, !colnames(sel_otu_table) %in% sel_ck, drop = FALSE]
  sel_otu_table_lg <- log1p(sel_otu_table)

  total_num <- ncol(ck_otu_table)
  sel_samples <- colnames(sel_otu_table)

  for (num in seq_len(ncol(sel_otu_table))) {
    sampleID <- colnames(sel_otu_table)[num]
    cat("Processing Sample:", sampleID, "\n")

    # Construct sample matrix including control group
    sample_table <- cbind(ck_otu_table, sel_otu_table_lg[, num, drop = FALSE])
    sample_correlation <- corMicro(sample_table, method = "pearson", r.threshold = 0, p.threshold = 1, sel_group = sampleID, vose = FALSE)[[1]]

    # Compute differential PCC and significance
    sample_correlation <- sample_correlation[rownames(ck_matrix), rownames(ck_matrix)]
    delta_PCC <- sample_correlation - ck_matrix
    Zvals <- delta_PCC / (1 - (ck_matrix^2)) * (total_num - 1)
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
}
