#' Calculate Positive Cohesion
#'
#' Implements the positive cohesion metric defined by Herren & McMahon (2017)
#' by computing corrected correlations via null model permutations.
#'
#' @param abund_table OTU abundance table (rows = taxa, columns = samples)
#' @param n_perm Number of permutations for null model
#' @param seed Random seed for reproducibility
#' @return A named numeric vector of sample-level positive cohesion values
#' @export

# 主函数：计算 positive cohesion（使用 column shuffle null model）
calculate_positive_cohesion <- function(abund_table,ncpus=1,
                                        n_perm = 200,seed = 123) {
  # abund_table: 行为物种，列为样本，相对丰度矩阵
  if (!is.data.frame(abund_table) && !is.matrix(abund_table)) stop("'abund_table' must be a matrix or data frame.")
  if (any(is.na(abund_table))) warning("Missing values found in abundance data.")
  if (any(abund_table < 0)) warning("Negative values found; make sure abundances are non-negative.")
  
  set.seed(seed)
  
  # 转置矩阵为数据框
  abund <- t(abund_table)                    # 行=样本,列=OTU
  n_taxa <- ncol(abund)
  
  if (ncol(abund) < 2 || nrow(abund) < 2) {
    stop("Not enough taxa or samples to compute correlation.")
  }
  
  # 计算观察相关矩阵（Pearson）
  real_corr <- cor(abund)

  # 初始化期望相关矩阵
  expected_corr <- matrix(0, n_taxa, n_taxa,
                          dimnames = list(colnames(abund), colnames(abund)))
  
  for (i in seq_len(n_taxa)) {
    # i =2
    null_corrs <- matrix(NA, nrow = n_perm, ncol = n_taxa-1)
    for (p in seq_len(n_perm)) {
      # p = 1
      perm <- abund
      perm[ , -i] <- apply(perm[ , -i], 2, sample)
      null_corrs[p, ] <- cor(perm)[i, -i]
    }
    
    exp_vec <- apply(null_corrs, 2, median, na.rm = TRUE)
    expected_corr[i, -i]  <- exp_vec
    expected_corr[-i,  i] <- exp_vec
  }
  
  # 获取 corrected correlation = observed - expected
  corrected_corr <- real_corr - expected_corr
  diag(corrected_corr) <- 0
  
  # 计算每个物种的 positive connectedness
  pos_connectedness <- apply(corrected_corr, 1, function(x) {
    mean(x[x > 0], na.rm = TRUE)
  })
  
  # 计算每个样本的 positive cohesion
  cohesion <- as.matrix(abund) %*% as.numeric(pos_connectedness)
  names(cohesion) <- rownames(abund)
  
  return(drop(cohesion))
}

# 使用示例：
# otu_table <- read.csv("otu_abundance.csv", row.names = 1)
# pos_cohesion <- calculate_positive_cohesion(otu_table)
# print(pos_cohesion)
