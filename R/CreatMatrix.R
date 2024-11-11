#' Construct Microbial Network Based on Correlations
#'
#' This function constructs a microbial network based on different correlation methods.
#'
#' @param table1 Filtered OTU table (rows as OTUs, columns as samples).
#' @param method Correlation method: "spearman", "pearson", "sparcc", or "SpiecEasi".
#' @param p.adj Multiple testing correction method for p-values, default is "BH" (Benjamini-Hochberg).
#' @param sel_group The group to analyze.
#' @param R Iteration parameter for SparCC network construction.
#' @param ncpus Number of CPU cores to use.
#' @param r.threshold Threshold for correlation coefficient.
#' @param p.threshold Threshold for p-value.
#' @param vose Logical, whether to output the calculated correlation and p-value matrices. Default is TRUE.
#' @param output_dir Character string specifying the output directory.
#' @return A list containing the correlation matrix and the edge list.
#' @importFrom Hmisc rcorr
#' @importFrom stats p.adjust
#' @importFrom utils write.csv write.table
#' @importFrom SpiecEasi spiec.easi
#' @export
corMicro <- function(table1,
                     method = "spearman",
                     p.adj = "BH",
                     sel_group = "total",
                     R = 10,
                     ncpus = 10,
                     r.threshold = 0.3,
                     p.threshold = 0.05,
                     vose = TRUE,
                     output_dir = "./"  # 添加 output_dir 参数，默认值为当前目录
) {
  if (!is.data.frame(table1)) stop("Error: 'table1' must be a data frame.")
  if (!method %in% c("spearman", "pearson", "sparcc", "SpiecEasi")) stop("Error: Invalid method selected.")
  if (!is.numeric(r.threshold) || !is.numeric(p.threshold)) stop("Error: 'r.threshold' and 'p.threshold' must be numeric.")
  if (!is.logical(vose)) stop("Error: 'vose' must be a logical value (TRUE or FALSE).")
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cor <- NULL

  if (method %in% c("pearson", "spearman")) {
    occor <- Hmisc::rcorr(log1p(t(table1)), type = method)
    occor.r <- occor$r
    occor.P <- occor$P
    diag(occor.r) <- 0
    diag(occor.P) <- 1

    p_values_vector <- occor.P[upper.tri(occor.P)]
    p_adj_vector <- p.adjust(p_values_vector, method = p.adj)
    occor.p <- occor.P  # 复制原始矩阵结构
    occor.p[upper.tri(occor.p)] <- p_adj_vector  # 填充调整后的 p 值
    occor.p[lower.tri(occor.p)] <- t(occor.p)[lower.tri(occor.p)]

    # 过滤相关性
    occor.r[abs(occor.r) < r.threshold | occor.p > p.threshold] <- 0
    diag(occor.r) <- 0
    occor.r[occor.r > 1] <- 1

    cor <- occor.r

    if (vose == TRUE) {
      write.csv(occor.r, file = file.path(output_dir, paste0(sel_group, "_ori_cor_matrix.csv")))
      write.csv(occor.p, file = file.path(output_dir, paste0(sel_group, "_ori_pval_matrix.csv")))
      write.csv(cor, file = file.path(output_dir, paste0(sel_group, "_cor_matrix.csv")))
    }
  } else if (method == "sparcc") {
    result <- sparcc.micro(data = t(table1), R = R, ncpus = ncpus)
    occor.r <- result[[1]]
    occor.p <- result[[2]]

    occor.r[abs(occor.r) < r.threshold | occor.p > p.threshold] <- 0
    diag(occor.r) <- 0
    cor <- occor.r

    if (vose == TRUE) {
      write.csv(occor.r, file = file.path(output_dir, paste0(sel_group, "_sparcc_ori_cor_matrix.csv")))
      write.csv(occor.p, file = file.path(output_dir, paste0(sel_group, "_sparcc_ori_pval_matrix.csv")))
      write.csv(cor, file = file.path(output_dir, paste0(sel_group, "_sparcc_cor_matrix.csv")))
    }
  } else if (method == "SpiecEasi") {
    occor <- SpiecEasi::spiec.easi(data = as.matrix(t(table1)), method = 'glasso', lambda.min.ratio = 0.01,
                                   nlambda = 20, pulsar.params = list(rep.num = 50, thresh = p.threshold))
    adjacency_unweight <- as.matrix(occor$refit$stars)
    rownames(adjacency_unweight) <- colnames(t(table1))
    colnames(adjacency_unweight) <- colnames(t(table1))

    cor <- adjacency_unweight

    write.table(adjacency_unweight, file = file.path(output_dir, paste0(sel_group, '_SpiecEasi_unweight_glasso.tsv')),
                col.names = NA, sep = '\t', quote = FALSE)
  }

  cor_list <- get_full_edge_list(cor)
  return(list(cor = cor, cor_list = cor_list))
}



#' SparCC Microbial Correlation
#'
#' This function calculates the SparCC correlation matrix and p-values.
#'
#' @param data Numeric matrix, input data with features in columns and samples in rows.
#' @param R Integer, number of bootstrap replicates.
#' @param ncpus Integer, number of CPU cores to use.
#' @return A list containing the correlation matrix and the p-value matrix.
#' @importFrom SpiecEasi sparcc sparccboot pval.sparccboot
#' @export
sparcc.micro <- function(data, R = 10, ncpus = 1) {
  if (!is.matrix(data)) stop("Error: 'data' must be a matrix.")

  spmatrix <- SpiecEasi::sparcc(data, iter = 20, inner_iter = 10, th = 0.1)

  sp.boot <- SpiecEasi::sparccboot(
    data,
    R = R,
    ncpus = ncpus
  )

  sp.p <- SpiecEasi::pval.sparccboot(sp.boot, sided = "both")
  cors <- sp.p$cors
  sp.p$pvals[is.na(sp.p$pvals)] <- 1
  pvals <- sp.p$pvals

  n <- dim(spmatrix$Cor)[1]
  sparCCpcors <- matrix(0, nrow = n, ncol = n)
  sparCCpcors[upper.tri(sparCCpcors)] <- cors
  sparCCpcors <- sparCCpcors + t(sparCCpcors)

  sparCCpval <- matrix(1, nrow = n, ncol = n)
  sparCCpval[upper.tri(sparCCpval)] <- pvals
  sparCCpval <- sparCCpval + t(sparCCpval)

  rownames(sparCCpcors) <- colnames(data)
  colnames(sparCCpcors) <- colnames(data)
  rownames(sparCCpval) <- colnames(data)
  colnames(sparCCpval) <- colnames(data)

  reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
  occor.r <- reordered_all_sparcc$r
  occor.p <- reordered_all_sparcc$p

  return(list(occor.r, occor.p))
}


#' Reorder Correlation and P-value Matrices
#'
#' This function reorders the correlation and p-value matrices based on hierarchical clustering.
#'
#' @param cormat Correlation matrix.
#' @param pmat P-value matrix.
#' @return A list containing reordered correlation and p-value matrices.
#' @importFrom stats as.dist hclust
#' @export
reorder_cor_and_p <- function(cormat, pmat) {
  dd <- as.dist((1 - cormat) / 2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}
