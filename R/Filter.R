#'  过滤、筛选OTU
#'  @param otu_table 过滤后的otu表
#'  @param Top 过滤原始表格，top>=1时按照top值选取平均丰度前top，1>top>0时，选取平均丰度高于top值
#'  @param Pre 物种流行度阈值
#'  @export
#行为otu,列为样品

filter_OTU <- function(otu_table, Top = NULL) {
  if (!is.null(Top) && Top != 0) {

    # 转为相对丰度
    relative_abundance <- sweep(otu_table, 2, colSums(otu_table), FUN = "/")
    relative_abundance$mean <- rowMeans(relative_abundance)
    relative_abundance$ID <- rownames(relative_abundance)

    # 根据相对丰度均值降序排列
    relative_abundance <- dplyr::arrange(relative_abundance, desc(mean))


    if (Top > 0 && Top < 1) {
      # 当 Top 为小数时，提取 mean 大于这个小数的 OTU
      subtab <- relative_abundance[relative_abundance$mean > Top, ]
    } else {
      # 当 Top 为整数时，提取前 Top 个 OTU
      subtab <- head(relative_abundance, Top)
    }

    # 返回筛选后的 OTU 表
    otu_table_filtered <- otu_table[subtab$ID, , drop = FALSE]
    return(otu_table_filtered)
  } else {
    # 如果 Top 为 0 或 NULL，返回原始 OTU 表
    return(otu_table)
  }
}

###行为otu，列为样品
###根据流行率筛选

filter_OTU2 <- function(otu_table, Pre = NULL) {

  if (Pre<0 || Pre>1){
    stop("Error: x must be between 0 and 1.")
  }
  if (!is.null(Pre) && Pre != 0) {
    prevalence <- rowSums(otu_table>0)
    otu_table_filtered <- otu_table[prevalence>= round(ncol(otu_table)*Pre),]
    return(otu_table_filtered)
  } else {
    # 如果 Top 为 0 或 NULL，返回原始 OTU 表
    return(otu_table)
  }
}
