#' 基于微生物相关性构建微生物网络
#'
#'  @param table1 过滤后的otu表
#'  @param method 相关性网络构建方法，spearman、pearson、sparcc、spieceasi四种方法
#'  @param p.adj 网络相关性的多重矫正方法
#'  @param sel_group 所需要分析的组
#'  @param R SSPCC构建网络迭代参数
#'  @param ncpus 运行设置的线程
#'  @param r.threshold 相关性R值阈值
#'  @param p.threshold 相关性P值阈值
#'  @param vose 是否输出相关文件
#'  @export

corMicro = function(table1 ,
                    method = "spearman",
                    p.adj = "BH",
                    sel_group = "total",
                    R = 10,
                    ncpus = 10,
                    r.threshold = 0.3,
                    p.threshold = 0.05,
                    vose = TRUE
                    ){


  if (method %in% c("pearson","spearman")) {

    occor <- rcorr(t(table1),type = method)
    occor.r <- occor$r
    occor.P <- occor$P
    diag(occor.r) = 0
    diag(occor.P) = 1

    p_values_vector <- occor.P[upper.tri(occor.P)]
    p_adj_vector <- p.adjust(p_values_vector, method = p.adj)
    occor.p <- occor.P  # 复制原矩阵结构
    occor.p[upper.tri(occor.p)] <- p_adj_vector  # 填回上三角的校正结果
    occor.p[lower.tri(occor.p)] <- t(occor.p)[lower.tri(occor.p)]

    ###输出矩阵
    occor.r[abs(occor.r)<r.threshold|occor.p>p.threshold] = 0
    diag(occor.r) = 0
    occor.r[occor.r>1] <- 1
    if (vose == TRUE){
      write.csv(occor.r,file = paste0(sel_group,"_ori_cor_martix.csv")) ###输出原始相关性矩阵文件
      write.csv(occor.p,file = paste0(sel_group,"_ori_pval_martix.csv")) ###输出原始相关性矩阵文件
      write.csv(occor.r,file = paste0(sel_group,"cor_martix.csv")) ###输出过滤后相关性矩阵文件
    }
    cor = occor.r
  }

  if (method %in% c("sparcc")) {

    result <- sparcc.micro(data = t(table1),R = R,ncpus = ncpus)
    occor.r = result[[1]]
    occor.p = result[[2]]

    occor.r[abs(occor.r)<r.threshold|occor.p>p.threshold] = 0
    diag(occor.r) = 0
    cor = occor.r
    if (vose == TRUE){
      write.csv(occor.r,file = paste0(sel_group,"_sparcc_ori_cor_martix.csv")) ###输出原始相关性矩阵文件
      write.csv(occor.p,file = paste0(sel_group,"_sparcc_ori_pval_martix.csv")) ###输出原始相关性矩阵文件
      write.csv(occor.r,file = paste0(sel_group,"_sparcc_cor_martix.csv")) ###输出过滤后相关性矩阵文件
    }

  }

  if(method %in% c("SpiecEasi")){
    occor <- spiec.easi(data=as.matrix(t(table1)),method='glasso',lambda.min.ratio=0.01,
                              nlambda=20,pulsar.params=list(rep.num=50,thresh = p.threshold))
    adjacency_unweight<-data.frame(as.matrix(occor$refit$stars))
    rownames(adjacency_unweight)<-colnames(t(table1))
    colnames(adjacency_unweight)<-colnames(t(table1))

    cor <- adjacency_unweight
    if (vose == TRUE){
      write.table(adjacency_unweight,paste0(sel_group,'_SpiecEasi_unweight.glasso.tsv'),col.names=NA,sep='\t',quote=FALSE)
    }

  }

  cor_list <- get_full_edge_list(cor)
  return(list(cor,cor_list))

}


sparcc.micro <- function(
    data = data,
    R = 10,
    ncpus = 1

){
  spmatrix <- SpiecEasi::sparcc(data,iter = 20,inner_iter = 10,th = 0.1)

  sp.boot <- SpiecEasi::sparccboot(
    data,
    R = R,
    ncpus = ncpus
  )

  sp.p <- SpiecEasi::pval.sparccboot(sp.boot, sided = "both")
  cors <- sp.p$cors
  sp.p$pvals[is.na(sp.p$pvals)] = 1
  pvals <- sp.p$pvals

  sparCCpcors <- diag(0.5, nrow = dim(spmatrix$Cor)[1], ncol = dim(spmatrix$Cor)[1])
  sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
  sparCCpcors <- sparCCpcors + t(sparCCpcors)

  sparCCpval <- diag(0.5, nrow = dim(spmatrix$Cor)[1], ncol = dim(spmatrix$Cor)[1])
  sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
  sparCCpval <- sparCCpval + t(sparCCpval)
  dim(sparCCpval)

  rownames(sparCCpcors) <- colnames(data)
  colnames(sparCCpcors) <- colnames(data)
  rownames(sparCCpval) <- colnames(data)
  colnames(sparCCpval) <- colnames(data)

  reordered_all_sparcc <- reorder_cor_and_p(sparCCpcors, sparCCpval)
  occor.r <- reordered_all_sparcc$r
  occor.p <- reordered_all_sparcc$p


  return(list(occor.r,occor.p))

}
reorder_cor_and_p <- function(cormat, pmat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  pmat <- pmat[hc$order, hc$order]
  list(r = cormat, p = pmat)
}
get_full_edge_list <- function(mat) {
  # 提取矩阵的所有元素（包括对角线、上三角和下三角）的行和列索引
  all_pairs <- which(!is.na(mat), arr.ind = TRUE)

  # 获取节点名称
  node1 <- rownames(mat)[all_pairs[,1]]
  node2 <- colnames(mat)[all_pairs[,2]]

  # 获取对应的权重
  weights <- mat[all_pairs]

  # 组合成一个数据框
  edge_list <- data.frame(Node1 = node1, Node2 = node2, Weight = weights)

  return(edge_list)
}


