#' 基于微生物相关性构建微生物网络以及后续相关分析
#'
#'  @param table1 研究对象的丰度矩阵（标准化后）
#'  @param tax 分类注释信息
#'  @param cor_table_list 待分析的相关性矩阵列表
#'  @param group_df 研究对象分组信息表，NULL自动按照样品名生成，FALSE默认为只有一组样品
#'  @param pvl.threshold 物种流行度阈值
#'  @param r.threshold 相关性R值阈值
#'  @param p.threshold 相关性P值阈值
#'  @param method 相关性网络构建方法，spearman、pearson、sparcc、spieceasi四种方法
#'  @param R SSPCC构建网络迭代参数
#'  @param ncpus 运行设置的线程
#'  @param clu_method 网络聚类方法
#'  @param p.adj 网络相关性的多重矫正方法
#'  @param step 迭代次数
#'  @param top 过滤原始表格，top>=1时按照top值选取平均丰度前top，1>top>0时，选取平均丰度高于top值
#'  @param rm.p.list 模拟随机毁灭的比例
#'  @param zipi 对节点进行zipi值计算
#'  @param node.cluster 每个模块点数量阈值，大于这个阈值才上色
#'  @param conpare_net 与随机空白网络进行比较
#'  @param features 是否进行鲁棒性、脆弱性模拟
#'  @examples
#'  network_pipeline(otu,group_df = FALSE,zipi = TRUE,compare_net = TRUE,features = TRUE)
#'  @return 输出网络矩阵以及一系列分析结果和可视化
#'  @export



###输入行为物种，列为样品###
network_pipeline <- function(
    table1 ,
    tax = NULL,
    cor_table_list = NULL,
    group_df = NULL,
    pvl.threshold = 0,
    r.threshold=0.6,
    p.threshold=0.05,
    method = "spearman",
    R = 10,
    ncpus = 10,
    clu_method="cluster_fast_greedy",
    p.adj = "fdr",
    step = 100,
    top = NULL,
    rm.p.list = seq(0.05, 0.95, by = 0.05),
    zipi = FALSE,
    compare_net = FALSE,
    features = FALSE
){

  rrob_sum_list <- list()
  rrob_detail_list <- list()


  ###检测是否包含group文件####group要求包含group与sample
  if(!is.null(group_df)){
    strings <- colnames(group_df)
    sample1 <- grepl("sample", strings, ignore.case = TRUE)
    group1 <- grepl("group", strings, ignore.case = TRUE)
    if (any(sample1)) {
      sample2 <- strings[sample1][[1]]
    } else {
      stop("Error: 'sample' column not found in group_df")
    }

    if (any(group1)) {
      group2 <- strings[group1]
    } else {
      stop(paste("Error:", vscol1, "column not found in group_df"))
    }
    ###获取group###
    group_list <- unique(group_df$group)
  }else{
    samples <- colnames(table1)
    # 创建数据框并生成分组名
    group_df <- data.frame(
      sample = samples,
      group = as.factor(gsub("[[:digit:]]|[[:punct:]]", "", samples))
    )
    group_list <- unique(group_df$group)
  }


  for (sel_group in group_list){
    ###提示迭代group###
    cat("Processing group:",sel_group,"\n")
    sel_species_df <- table1[,group_df$sample[group_df$group==sel_group]] #提取指定group丰度表
    sel_species_df1 <- sel_species_df[rowSums(sel_species_df>0)>0,] #提取总行不为0的丰度表
    sel_species_df1 <- filter_OTU2(sel_species_df1,Pre = pvl.threshold) #按流行率过滤
    sel_species_df1 <- filter_OTU(sel_species_df1,Top = top) #按照平均丰度阈值或者前多少个平均丰度TOP过滤

    if(!is.null(cor_table_list)){
      cor <- cor_table_list[[sel_group]]
    }else{
      cor_pval <- corMicro(otu = sel_species_df1,p.threshold = p.threshold,method = method,R = R,ncpus = ncpus)[[1]]
    }


    ###构建网络矩阵
    g <-  graph_from_adjacency_matrix(as.matrix(cor), weighted = TRUE, mode = 'undirected')#构建邻接矩阵
    E(g)$correlation <- E(g)$weight #赋予正负相关性
    E(g)$weight <- abs(E(g)$weight) #不考虑正负的相关性
    write.graph(g, paste0(sel_group,'.net.graphml'), format = 'graphml')
    write.graph(g, paste0(sel_group,'.net.gml'), format = 'gml')
    ###计算属性
    pro <- node_properties(g,clu_method=clu_method,zipi=zipi,tag = sel_group) #计算网络属性
    vul <- Vulnerability(cor)
    ###可添加一个度丰度的图

    local_pro <- pro[[1]] #局部属性
    global_pro <- pro[[2]] #全局属性
    global_pro$vul <- vul
    write.csv(local_pro, paste(sel_group,"_local_properties.csv"), row.names = FALSE)
    write.csv(global_pro, paste(sel_group,"_global_properties.csv"), row.names = FALSE)

    ###可视化
    if(!is.null(tax)){
      plotnetwork_tax(g = g,tag = sel_group,tax = tax,df = sel_species_df1)
    }else{
      plotnetwork(g = g,clu_method = clu_method,tag = sel_group,df = sel_species_df1,node.cluster = node.cluster)
    }

    ###空模型比对
    if(compare_net){
      compare_rmc <- grobal_pro_compare(graph = g,step = 100,netName = "Control",ncpus = ncpus)
      compare_table <- compare_rmc[[1]]
      compare_p <- compare_rmc[[2]]
      write.csv(compare_table, paste(sel_group,"_compare_rmc.csv"), row.names = TRUE)
      ggsave(filename= paste(sel_group,"_compare_rmc_Degree_Distribution.pdf"),plot=compare_p,width = 10,height = 6,dpi = 1000)
    }

    ###鲁棒性

    if(features){
      random_rob <- Robustness.Random.removal(cor = cor,otu = sel_species_df1,tag = sel_group,rm.p.list = rm.p.list,nperm = step,ncpus = ncpus)

      rrob_sum_list[[paste0(sel_group,"_rrob_sum")]] <- random_rob[1]
      rrob_detail_list[[paste0(sel_group,"_rrob_detail")]] <- random_rob[2]

      target_rob <- Robustness.Targeted.removal(cor = cor,otu = sel_species_df1,tag = sel_group,loc_pro  = local_pro,bet_centrality = FALSE,clo_centrality = FALSE,
                                                pagerank = FALSE,eig_centrality = FALSE,kcore = FALSE)
    }
  }
  ###组组之间比较
  compare_rob_vul(rrob_sum_list,rrob_detail_list)
}


