#' 与随机生成的网络进行全局属性比较
#'
#' @param graph `igraph`对象，待比较的网络
#' @param type 字符串，随机网络的类型，默认为 "gnm"
#' @param step 整数，随机网络的迭代次数，默认为 100
#' @param netName 字符串，网络名称，默认为 "CK"
#' @param ncpus 整数，运行并行计算时的 CPU 数量
#' @return 列表，包含实际网络与随机网络的全局属性对比结果和拟合模型的 p 值
#' @importFrom igraph delete_vertices degree V E
#' @importFrom parallel makeCluster stopCluster clusterExport clusterEvalQ parLapply
#' @export

grobal_pro_compare <- function(graph, type = "gnm", step = 100, netName = sel_group, ncpus = 1) {
  # 参数检查
  if (!inherits(graph, "igraph")) stop("Error: 'graph' must be an 'igraph' object.")
  if (!is.character(type) || !(type %in% c("gnm", "other_type"))) stop("Error: 'type' must be 'gnm' or another valid type.")
  if (!is.numeric(step) || step <= 0) stop("Error: 'step' must be a positive integer.")
  if (!is.character(netName)) stop("Error: 'netName' must be a character string.")
  if (!is.numeric(ncpus) || ncpus <= 0) stop("Error: 'ncpus' must be a positive integer.")

  # graph <- g
  # 删除度为0的节点
  graph <- delete_vertices(graph, which(degree(graph) == 0))

  # 并行计算随机网络属性
  cluster <- parallel::makeCluster(ncpus)
  parallel::clusterExport(cluster, c("graph", "netName","type", "node_properties", "generate_random_network", "modularity_igraph"), envir = environment())
  parallel::clusterEvalQ(cluster, library(igraph))  # 在每个子进程中加载 igraph 包

  boot <- function(i) {

    rand_g <- generate_random_network(graph, type)
    properties <- t(node_properties(rand_g,tag = netName)[[2]])
    as.data.frame(properties)
  }

  rand_g_netpro_result <-  parallel::parLapply(cluster, 1:step, boot)
  parallel::stopCluster(cluster)

  # 转换结果为矩阵并计算均值和标准差
  rand_g_netpro_result <- do.call(cbind, rand_g_netpro_result)
  # colnames(rand_g_netpro_result) <- rand_g_netpro_result[1,]
  rand_g_netpro_result <- as.data.frame(lapply(rand_g_netpro_result[-1,], as.numeric))
  result_summary <- data.frame(
    Means = rowMeans(rand_g_netpro_result, na.rm = TRUE),
    SD = apply(rand_g_netpro_result, 1, sd, na.rm = TRUE)
  )

  # 计算实际网络的属性
  netpro_result <- t(node_properties(graph)[[2]])
  colnames(netpro_result) <- netName
  netpro_result <- as.data.frame(netpro_result[-1, , drop = FALSE])

  # 合并实际网络与随机网络的结果
  sum_net <- cbind( result_summary, netpro_result)

  # 网络类型推断
  fitinf <- fit_poweRlaw(graph)
  if(is.na(fitinf[[1]])) {

    message("Fitting failed: Unable to determine network type due to fitting issues.")
  } else if (fitinf[[1]] <= 0.05) {
    message("This network may be a Scale-free network!")
  } else if (sum_net[3, netName] > sum_net[3, "Means"]) {
    message("This network may be a Small-World network!")
  } else {
    message("This network may be a random network or some other type of network!")
  }

  return(list(sum_net, fitinf[[2]]))
}

#' 生成随机网络
#'
#' @param graph `igraph`对象，输入的网络
#' @param type 字符串，随机网络的类型
#' @return `igraph`对象，生成的随机网络
#' @importFrom igraph erdos.renyi.game V E
generate_random_network <- function(graph, type = "gnm") {
  if (type == "gnm") {
    rand_g <- igraph::erdos.renyi.game(length(igraph::V(graph)), length(igraph::E(graph)), type = "gnm")
  } else {
    # 其他类型的随机网络生成逻辑
    rand_g <- igraph::erdos.renyi.game(length(igraph::V(graph)), length(igraph::E(graph)), type = "gnm") # 默认 gnm
  }
  return(rand_g)
}

#' 计算网络的幂律拟合信息
#'
#' @param graph `igraph`对象，输入的网络
#' @return 列表，包含 p 值和拟合模型的图形
#' @importFrom igraph degree
#' @importFrom stats lm coef nls fitted
#' @importFrom ggplot2 ggplot aes geom_point stat_smooth labs theme_minimal annotate
fit_poweRlaw <- function(graph) {
  degree_dist <- table(igraph::degree(graph))
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)
  dat <- data.frame(degree = degree_num, count = degree_count)
  dat <- dat[dat$degree > 0 & dat$count > 0, ]

  # 拟合幂律分布
  log_degree <- log(dat$degree)
  log_count <- log(dat$count)
  linear_model <- lm(log_count ~ log_degree)
  coefficients <- coef(linear_model)
  initial_a <- exp(coefficients[1])
  initial_b <- coefficients[2]

  # 非线性拟合
  mod <- tryCatch({
    nls(count ~ a * degree^b, data = dat, start = list(a = initial_a, b = initial_b), control = list(maxiter = 1000))
  }, error = function(e) {
    message("Nonlinear fitting failed: ", e$message)
    NULL  # 返回 NULL 以便在上层逻辑中处理拟合失败的情况
  })

  if (is.null(mod)) return(list(p_value = NA, plot = NULL))

  a <- coef(mod)[1]
  b <- coef(mod)[2]

  # 计算 R^2 和 p 值
  fit <- fitted(mod)
  SSre <- sum((dat$count - fit)^2)
  SStot <- sum((dat$count - mean(dat$count))^2)
  R2 <- 1 - SSre / SStot

  # 计算 p 值
  p_num <- 1
  for (i in 1:999) {
    dat_rand <- dat
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count - fit)^2)
    SStot_rand <- sum((dat_rand$count - mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand / SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / 1000

  # 绘制拟合曲线
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = degree, y = count)) +
    ggplot2::geom_point(color = 'blue') +
    ggplot2::stat_smooth(method = 'nls', formula = y ~ a * x^b, method.args = list(start = list(a = a, b = b)), se = FALSE) +
    ggplot2::labs(x = 'Degree', y = 'Count') +
    ggplot2::theme_minimal() +
    ggplot2::annotate("text", x = max(dat$degree) * 0.8, y = max(dat$count) * 0.9,
             label = paste("R^2 =", round(R2, 3), "\nP =", round(p_value, 3)))

  return(list(p_value, p))
}
