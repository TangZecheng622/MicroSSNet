#'  与随机生成的网络进行全局属性比较
#'
#'  @param graph 输入需要比较的网络
#'  @param type 随机网络的类型
#'  @param step 迭代次数
#'  @param netName 随机网络的名字
#'  @param ncpus 运行设置的线程
#'  @export

grobal_pro_compare = function(
    graph = graph,
    type = "gnm",
    step = 100,
    netName = "CK",
    ncpus = ncpus
){


  graph = delete_vertices(graph, which(degree(graph)==0) )
  ###boot###
  num = c(1:step)


  boot <- function(i, graph,type){

    # library(igraph)
    # source('C:/Users/86135/Desktop/script/R/network/4properties.R')
    # source('C:/Users/86135/Desktop/script/R/network/modularity.R')
    rand_g <- erdos.renyi.game(length(V(graph)), length(E(graph)), type = type)
    tem_netpro_result <- t(node_properties(rand_g)[[2]])
    tem_netpro_result <- as.data.frame(tem_netpro_result)
    colnames(tem_netpro_result) <- "value"
    tem_netpro_result$value <- as.numeric(tem_netpro_result$value)
    tem_netpro_result <- as.matrix(tem_netpro_result)

  }

  cl <- makeCluster(ncpus)
  clusterExport(cl, c("graph","type"),envir = environment())  # 将对象'g'导出到每个节点
  dfp <- parLapply(cl, num, boot, graph = graph , type = type)
  stopCluster(cl)

  rand_g_netpro_result <- do.call(cbind, dfp)




  # for (i in 1:step){
  #   # Random null model
  #   rand_g <- igraph::erdos.renyi.game(length(igraph::V(graph)), length(igraph::E(graph)), type = type)
  #   tem_netpro_result <- t(node_properties(rand_g)[[2]])
  #   tem_netpro_result <- as.data.frame(tem_netpro_result)
  #   colnames(tem_netpro_result) <- "value"
  #   tem_netpro_result$value <- as.numeric(tem_netpro_result$value)
  #   tem_netpro_result <- as.matrix(tem_netpro_result)
  #   rand_g_netpro_result<-cbind(rand_g_netpro_result,tem_netpro_result)
  # }

  result_summary <- cbind(rowMeans(rand_g_netpro_result,na.rm = TRUE), apply(rand_g_netpro_result, 1, function(x) sd(x, na.rm = TRUE)))
  colnames(result_summary) <- c("Means", "SD")
  head(result_summary)

  netpro_result <- t(node_properties(graph)[[2]])
  colnames(netpro_result) <- netName
  sum_net <- cbind(netpro_result, result_summary)

  fitinf <- fit_poweRlaw(graph)
  if (fitinf[[1]] <= 0.05) {
    print("This network may be a Scale-free network!")
  } else if (sum_net[3, netName] > sum_net[3, "Means"]) {
    print("This network may be a Small-World network!")
  } else {
    print("This network may be a random network or some other type of network!")
  }

  return(list(sum_net,fitinf[[2]]))
}


fit_poweRlaw = function(graph){
  degree_dist <- table(igraph::degree(graph))
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)

  dat <- data.frame(degree = degree_num, count = degree_count)

  par(mfrow = c(1, 3))
  hist(degree_dist, ylab = 'Number of degrees', xlab = 'Frequency',
       main = 'Degree distribution')
  plot(degree_num,degree_count,  xlab = 'Degree', ylab = 'Count',
       main = 'Degree distribution')
  plot(degree_num,degree_count,  log = 'xy', xlab = 'Log-degree',
       ylab = 'Log-count', main = 'Log-log degree distribution')

  log_degree <- log(dat$degree)
  log_count <- log(dat$count)
  linear_model <- lm(log_count ~ log_degree)
  summary(linear_model)

  coefficients <- coef(linear_model)
  initial_a <- exp(coefficients[1])
  initial_b <- coefficients[2]
  mod <- nls(count ~ a*degree^b, data = dat, start = list(a = initial_a, b = initial_b),
             control = list(maxiter = 1000))
  summary(mod)
  a <- round(coef(mod)[1], 3)
  b <- round(coef(mod)[2], 3)

  fit <- fitted(mod)
  SSre <- sum((dat$count-fit)^2)
  SStot <- sum((dat$count-mean(dat$count))^2)
  R2 <- round(1 - SSre/SStot, 3)
  print(R2)

  p_num <- 1
  dat_rand <- dat
  for (i in 1:999) {
    dat_rand$count <- sample(dat_rand$count)
    SSre_rand <- sum((dat_rand$count-fit)^2)
    SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
    R2_rand <- 1 - SSre_rand/SStot_rand
    if (R2_rand > R2) p_num <- p_num + 1
  }
  p_value <- p_num / (999+1)
  print(p_value)


  p <- ggplot(dat, aes(x = degree, y = count)) +
    geom_point(color = 'blue') +
    theme(panel.grid.major = element_line(color = 'gray'), panel.background = element_rect(color = 'black', fill = 'transparent')) +
    stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = a, b = b),control = list(maxiter = 1000)), se = FALSE) +
    labs(x = 'Degree', y = 'Count')

  p <- p + xlim(min(dat$degree, 0)*1.2, max(dat$degree)*1.2)


  #添加公式拟合的注释
  label <- data.frame(
    formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
    R2 = sprintf('italic(R^2) == %.3f', R2),
    p_value = sprintf('italic(P) < %.3f', p_value)
  )

  x_max <- round(8*max(dat$degree)/9)
  y_max <- round(9*max(dat$count)/10)

  p <- p + geom_text(x = x_max, y = y_max, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
    geom_text(x = x_max, y = (10*y_max/11), aes(label = R2), data = label, parse = TRUE, hjust = 0) +
    geom_text(x = x_max, y = (9*y_max/11), aes(label = p_value), data = label, parse = TRUE, hjust = 0)

  return(list(p_value,p))
}


