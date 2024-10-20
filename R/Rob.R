#'  网络模块化
#'  @param cor 输入网络相关性矩阵
#'  @param otu 输入物种丰度表
#'  @param tag 网络名字
#'  @param nperm 迭代次数
#'  @param rm.p.list 模拟随机毁灭的比例
#'  @param ncpus 运行设置的线程
#'  @param loc_pro 网络属性
#'  @param zipi 是否根据zipi值作为target
#'  @param bet_centrality 是否根据bet_centralityi值作为target
#'  @param clo_centrality 是否根据clo_centrality值作为target
#'  @param eig_centrality 是否根据eig_centrality值作为target
#'  @param pagerank 是否根据pagerank值作为target
#'  @param kcore 是否根据kcore值作为target
#'  @param degree 是否根据degree值作为target
#'  @export
Robustness.Random.removal <- function(
    cor,
    otu,
    tag=tag,
    nperm = 100,
    rm.p.list = seq(0.05, 1, by = 0.05),
    ncpus = ncpus
){
  cor1 <- cor
  diag(cor1) <- 0
  comm <- t(otu)


  sp.ra <- colMeans(comm)
  cor1_raw <- cor1[colSums(abs(cor1))>0,colSums(abs(cor1))>0]
  sp.ra2 <- sp.ra[colSums(abs(cor1)) > 0]
  sum(row.names(cor1_raw)==names(sp.ra2))  #check if matched


  Vulnerability <- function(
    cor){
    cor1 <- cor
    diag(cor1)<-0
    cor1[abs(cor1)>0]<-1
    g = graph_from_adjacency_matrix(as.matrix(cor1), mode="undirected", weighted = NULL, diag = FALSE, add.colnames = NULL)

    #移除孤立节点
    iso_node_id = which(degree(g)==0)
    g2 = delete_vertices(g, iso_node_id)
    length(V(g2));length(E(g2))
    network.efficiency <- function(graph){
      if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
      dd <- 1/distances(graph)
      diag(dd) <- NA
      efficiency <- mean(dd, na.rm=T)
      #denom <- nrow(dd)*(ncol(dd)-1)
      #sum(dd, na.rm=T)/denom
      return(efficiency)
    }

    info.centrality.vertex <- function(graph, net=NULL, verbose=F){
      if(is_igraph(graph)==F) warning("Please use a valid iGraph object")
      if(is.null(net)) net <- network.efficiency(graph)
      if(is.numeric(net)==F){
        warning("Please ensure net is a scalar numeric")
        net <- network.efficiency(graph)
      }
      count <- c()
      for(i in 1:length(V(graph))){
        count <- c(count, (net-network.efficiency(delete.vertices(graph, i)))/net)
        if(verbose){
          print(paste("node",i,"current\ info\ score", count[i], collapse="\t"))
        }
      }
      return(count)
    }
    #计算每个节点的漏洞
    node.vul<-info.centrality.vertex(g2)
    #计算出该网络的最大节点易损性
    max(node.vul)

  }

  boot <- function(i, netRaw, rm.percent, sp.ra, abundance.weighted) {
    library(igraph)
    id.rm <- sample(1:nrow(netRaw), round(nrow(netRaw) * rm.percent))
    net.Raw <- netRaw  # 不要改变原始 netRaw
    net.Raw[id.rm, ] <- 0
    net.Raw[, id.rm] <- 0  # 移除这些物种的所有连结
    if (abundance.weighted) {
      net.stength <- net.Raw * sp.ra
    } else {
      net.stength <- net.Raw
    }
    sp.meanInteration <- colMeans(net.stength)
    id.rm2 <- which(sp.meanInteration <= 0)  # 移除与其他物种无正向交互或无交互的物种
    remain.percent <- (nrow(netRaw) - length(id.rm2)) / nrow(netRaw)
    vul <- Vulnerability(cor = net.stength)
    return(list(remain = remain.percent, vul = vul))

  }


  rmsimu3 <- function(netRaw, rm.p.list, sp.ra, abundance.weighted = TRUE,nperm = 100,ncpus = ncpus) {
    dat1_results <- data.frame(
      remain.mean = numeric(0),
      remain.sd = numeric(0),
      remain.se = numeric(0)
    )
    colnames(dat1_results) <- c("remain.mean", "remain.sd", "remain.se")

    dat2_results <- data.frame()
    cl <- makeCluster(ncpus)
    clusterExport(cl, c("netRaw","sp.ra","abundance.weighted","boot","Vulnerability"),envir = environment())  # 将对象'g'导出到每个节点

    for (rm.percent in rm.p.list){
      num = c(1:nperm)
      results_list <- parLapply(cl, num, function(i) boot(i, netRaw, rm.percent, sp.ra, abundance.weighted))
      Proportion.remain <- sapply(results_list, `[[`, "remain")
      Vulnerability <- sapply(results_list, `[[`, "vul")

      remain.mean <- mean(Proportion.remain)
      remain.sd <- sd(Proportion.remain)
      remain.se <- sd(Proportion.remain) / (nperm^0.5)
      # remain_result <- c(remain.mean, remain.sd, remain.se)
      # names(remain_result) <- c("remain.mean", "remain.sd", "remain.se")
      remain_result <- data.frame(
        remain.mean = remain.mean,
        remain.sd = remain.sd,
        remain.se = remain.se
      )
      dat1_results <- rbind(dat1_results,remain_result)



      dat2 <- data.frame(
        group = rep(tag,nperm),
        iteration = 1:nperm,
        Proportion.removed = rep(rm.percent, nperm),
        Proportion.remain = Proportion.remain,
        Vulnerability = Vulnerability
        )
      dat2_results <- rbind(dat2_results,dat2)

    }


      stopCluster(cl)
      return(list(dat1_results,dat2_results))

  }

  Unweighted.simu3 <- rmsimu3(netRaw = cor1_raw, rm.p.list = rm.p.list, sp.ra = sp.ra2, abundance.weighted = FALSE, nperm = nperm,ncpus = ncpus)
  Weighted.simu3 <- rmsimu3(netRaw = cor1_raw, rm.p.list = rm.p.list, sp.ra = sp.ra2, abundance.weighted = TRUE, nperm = nperm,ncpus = ncpus)

  rep_times <- length(rm.p.list)
  # 整理模拟结果
  dat1 <- data.frame(
    Proportion.removed = rep(rm.p.list, 2),
    rbind(Weighted.simu3[[1]], Unweighted.simu3[[1]]),
    weighted = rep(c("weighted", "unweighted"), each = rep_times),
    group = rep(tag,2*rep_times))

  write.csv(dat1, paste("./",tag,"_random_deletion_results_sum.csv"))


  dat2 <- data.frame(
    rbind(Weighted.simu3[[2]],Unweighted.simu3[[2]]),
    weighted = rep(c("weighted", "unweighted"), each = rep_times*nperm))
  write.csv(dat2, paste("./",tag,"_random_deletion_results_detail.csv"))

  return(list(dat1,dat2))

}

Robustness.Targeted.removal <- function(
    cor,
    otu,
    tag=tag,
    loc_pro=local_pro,
    degree = TRUE,
    zipi = TRUE,
    bet_centrality = TRUE,
    clo_centrality = TRUE,
    eig_centrality = TRUE,
    pagerank = TRUE,
    kcore = TRUE

){
  cor1 <- cor
  diag(cor1) <- 0
  comm <- t(otu)

  sp.ra <- colMeans(comm)
  cor1_raw <- cor1[colSums(abs(cor1))>0,colSums(abs(cor1))>0]
  sp.ra2 <- sp.ra[colSums(abs(cor1)) > 0]
  sum(row.names(cor1_raw)==names(sp.ra2))  #check if matched

  property <- loc_pro

  rand.remov2.once<-function(netRaw, rm.num,
                             keystonelist, sp.ra, abundance.weighted=T){
    rm.num2<-ifelse(rm.num > length(keystonelist), length(keystonelist), rm.num)
    id.rm<-sample(keystonelist, rm.num2)
    net.Raw = netRaw
    dim(net.Raw)
    net.new = net.Raw[!names(sp.ra) %in% id.rm, !names(sp.ra) %in% id.rm]   ##remove all the links to these species
    if (nrow(net.new)<2){
      0
    } else {
      sp.ra.new=sp.ra[!names(sp.ra) %in% id.rm]

      if (abundance.weighted){
        net.stength= net.new*sp.ra.new
      } else {
        net.stength= net.new
      }

      sp.meanInteration<-colMeans(net.stength)


      while ( length(sp.meanInteration)>1 & min(sp.meanInteration) <=0){
        id.remain<- which(sp.meanInteration>0)
        net.new=net.new[id.remain,id.remain]
        sp.ra.new=sp.ra.new[id.remain]

        if (abundance.weighted){
          net.stength= net.new*sp.ra.new
        } else {
          net.stength= net.new
        }

        if (length(net.stength)>1){
          sp.meanInteration<-colMeans(net.stength)
        } else{
          sp.meanInteration<-0
        }

      }

      remain.percent<-length(sp.ra.new)/length(sp.ra)

      remain.percent}
  }

  rmsimu<-function(netRaw, rm.p.list, keystonelist,sp.ra, abundance.weighted=T,nperm=100){
    t(sapply(rm.p.list,function(x){
      remains=sapply(1:nperm,function(i){
        rand.remov2.once(netRaw=netRaw, rm.num=x, keystonelist=keystonelist, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
      })
      remain.mean=mean(remains)
      remain.sd=sd(remains)
      remain.se=sd(remains)/(nperm^0.5)
      result<-c(remain.mean,remain.sd,remain.se)
      names(result)<-c("remain.mean","remain.sd","remain.se")
      result
    }))
  }



  # Weighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
  # Unweighted.simu<-rmsimu(netRaw=network.raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)
  #
  # dat1<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
  #                  weighted=rep(c("weighted","unweighted"),each=length(module.hub)))
  #
  # currentdat = dat1
  # write.csv(currentdat,"targeted_deletion_results.csv")
  #
  # return(dat1)

  if(zipi){
    model <- filter(loc_pro,type == "Module hubs")
    module.hub <- as.character(row.names(model))
    if(length(module.hub)==0){
      print("The number of module hubs provided by zipi is zero, so you might want to explore other options.")
    }else{
      Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
      Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

      dat_zipi<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                           weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

      write.csv(dat_zipi,paste("./",tag,"_targeted_deletion_results_zipi.csv"))

    }
  }

  if(degree){

    model <- property %>%
      as.data.frame() %>%
      filter(!is.na(Degree)) %>%
      arrange(desc(Degree))
    tem <- round(length(model$Degree)*0.05,0)
    module.hub = row.names(model)[1:tem]

    Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat_degree<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                           weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

    write.csv(dat_degree,paste("./",tag,"_targeted_deletion_results_degree.csv"))

    ggplot(dat_degree[dat_degree$weighted=="weighted",], aes(x=Number.hub.removed, y=remain.mean)) +
      geom_line()+
      geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
      xlab("Number of module hubs removed")+
      ylab("Proportion of species remained")+
      theme_light()

    ggplot(dat_degree[dat_degree$weighted=="unweighted",], aes(x=Number.hub.removed, y=remain.mean)) +
      geom_line()+
      geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
      xlab("Number of module hubs removed")+
      ylab("Proportion of species remained")+
      theme_light()

  }

  if(bet_centrality){

    model <- property %>%
      as.data.frame() %>%
      filter(!is.na(Betweeness_Centrality)) %>%
      arrange(desc(Betweeness_Centrality))
    tem <- round(length(model$Betweeness_Centrality)*0.05,0)
    module.hub = row.names(model)[1:tem]

    Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat_Betweeness_Centrality<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                                          weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

    write.csv(dat_Betweeness_Centrality,paste("./",tag,"_targeted_deletion_results_Betweeness_Centrality.csv"))

  }

  if(clo_centrality){

    model <- property %>%
      as.data.frame() %>%
      filter(!is.na(Closeness_Centrality)) %>%
      arrange(desc(Closeness_Centrality))
    tem <- round(length(model$Closeness_Centrality)*0.05,0)
    module.hub = row.names(model)[1:tem]

    Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat_Closeness_Centrality<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                                         weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

    write.csv(dat_Closeness_Centrality,paste("./",tag,"_targeted_deletion_results_Closeness_Centrality.csv"))

  }

  if(eig_centrality){

    model <- property %>%
      as.data.frame() %>%
      filter(!is.na(Eigenvector_Centrality)) %>%
      arrange(desc(Eigenvector_Centrality))
    tem <- round(length(model$Eigenvector_Centrality)*0.05,0)
    module.hub = row.names(model)[1:tem]

    Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat_Eigenvector_Centrality<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                                           weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

    write.csv(dat_Eigenvector_Centrality,paste("./",tag,"_targeted_deletion_results_Eigenvector_Centrality.csv"))

  }

  if(pagerank){

    model <- property %>%
      as.data.frame() %>%
      filter(!is.na(PageRank)) %>%
      arrange(desc(PageRank))
    tem <- round(length(model$PageRank)*0.05,0)
    module.hub = row.names(model)[1:tem]

    Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat_PageRank<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                             weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

    write.csv(dat_PageRank,paste("./",tag,"_targeted_deletion_results_PageRank.csv"))

  }

  if(kcore){

    model <- property %>%
      as.data.frame() %>%
      filter(!is.na(K_core)) %>%
      arrange(desc(K_core))
    tem <- round(length(model$K_core)*0.05,0)
    module.hub = row.names(model)[1:tem]

    Weighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub),keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=T,nperm=100)
    Unweighted.simu<-rmsimu(netRaw=cor1_raw, rm.p.list=1:length(module.hub), keystonelist=module.hub, sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

    dat_K_core<-data.frame(Number.hub.removed=rep(1:length(module.hub),2),rbind(Weighted.simu,Unweighted.simu),
                           weighted=rep(c("weighted","unweighted"),each=length(module.hub)))

    write.csv(dat_K_core,paste("./",tag,"_targeted_deletion_results_K_core.csv"))

  }


}
