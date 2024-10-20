#'  鲁棒值脆弱值的可视化
#'
#'  @param rrob_sum_list 生成每组的鲁棒值的平均值
#'  @param rrob_detail_list 生成每组的鲁棒值的详细值
#'  @export
compare_rob_vul <- function(rrob_sum_list,rrob_detail_list){
  if(length(rrob_sum_list)|length(rrob_detail_list)==0){
    return(invisible(NULL))
  }
  rrob_sum <- data.frame()
  for (i in rrob_sum_list){

    rrob_sum <- rbind(rrob_sum,as.data.frame(i))

  }
  rrob_detail <- data.frame()
  for (i in rrob_detail_list){

    rrob_detail <- rbind(rrob_detail,as.data.frame(i))
  }
  rrob_detail <- subset(rrob_detail,Proportion.removed != 1)


  ###网络鲁棒性模拟

  p_rrob_sum = ggplot(rrob_sum, aes(x = Proportion.removed, y = remain.mean, group = interaction(weighted, group), color = group)) +
    geom_line() +
    geom_pointrange(aes(ymin = remain.mean - remain.sd, ymax = remain.mean + remain.sd), size = 0.2) +
    facet_wrap(~ weighted, ncol = 2) +
    xlab("Proportion of species removed") +
    ylab("Proportion of species remained") +
    theme_bw() +
    #    scale_color_manual(values = colr)+  # 这里可以根据组的数量调整颜色+
    ggtitle("The impact of random removal of SGB on network stability")
  p_rrob_sum =p_rrob_sum +scale_color_aaas()
  ggsave(paste0(sel_group,"random_removal_Mean_SD.0.05-1.pdf"),plot = p_rrob_sum,height = 6,width = 12)

  #  compare_groups=list(c('KM','DQ'),c('DQ','NM'),c('KM','NM'))
  compare_groups <- lapply(combn(group_list, 2, simplify = FALSE), function(x) as.character(x))
  rrob_detail$group=factor(rrob_detail$group,c(group_list))
  unique(rrob_detail$Proportion.removed)


  prrob_detail=ggplot(rrob_detail,aes(x=group,y=Proportion.remain,col=group,fill=group))+
    facet_wrap(~ Proportion.removed,nrow = 1)+
    geom_boxplot(
      width = .2, fill = "white",
      size = 1, outlier.shape = NA
    ) +
    ggdist::stat_halfeye(
      adjust = .33, ## bandwidth
      width = .67,
      color = NA, ## remove slab interval
      position = position_nudge(x = .2)
    ) +
    # gghalves::geom_half_point(
    #   side = "l",
    #   range_scale = .3,
    #   alpha = .5, size = 2
    # )+
    stat_compare_means(color="gray",comparisons=compare_groups,
                       label = "p.signif")+
    # scale_color_manual(values = colr)+
    # scale_fill_manual(values = colr)+
    #facet_wrap(~Index,scales = 'free',nrow = 1)+
    theme_bw()+
    ylab("Robustness") +
    labs(title = "Robustness_Random 5%~95% of SGBs's network ")

  prrob_detail <- prrob_detail+scale_color_lancet()
  prrob_detail <- prrob_detail+scale_fill_lancet()
  ggsave(paste0(sel_group,"random_removal_SGB_robustness.6x6.pdf"),plot = prrob_detail,height = 6,width = 12)

  pvul_detail=ggplot(rrob_detail,aes(x=group,y=Vulnerability,col=group,fill=group))+
    facet_wrap(~ Proportion.removed,nrow = 1)+
    geom_boxplot(
      width = .2, fill = "white",
      size = 1, outlier.shape = NA
    ) +
    ggdist::stat_halfeye(
      adjust = .33, ## bandwidth
      width = .67,
      color = NA, ## remove slab interval
      position = position_nudge(x = .2)
    ) +
    # gghalves::geom_half_point(
    #   side = "l",
    #   range_scale = .3,
    #   alpha = .5, size = 2
    # )+
    stat_compare_means(color="gray",comparisons=compare_groups,
                       label = "p.signif")+
    # scale_color_manual(values = colr)+
    # scale_fill_manual(values = colr)+
    #facet_wrap(~Index,scales = 'free',nrow = 1)+
    theme_bw()+
    ylab("Vulnerability") +
    labs(title = "Vulnerability_Random 5%~95% of SGBs's network ")

  pvul_detail <- pvul_detail+scale_color_lancet()
  pvul_detail <- pvul_detail+scale_fill_lancet()
  ggsave(paste0(sel_group,"random_removal_SGB_Vulnerability.6x6.pdf"),plot = pvul_detail,height = 6,width = 12)


}
