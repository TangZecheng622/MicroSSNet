#' Compare Robustness and Vulnerability Across Groups
#'
#' This function generates plots comparing network robustness, vulnerability, and complexity across multiple groups.
#'
#' @param rrob_sum_list A list of data frames containing robustness summary for each group.
#' @param rrob_detail_list A list of data frames containing robustness details for each group.
#' @param group_list A character vector of group names.
#' @param sel_group A character string for output file labeling.
#' @param output_dir Character string specifying the output directory.
#' @importFrom ggplot2 ggsave
#' @importFrom grDevices dev.off
#' @importFrom utils combn
#' @export
compare_rob_vul <- function(rrob_sum_list, rrob_detail_list, group_list, sel_group = "Group", output_dir = "./") {
  if (length(rrob_sum_list) == 0 || length(rrob_detail_list) == 0) {
    message("Warning: Input lists are empty.")
    return(invisible(NULL))
  }

  # Combine data
  rrob_sum <- do.call(rbind, rrob_sum_list)
  rrob_detail <- do.call(rbind, rrob_detail_list)

  # Prepare output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  compare_groups <- if (length(group_list) > 1 & length(group_list) <= 3 & length(rrob_sum_list) > 1 & length(rrob_sum_list) <= 3) {
    lapply(utils::combn(group_list, 2, simplify = FALSE), function(x) as.character(x))
  } else {
    print("If there is only one group, or if there are more than three groups, no significant comparison will be made.")
    NULL  # 只有一个组时，设为 NULL
  }

  rrob_detail$group <- factor(rrob_detail$group, levels = group_list)
  write.csv(x = rrob_detail,file = paste0(sel_group,"_Random_Removal_Detail.csv"))
  write.csv(x = rrob_sum,file = paste0(sel_group,"_Random_Removal_Summary.csv"))
  # Plot robustness summary
  p_rrob_sum <- plot_robustness_summary(rrob_sum)
  ggplot2::ggsave(filename = file.path(output_dir, paste0(sel_group, "_random_removal_Mean_SD.pdf")),
                  plot = p_rrob_sum, height = 6, width = 12)

  # Plot robustness detail
  prrob_detail <- plot_robustness_detail(rrob_detail, compare_groups)
  ggplot2::ggsave(filename = file.path(output_dir, paste0(sel_group, "_random_removal_Robustness.pdf")),
                  plot = prrob_detail, height = 6, width = 6)

  # Plot vulnerability detail
  pvul_detail <- plot_vulnerability_detail(rrob_detail, compare_groups)
  ggplot2::ggsave(filename = file.path(output_dir, paste0(sel_group, "_random_removal_Vulnerability.pdf")),
                  plot = pvul_detail, height = 6, width = 6)

  # Plot complexity detail
  cpx_detail <- plot_complexity_detail(rrob_detail, compare_groups)
  ggplot2::ggsave(filename = file.path(output_dir, paste0(sel_group, "_random_removal_Complexity.pdf")),
                  plot = cpx_detail, height = 6, width = 6)
}


#' Plot Robustness Summary
#'
#' @param rrob_sum Data frame containing robustness summary data.
#' @return A ggplot object.
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_wrap labs theme_minimal
#' @importFrom ggsci scale_color_aaas
#' @export
plot_robustness_summary <- function(rrob_sum) {
  # Reshape data for plotting
  df_long <- tidyr::pivot_longer(rrob_sum,
                                 cols = c("Robustness.mean", "Vulnerability.mean", "Complexity.mean"),
                                 names_to = "Metric",
                                 values_to = "Mean")

  ggplot2::ggplot(df_long, ggplot2::aes(x = Proportion.removed, y = Mean, color = group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~ Metric, scales = "free_y") +
    ggplot2::labs(x = "Proportion Removed", y = "Mean Value", title = "Robustness, Vulnerability, and Complexity Summary") +
    ggplot2::theme_minimal() +
    ggsci::scale_color_aaas()
}
#' Plot Robustness Detail
#'
#' This function generates a detailed plot of network robustness under random node removal.
#'
#' @param rrob_detail Data frame containing the detailed robustness data.
#' @param compare_groups List of group comparisons for statistical testing (optional).
#' @return A `ggplot` object representing the robustness plot.
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw ylab labs position_dodge
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggdist stat_halfeye
#' @importFrom ggsci scale_color_lancet scale_fill_lancet
#' @export
plot_robustness_detail <- function(rrob_detail, compare_groups = NULL) {
  if (!is.data.frame(rrob_detail)) stop("Error: 'rrob_detail' must be a data frame.")
  if (!is.null(compare_groups) && !is.list(compare_groups)) stop("Error: 'compare_groups' must be a list or NULL.")

  rrob_detail$Proportion.removed <- as.factor(rrob_detail$Proportion.removed)

  p <- ggplot2::ggplot(rrob_detail, ggplot2::aes(x = Proportion.removed, y = Robustness, color = group, fill = group)) +
    ggplot2::geom_boxplot(width = 0.25, fill = "white", size = 0.1, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.6)) +
    ggdist::stat_halfeye(adjust = 0.33, width = 0.67, color = NA, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Robustness") +
    ggplot2::labs(title = "Robustness of Network under Random Node Removal") +
    ggsci::scale_color_lancet() +
    ggsci::scale_fill_lancet()

  # Add statistical comparisons if compare_groups is provided
  if (!is.null(compare_groups)) {
    p <- p + ggpubr::stat_compare_means(
      # comparisons = compare_groups,
      ggplot2::aes(group = group),
      label = "p.signif",
      method = "wilcox.test",
      hide.ns = FALSE
    )
  }

  return(p)
}
#' Plot Vulnerability Detail
#'
#' This function generates a detailed plot of network vulnerability under random node removal.
#'
#' @param rrob_detail Data frame containing the detailed vulnerability data.
#' @param compare_groups List of group comparisons for statistical testing (optional).
#' @return A `ggplot` object representing the vulnerability plot.
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw ylab labs position_dodge
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggdist stat_halfeye
#' @importFrom ggsci scale_color_lancet scale_fill_lancet
#' @export
plot_vulnerability_detail <- function(rrob_detail, compare_groups = NULL) {
  if (!is.data.frame(rrob_detail)) stop("Error: 'rrob_detail' must be a data frame.")
  if (!is.null(compare_groups) && !is.list(compare_groups)) stop("Error: 'compare_groups' must be a list or NULL.")

  rrob_detail$Proportion.removed <- as.factor(rrob_detail$Proportion.removed)

  p <- ggplot2::ggplot(rrob_detail, ggplot2::aes(x = Proportion.removed, y = Vulnerability, color = group, fill = group)) +
    ggplot2::geom_boxplot(width = 0.25, fill = "white", size = 0.1, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.6)) +
    ggdist::stat_halfeye(adjust = 0.33, width = 0.67, color = NA, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Vulnerability") +
    ggplot2::labs(title = "Vulnerability of Network under Random Node Removal") +
    ggsci::scale_color_lancet() +
    ggsci::scale_fill_lancet()

  # Add statistical comparisons if compare_groups is provided
  if (!is.null(compare_groups)) {
    p <- p + ggpubr::stat_compare_means(
      # comparisons = compare_groups,
      ggplot2::aes(group = group),
      label = "p.signif",
      method = "wilcox.test",
      hide.ns = FALSE
    )
  }

  return(p)
}
#' Plot Complexity Detail
#'
#' This function generates a detailed plot of network complexity under random node removal.
#'
#' @param rrob_detail Data frame containing the detailed complexity data.
#' @param compare_groups List of group comparisons for statistical testing (optional).
#' @return A `ggplot` object representing the complexity plot.
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw ylab labs position_dodge
#' @importFrom ggpubr stat_compare_means
#' @importFrom ggdist stat_halfeye
#' @importFrom ggsci scale_color_lancet scale_fill_lancet
#' @export
plot_complexity_detail <- function(rrob_detail, compare_groups = NULL) {
  if (!is.data.frame(rrob_detail)) stop("Error: 'rrob_detail' must be a data frame.")
  if (!is.null(compare_groups) && !is.list(compare_groups)) stop("Error: 'compare_groups' must be a list or NULL.")

  rrob_detail$Proportion.removed <- as.factor(rrob_detail$Proportion.removed)

  p <- ggplot2::ggplot(rrob_detail, ggplot2::aes(x = Proportion.removed, y = Complexity, color = group, fill = group)) +
    ggplot2::geom_boxplot(width = 0.25, fill = "white", size = 0.1, outlier.shape = NA, position = ggplot2::position_dodge(width = 0.6)) +
    ggdist::stat_halfeye(adjust = 0.33, width = 0.67, color = NA, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Complexity") +
    ggplot2::labs(title = "Complexity of Network under Random Node Removal") +
    ggsci::scale_color_lancet() +
    ggsci::scale_fill_lancet()

  # Add statistical comparisons if compare_groups is provided
  if (!is.null(compare_groups)) {
    p <- p + ggpubr::stat_compare_means(
      # comparisons = compare_groups,
      ggplot2::aes(group = group),
      label = "p.signif",
      method = "wilcox.test",
      hide.ns = FALSE

    )
  }

  return(p)
}
