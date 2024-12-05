#' Volcano Plot Function
#'
#' Generate a volcano plot for visualizing differential expression across groups.
#'
#' @param total_res Data frame containing differential analysis results.
#' @param typename Character string indicating the type of data ("Nodes" or "Edges").
#' @return A ggplot2 object of the volcano plot.
#' @importFrom ggplot2 ggplot aes geom_point theme_bw labs theme element_text scale_color_manual geom_col geom_jitter guides
#' @importFrom dplyr group_by summarise
#' @importFrom ggsci scale_fill_npg
#' @importFrom scales alpha
#' @importFrom magrittr %>%
#' @export
PlotVolcano <- function(total_res, typename) {
  # Check if 'group' column exists
  if ("group" %in% colnames(total_res)) {
    num_groups <- length(unique(total_res$group))
  } else {
    num_groups <- 1
    total_res$group <- "Group"
  }

  # Plotting
  if (num_groups == 1) {
    p <- ggplot2::ggplot(total_res, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = Regulation)) +
      ggplot2::geom_point(alpha = 0.6, size = 1.5) +
      ggplot2::scale_color_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey")) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste("Volcano Plot of Differential", typename),
        x = "Log2 Fold Change",
        y = "-Log10 Adjusted P-value"
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5),
        legend.title = ggplot2::element_blank()
      )
  } else {
    # Calculate background bars for each group
    res_bg <- total_res %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        max_log2FC = max(logFC, na.rm = TRUE),
        min_log2FC = min(logFC, na.rm = TRUE)
      )

    # Plotting
    p <- ggplot2::ggplot() +
      ggplot2::geom_col(data = res_bg, ggplot2::aes(x = group, y = max_log2FC), fill = "grey85", width = 0.8, alpha = 0.5) +
      ggplot2::geom_col(data = res_bg, ggplot2::aes(x = group, y = min_log2FC), fill = "grey85", width = 0.8, alpha = 0.5) +
      ggplot2::geom_jitter(data = total_res, ggplot2::aes(x = group, y = logFC, color = Regulation), size = 3, width = 0.4, alpha = 0.7) +
      ggplot2::scale_color_manual(values = c("Up" = "#e42313", "Down" = "#0061d5", "NS" = "#8b8c8d")) +
      ggplot2::geom_col(data = res_bg, ggplot2::aes(x = group, y = 1,fill = group), width = 0.8) +
      ggplot2::geom_col(data = res_bg, ggplot2::aes(x = group, y = -1,fill = group), width = 0.8) +
      ggsci::scale_fill_npg() +
      ggplot2::geom_text(data = res_bg,ggplot2::aes(x = group, y = 0, label = group),size = 4,color = "#dbebfa")+
      ggplot2::theme_classic() +
      ggplot2::theme(
        axis.line.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_line(linewidth = 0.8),
        axis.text.y = ggplot2::element_text(size = 12, color = "black"),
        axis.title = ggplot2::element_text(size = 14, color = "black"),
        axis.ticks.y = ggplot2::element_line(linewidth = 0.8)
      ) +
      ggplot2::labs(x = "Group", y = "Log2 Fold Change", fill = NULL, color = NULL) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 6, alpha = 1)))
  }

  return(p)
}
