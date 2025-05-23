#' PCA and PCoA Plotting Function
#'
#' This function generates 2D PCA or PCoA plots for structured data, allowing for centroid and ellipse visualizations.
#'
#' @param table Data frame. The input data where columns represent samples.
#' @param group_df Data frame. Group metadata with sample names and a grouping column.
#' @param vscol Character. Column name in `group_df` indicating the group variable.
#' @param save Logical. If TRUE, saves the plot as PDF.
#' @param showType Character. Options for additional plot features: `"centroid"` or `"ellipse"`. Default is NULL.
#' @param mode Character. Analysis mode: `"PCA"` or `"PCOA"`. Default is `"PCOA"`.
#' @param offset Logical. If TRUE, combines two columns into "edges" in the input table. Default is FALSE.
#' @param method Methods for calculating distance."manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", or "robust.aitchison".
#' @importFrom vegan adonis2 vegdist
#' @importFrom ggplot2 ggplot aes geom_point geom_segment geom_label labs theme_bw theme element_text stat_ellipse
#' @importFrom dplyr summarise group_by left_join
#' @importFrom ggsci scale_color_aaas scale_fill_aaas
#' @importFrom tidyr unite
#' @importFrom stats prcomp dist cmdscale as.formula var
#' @importFrom grDevices pdf dev.off
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @export
#'
PCA_draw <- function(table, group_df, vscol, save = TRUE, showType = NULL, mode = "PCOA", offset = FALSE,method = "euclidean") {
  # Load required packages
  required_packages <- c("vegan", "ggplot2", "dplyr", "ggsci", "tidyr")
  lapply(required_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed.")
      )}
  })

  # Data preparation
  if (offset) {
    table <- table %>%
      tidyr::unite("edges", 1, 2, sep = "_")
    row_names <- as.character(table$edges)
    table <- table[, -1]
    rownames(table) <- row_names
    df_t <- as.matrix(t(table))
#    df_t <- df_t[rowSums(df_t>0)>0,]
    # df_t <- df_t[, apply(df_t, 2, var) != 0]
    filename <- "edges_weight"
  } else {
    df_t <- as.matrix(t(table))
    # df_t <- df_t[, apply(df_t, 2, var) != 0]
#    df_t <- df_t[rowSums(df_t>0)>0,]
    filename <- "nodes_weight"
  }

  df_t <- df_t[, colSums(abs(df_t) > 0) > 0]

  group_df <- as.data.frame(group_df)
  group_df[[vscol]] <- factor(group_df[[vscol]])
  rownames(group_df) <- group_df$sample

  if (mode == "PCA") {
    col_variance <- apply(df_t, 2, var)
    zero_variance_cols <- which(col_variance == 0)

    if (length(zero_variance_cols) > 0) {
      df_t <- df_t[, -zero_variance_cols]
      cat("Removed columns with zero variance:", names(df_t)[zero_variance_cols], "\n")
    }

    pca_res <- prcomp(df_t, scale. = TRUE, center = TRUE)
    pca_data <- as.data.frame(pca_res$x)
    pca_data$sample <- rownames(pca_data)
    pca_data <- merge(pca_data, group_df, by = "sample")

    pca_var <- pca_res$sdev^2
    pca_var_perc <- round(pca_var / sum(pca_var) * 100, 2)

    groupVS_String <- paste(unique(pca_data[[vscol]]), collapse = "_Vs_")

    if (is.null(showType)) showType <- "default"

    if (showType == "centroid") {
      pca_means <- dplyr::summarize(
        dplyr::group_by(pca_data, !!sym(vscol)),
        PC1_mean = mean(PC1, na.rm = TRUE),
        PC2_mean = mean(PC2, na.rm = TRUE)
      )
      pca_data <- left_join(pca_data, pca_means, by = vscol)

      plot_pca <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[vscol]])) +
        ggplot2::geom_segment(ggplot2::aes(xend = PC1_mean, yend = PC2_mean), alpha = 0.6, size = 1.2, show.legend = FALSE) +
        ggplot2::geom_point(size = 1) +
        ggplot2::geom_label(data = pca_means, ggplot2::aes(x = PC1_mean, y = PC2_mean, label = .data[[vscol]]), size = 3, fontface = "bold") +
        ggplot2::labs(
          title = paste("PCA Analysis", groupVS_String),
          x = paste("PC1 (", pca_var_perc[1], "%)", sep = ""),
          y = paste("PC2 (", pca_var_perc[2], "%)", sep = "")
        ) +
        ggplot2::theme_bw() +
        ggsci::scale_color_aaas() +
        ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (showType == "ellipse") {
      plot_pca <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[vscol]])) +
        ggplot2::geom_point(size = 1) +
        ggplot2::stat_ellipse(aes(fill = .data[[vscol]]), level = 0.70, geom = "polygon", alpha = 0.2, color = NA) +
        ggplot2::labs(
          title = paste("PCA Analysis", groupVS_String),
          x = paste("PC1 (", pca_var_perc[1], "%)", sep = ""),
          y = paste("PC2 (", pca_var_perc[2], "%)", sep = "")
        ) +
        ggplot2::theme_bw() +
        ggsci::scale_color_aaas() +
        ggsci::scale_fill_aaas() +
        ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      plot_pca <- ggplot2::ggplot(pca_data, ggplot2::aes(x = PC1, y = PC2, color = .data[[vscol]])) +
        ggplot2::geom_point(size = 1) +
        ggplot2::labs(
          title = paste("PCA Analysis", groupVS_String),
          x = paste("PC1 (", pca_var_perc[1], "%)", sep = ""),
          y = paste("PC2 (", pca_var_perc[2], "%)", sep = "")
        ) +
        ggplot2::theme_bw() +
        ggsci::scale_color_aaas() +
        ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
    }

    if (save) {
      output_path <- paste0("./SSN/PCA_Result/", groupVS_String, "_", filename, "_PCA_Plot_2d.pdf")
      output_file <- paste0("./SSN/PCA_Result/", groupVS_String, "_", filename, "_PCA.csv")
      if (!dir.exists(dirname(output_path))) {
        dir.create(dirname(output_path), recursive = TRUE)
      }
      write.csv(pca_data,output_file)
      pdf(output_path, width = 6, height = 6)
      print(plot_pca)
      dev.off()
    } else {
      print(plot_pca)
    }
  } else if (mode == "PCOA") {
    dist_matrix <- vegan::vegdist(df_t,method = method)
    # dist_matrix <- vegan::vegdist(df_t,method = "bray")
    # dist_matrix <- dist(df_t)
    pcoa_res <- cmdscale(dist_matrix, k = 2,eig = TRUE)
    pcoa_points <- pcoa_res$points
    group_df <- group_df[match(rownames(pcoa_points), group_df$sample), ]
    eig_vals <- pcoa_res$eig
    total_variance <- sum(eig_vals)
    explained_variance <- round(eig_vals / total_variance * 100 ,2) # 转化为百分比

    pcoa_data <- data.frame(PCoA1 = pcoa_points[, 1], PCoA2 = pcoa_points[, 2], vscol = group_df[[vscol]])
    colnames(pcoa_data)[3] <- vscol

    formula_dynamic <- as.formula(paste("dist_matrix ~", vscol))
    adonis_res <- vegan::adonis2(formula_dynamic, data = group_df)
    p_val <- adonis_res$`Pr(>F)`[1]
    r_squared <- adonis_res$R2[1]

    groupVS_String <- paste(unique(pcoa_data[[vscol]]), collapse = "_Vs_")

    if (is.null(showType)) showType <- "default"


    if (showType == "centroid") {
      pcoa_means <- dplyr::summarize(
        dplyr::group_by(pcoa_data, !!sym(vscol)),
        PCoA1_mean = mean(PCoA1, na.rm = TRUE),
        PCoA2_mean = mean(PCoA2, na.rm = TRUE)
      )
      pcoa_data <- left_join(pcoa_data, pcoa_means, by = vscol)

      plot_pcoa <- ggplot2::ggplot(pcoa_data, ggplot2::aes(x = PCoA1, y = PCoA2, color = .data[[vscol]])) +
        ggplot2::geom_segment(ggplot2::aes(xend = PCoA1_mean, yend = PCoA2_mean), alpha = 0.6, size = 1.2, show.legend = FALSE) +
        ggplot2::geom_point(size = 1) +
        ggplot2::geom_label(data = pcoa_means, ggplot2::aes(x = PCoA1_mean, y = PCoA2_mean, label = .data[[vscol]]), size = 3, fontface = "bold") +
        ggplot2::labs(
          title = sprintf("PCoA Analysis %s\nPERMANOVA: P-value = %.4f,  R^2 = %.4f", groupVS_String, p_val, r_squared),
          x = paste("PCoA1 (", explained_variance[1], "%)", sep = ""),
          y = paste("PCoA2 (", explained_variance[2], "%)", sep = "")
        ) +
        ggplot2::theme_bw() +
        ggsci::scale_color_aaas() +
        ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5))
    } else if (showType == "ellipse") {
      plot_pcoa <- ggplot2::ggplot(pcoa_data, ggplot2::aes(x = PCoA1, y = PCoA2, color = .data[[vscol]])) +
        ggplot2::geom_point(size = 1) +
        ggplot2::stat_ellipse(aes(fill = .data[[vscol]]), level = 0.70, geom = "polygon", alpha = 0.2, color = NA) +
        ggplot2::labs(
          title = sprintf("PCoA Analysis %s\nPERMANOVA: P-value = %.4f,  R^2 = %.4f", groupVS_String, p_val, r_squared),
          x = paste("PCoA1 (", explained_variance[1], "%)", sep = ""),
          y = paste("PCoA2 (", explained_variance[2], "%)", sep = "")
        ) +
        ggplot2::theme_bw() +
        ggsci::scale_color_aaas() +
        ggsci::scale_fill_aaas() +
        ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
    } else {
      plot_pcoa <- ggplot2::ggplot(pcoa_data, ggplot2::aes(x = PCoA1, y = PCoA2, color = .data[[vscol]])) +
        ggplot2::geom_point(size = 1) +
        ggplot2::labs(
          title = sprintf("PCoA Analysis %s\nPERMANOVA: P-value = %.4f,  R^2 = %.4f", groupVS_String, p_val, r_squared),
          x = paste("PCoA1 (", explained_variance[1], "%)", sep = ""),
          y = paste("PCoA2 (", explained_variance[2], "%)", sep = "")
        ) +
        ggplot2::theme_bw() +
        ggsci::scale_color_aaas() +
        ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5))
    }

    if (save) {
      output_path <- paste0("./SSN/PCoA_Result/", groupVS_String, "_", filename, "_PCoA_Plot_2d.pdf")
      output_file <- paste0("./SSN/PCoA_Result/", groupVS_String, "_", filename, "_PCoA.csv")
      if (!dir.exists(dirname(output_path))) {
        dir.create(dirname(output_path), recursive = TRUE)
      }
      write.csv(pcoa_data,output_file)
      pdf(output_path, width = 6, height = 6)
      print(plot_pcoa)
      dev.off()
    } else {
      print(plot_pcoa)
    }
  } else {
    stop("Invalid mode. Choose either 'PCA' or 'PCOA'.")
  }
}
