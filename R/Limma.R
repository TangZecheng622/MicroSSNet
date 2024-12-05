#' Differential Analysis of Nodes or Edges
#'
#' Perform differential analysis on network nodes or edges using the limma package with empirical Bayes.
#'
#' @param network Data frame containing network data, with samples as columns.
#' @param sample_metadata Data frame containing sample metadata with grouping information.
#' @param offset Logical; if TRUE, treats the first two columns of `network` as edge identifiers.
#' @param log_transform Logical; if TRUE, applies log transformation to the network data.
#' @param vscol Character; name of the column in `sample_metadata` used for grouping.
#' @param save Logical; if TRUE, saves the volcano plot as a PDF.
#' @param plot Logical; if TRUE, generates a volcano plot.
#' @return A list of differential expression results for each comparison and an optional volcano plot.
#' @importFrom limma lmFit contrasts.fit eBayes topTable makeContrasts
#' @importFrom stats as.formula model.matrix
#' @importFrom data.table fwrite
#' @importFrom utils write.table
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggsave
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export
Differential_nodes <- function(network, sample_metadata, offset = FALSE, log_transform = FALSE, vscol, save = TRUE, plot = TRUE) {
  # Prepare sample metadata
  sample_metadata <- sample_metadata[sample_metadata$sample %in% colnames(network), , drop = FALSE]
  rownames(sample_metadata) <- sample_metadata$sample
  sample_metadata$sample <- NULL
  sample_metadata[[vscol]] <- factor(sample_metadata[[vscol]], levels = unique(sample_metadata[[vscol]]))

  # Define design matrix for limma differential analysis
  design_formula <- stats::as.formula(paste("~ 0 +", vscol))
  design <- stats::model.matrix(design_formula, data = sample_metadata)
  vscolname <- levels(sample_metadata[[vscol]])
  colnames(design) <- vscolname

  # Set row names and type identifier
  if (offset) {
    rownames(network) <- paste(network[[1]], network[[2]], sep = "-")
    network <- network[, -c(1, 2)]  # Remove edge identifiers if offset is TRUE
    typename <- "Edges"
  } else {
    typename <- "Nodes"
  }

  # Log transformation if specified
  if (log_transform) {
    network <- log2(network + 1)
  }


  # Initialize list to store results of each comparison
  res_list <- list()

  # Iterate over all pairwise comparisons
  for (i in seq_along(vscolname)[-length(vscolname)]) {
    for (j in seq((i + 1), length(vscolname))) {

      # Define groups and create contrast matrix
      group1 <- vscolname[i]
      group2 <- vscolname[j]
      comparison_name <- paste0(group1, "_vs_", group2)
      cont.matrix <- limma::makeContrasts(contrasts = paste(group1, "-", group2), levels = design)

      # Perform differential analysis using limma
      fit <- limma::lmFit(network, design)
      fit <- limma::contrasts.fit(fit, cont.matrix)
      fit <- limma::eBayes(fit)

      # Store results
      res <- limma::topTable(fit, adjust.method = "BH", number = nrow(network))
      res_list[[comparison_name]] <- res
    }
  }

  # Initialize empty data frame for combined results
  total_res <- data.frame(name = character(), logFC = numeric(), adj.P.Val = numeric(), group = character())

  # Combine results from each comparison
  for (comp in names(res_list)) {
    res <- res_list[[comp]][, c("logFC", "adj.P.Val")]
    res$name <- rownames(res)
    row.names(res) <- NULL
    res$group <- comp
    total_res <- rbind(total_res, res)
  }

  # Add differential expression labels
  total_res <- total_res %>%
    dplyr::mutate(
      Regulation = factor(
        ifelse(adj.P.Val < 0.05 & abs(logFC) >= 1,
               ifelse(logFC >= 1, 'Up', 'Down'), 'NS'),
        levels = c("Up", "Down", "NS")
      )
    )

  if (save) {
    output_path <- paste0("./SSN/Limma_Result/", typename, "_Limma_Result.tsv")
    if (!dir.exists(dirname(output_path))) {
      dir.create(dirname(output_path), recursive = TRUE)
    }
    data.table::fwrite(total_res, output_path, sep = "\t")
  }

  # Plot volcano plots if requested
  if (plot) {
    pv <- PlotVolcano(total_res, typename)

    if (save) {
      output_path <- paste0("./SSN/Limma_Result/", typename, "_Limma_Volcano.pdf")
      if (!dir.exists(dirname(output_path))) {
        dir.create(dirname(output_path), recursive = TRUE)
      }
      grDevices::pdf(output_path, width = 10, height = 8)
      print(pv)
      grDevices::dev.off()
    } else {
      print(pv)
    }
  }
  return(list(LimmaResult = total_res,LimmaPlot = plot))
}
