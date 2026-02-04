#' @keywords internal
coerce_to_table_and_group <- function(
    table1,
    group_df = NULL,
    phyloseq_tax_rank = NULL,
    phyloseq_transform = c("none", "relative")
) {
  phyloseq_transform <- match.arg(phyloseq_transform)

  if (methods::is(table1, "phyloseq")) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("Error: 'phyloseq' package is required to use phyloseq input. Please install it first.")
    }
    ps <- table1

    # ---- tax_glom strict checks ----
    if (!is.null(phyloseq_tax_rank)) {
      tt <- tryCatch(phyloseq::tax_table(ps), error = function(e) NULL)
      if (is.null(tt) || ncol(tt) == 0) {
        stop(
          sprintf("Error: phyloseq_tax_rank='%s' was provided, but tax_table(ps) is missing/empty. ",
                  phyloseq_tax_rank),
          "Please add taxonomy (tax_table) to the phyloseq object or set phyloseq_tax_rank=NULL."
        )
      }
      tt_df <- as.data.frame(tt)
      if (!(phyloseq_tax_rank %in% colnames(tt_df))) {
        stop(
          sprintf(
            "Error: phyloseq_tax_rank='%s' not found in tax_table(ps). Available ranks are: %s",
            phyloseq_tax_rank,
            paste(colnames(tt_df), collapse = ", ")
          )
        )
      }
      ps <- phyloseq::tax_glom(ps, taxrank = phyloseq_tax_rank, NArm = FALSE)
    }

    # ---- extract OTU table, ensure taxa x samples ----
    otu <- phyloseq::otu_table(ps)
    tab <- if (phyloseq::taxa_are_rows(otu)) as.data.frame(otu) else as.data.frame(t(otu))

    # preserve rownames while coercing numeric
    rn <- rownames(tab)
    tab <- as.data.frame(lapply(tab, function(x) as.numeric(as.character(x))))
    rownames(tab) <- rn

    # ---- group_df from sample_data if missing ----
    if (is.null(group_df)) {
      sd <- tryCatch(phyloseq::sample_data(ps), error = function(e) NULL)
      if (is.null(sd) || nrow(as.data.frame(sd)) == 0) {
        stop(
          "Error: phyloseq input provided but sample_data(ps) is missing/empty, and group_df is NULL. ",
          "Please either provide group_df explicitly or add sample_data to the phyloseq object."
        )
      }
      group_df <- as.data.frame(sd, stringsAsFactors = FALSE)
      class(group_df) <- "data.frame"


      # guarantee a 'sample' column (needed by process_group_and_table)
      if (!any(grepl("sample", colnames(group_df), ignore.case = TRUE))) {
        group_df$sample <- phyloseq::sample_names(ps)
      }
      rownames(group_df) <- group_df$sample
    } else {
      # if user provided group_df, enforce it has a sample-like column
      if (!any(grepl("sample", colnames(group_df), ignore.case = TRUE))) {
        stop("Error: group_df was provided but no column name contains 'sample'. Please include a sample column (e.g., 'sample').")
      }
    }

    # ---- optional transform ----
    if (phyloseq_transform == "relative") {
      cs <- colSums(tab)
      cs[cs == 0] <- 1
      tab <- sweep(tab, 2, cs, "/")
    }

    return(list(table1 = tab, group_df = group_df))
  }

  if (!is.data.frame(table1)) {
    stop("Error: 'table1' must be a data.frame abundance table or a phyloseq object.")
  }

  # also enforce group_df sample column if provided (optional but nice)
  if (!is.null(group_df) && !any(grepl("sample", colnames(group_df), ignore.case = TRUE))) {
    stop("Error: group_df was provided but no column name contains 'sample'. Please include a sample column (e.g., 'sample').")
  }

  return(list(table1 = table1, group_df = group_df))
}
