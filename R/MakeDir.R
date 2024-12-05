#' Create an Output Directory with Logging
#'
#' This function ensures the specified directory exists, creating it if necessary,
#' and provides logging to indicate its actions.
#'
#' @param output_dir A character string specifying the path of the directory to be created.
#' @param recursive Logical. Should parent directories be created if they don't exist? Default is `TRUE`.
#' @param verbose Logical. Should messages be printed to the console? Default is `FALSE`.
#' @return Logical. `TRUE` if the directory exists or was created successfully, otherwise `FALSE`.
#' @examples
#' # Create a directory named "output" in the current working directory
#' create_output_dir("output", verbose = TRUE)
#'
#' # Create a nested directory structure
#' create_output_dir("nested/directory/structure", verbose = TRUE)
#'
#' @export
create_output_dir <- function(output_dir, recursive = TRUE, verbose = FALSE) {
  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("`output_dir` must be a single character string specifying a valid directory path.")
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = recursive, showWarnings = FALSE)
    if (verbose) {
      message("Directory created: ", output_dir)
    }
  } else if (verbose) {
    message("Directory already exists: ", output_dir)
  }
  dir.exists(output_dir)
}
