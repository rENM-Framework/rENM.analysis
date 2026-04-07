#' Save a suitability trend map for a species (via plot_trend)
#'
#' Generates a climatic suitability trend map using plot_trend() and saves
#' the result as a PNG file within the standard rENM project directory
#' structure for the specified species.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' The function accepts either a file path to a raster or an in-memory
#' raster object containing trend values. File-based inputs are expected
#' to follow the naming pattern:
#'
#' \code{<alpha_code>-<NAME>.tif}
#' or
#' \code{<alpha_code>-<NAME>.asc}
#'
#' where <NAME> represents a descriptive identifier for the trend surface.
#'
#' \strong{Outputs}:
#' The PNG file is written to:
#'
#' \code{<rENM_project_dir()>/runs/<alpha_code>/}
#' \code{Trends/suitability/}
#'
#' using the filename:
#'
#' \code{<alpha_code>-<NAME>.png}
#'
#' where <NAME> is derived from the input raster filename. For example,
#' an input file named:
#'
#' CASP-Suitability-Difference-Trend.tif
#'
#' will produce:
#'
#' CASP-Suitability-Difference-Trend.png
#'
#' If the function cannot confidently extract <NAME> from the input path
#' (for example, when a raster object is supplied instead of a file path),
#' it falls back to:
#'
#' \code{<alpha_code>-Suitability-Trend.png}
#'
#' \strong{Plot generation}:
#' The function calls \code{plot_trend()} to generate a ggplot object. If
#' \code{plot_trend()} returns a list, the \code{plot} element is used.
#'
#' \strong{Styling}:
#' The saved PNG uses a white background and consistent sizing
#' (1000 x 750 pixels at 72 dpi) with controlled font sizes for titles,
#' legends, and axes.
#'
#' \strong{Data requirements}:
#' The input raster must represent a suitability trend surface compatible
#' with \code{plot_trend()}.
#'
#' @param alpha_code Character. Four-letter banding code for the species.
#' @param raster_file Character, SpatRaster, or Raster. Path to a raster
#'   file (.tif or .asc) or an in-memory raster object containing trend
#'   values.
#'
#' @return Invisibly returns a Character scalar representing the file path
#' to the saved PNG output.
#'
#' Side effects:
#' \itemize{
#'   \item Writes a PNG file to the Trends/suitability directory
#'   \item Creates the output directory if it does not exist
#'   \item Writes progress messages to the console
#' }
#'
#' @importFrom ggplot2 ggsave theme element_text
#' @importFrom utils tail
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' \dontrun{
#'   save_trend_plot(
#'     alpha_code  = "CASP",
#'     raster_file = file.path(
#'       rENM_project_dir(),
#'       "runs", "CASP", "Trends", "suitability",
#'       "CASP-Suitability-Trend.tif"
#'     )
#'   )
#' }
#'
#' @seealso [plot_trend()], [find_suitability_trend()]
#'
#' @export
save_trend_plot <- function(alpha_code, raster_file) {

  if (!is.character(alpha_code) || length(alpha_code) != 1L || nchar(alpha_code) != 4L) {
    stop("`alpha_code` must be a single four-letter string (for example, 'CASP').",
         call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.", call. = FALSE)
  }
  if (!exists("plot_trend", mode = "function")) {
    stop("Internal helper `plot_trend()` not found. Is rENM.core loaded?", call. = FALSE)
  }

  code <- toupper(alpha_code)

  ## ------------------------- Generate plot -------------------------------- ##
  res <- plot_trend(code, raster_file)
  gp  <- if (is.list(res) && "plot" %in% names(res)) res$plot else res

  if (!inherits(gp, "ggplot")) {
    stop("plot_trend() did not return a ggplot object (or list(plot = ...)).",
         call. = FALSE)
  }

  ## -------------------------- Styling ------------------------------------- ##
  gp <- gp +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(size = 28, face = "bold"),
      legend.title = ggplot2::element_text(size = 15),
      legend.text  = ggplot2::element_text(size = 15),
      axis.title   = ggplot2::element_text(size = 20, face = "plain"),
      axis.text    = ggplot2::element_text(size = 18)
    )

  ## -------------------------- Paths --------------------------------------- ##
  root_dir <- file.path(rENM_project_dir(), "runs", code)
  out_dir  <- file.path(root_dir, "Trends", "suitability")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  ## -------------------------- Name extraction ----------------------------- ##
  name_hint <- "Suitability-Trend"

  if (is.character(raster_file) && length(raster_file) == 1L) {
    fbase <- basename(raster_file)

    stem <- tryCatch(
      utils::tail(strsplit(fbase, "\\.", fixed = FALSE)[[1]], 2L)[1],
      error = function(e) tools::file_path_sans_ext(fbase)
    )

    if (is.na(stem) || !nzchar(stem)) {
      stem <- tools::file_path_sans_ext(fbase)
    }

    if (grepl(paste0("^", code, "-"), stem, ignore.case = TRUE)) {
      name_part <- sub(paste0("^", code, "-"), "", stem, ignore.case = TRUE)
      if (nzchar(name_part)) name_hint <- name_part
    } else {
      parts <- strsplit(stem, "-", fixed = TRUE)[[1]]
      if (length(parts) >= 2L) {
        name_hint <- paste(parts[-1L], collapse = "-")
      } else if (nzchar(stem)) {
        name_hint <- stem
      }
    }

    name_hint <- gsub("[^A-Za-z0-9._-]+", "_", name_hint)
    if (!nzchar(name_hint)) name_hint <- "Suitability-Trend"
  }

  ## -------------------------- Output -------------------------------------- ##
  out_file <- file.path(out_dir, sprintf("%s-%s.png", code, name_hint))

  ggplot2::ggsave(
    filename  = out_file,
    plot      = gp,
    width     = 1000,
    height    = 750,
    units     = "px",
    dpi       = 72,
    pointsize = 15,
    bg        = "white"
  )

  ## -------------------------- Console ------------------------------------- ##
  cat("------------------------------------------------------------------------\n")
  cat(sprintf("save_trend_plot(): Saved plot for %s\n", code))
  cat(sprintf("Derived name: %s\n", name_hint))
  cat(sprintf("Output file: %s\n", out_file))

  invisible(out_file)
}
