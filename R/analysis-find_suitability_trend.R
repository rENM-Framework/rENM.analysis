#' Compute Theil-Sen suitability trend, p, and z across the time series
#'
#' Computes per-cell Theil-Sen trends from a time series of predicted
#' suitability rasters and derives Mann-Kendall statistics (p and Z)
#' to assess temporal trends in climatic suitability.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}:
#' Input rasters are expected at:
#'
#' \code{<rENM_project_dir()>/runs/<alpha_code>/}
#' \code{TimeSeries/<year>/model/}
#'
#' Input rasters are expected to be named:
#' \code{<alpha_code>-<year>-Prediction.asc}
#'
#' for years 1980, 1985, ..., 2020 (5-year increments).
#'
#' Only years whose Prediction.asc files exist are included in the
#' analysis. An error is raised if fewer than three valid years are
#' available.
#'
#' \strong{Outputs}:
#' Results are written to:
#'
#' \code{<rENM_project_dir()>/runs/<alpha_code>/}
#' \code{Trends/suitability/}
#'
#' as both AAIGrid (.asc) and GeoTIFF (.tif) with filenames:
#'
#' \itemize{
#'   \item <alpha_code>-Suitability-Trend.\{asc,tif\}
#'         (Theil-Sen slope per year)
#'   \item <alpha_code>-Suitability-Trend-p.\{asc,tif\}
#'         (Mann-Kendall p-value)
#'   \item <alpha_code>-Suitability-Trend-z.\{asc,tif\}
#'         (Mann-Kendall Z statistic)
#' }
#'
#' Any .xml or .prj sidecar files in the output directory are removed.
#' A formatted processing summary is appended to:
#'
#' \code{<rENM_project_dir()>/runs/<alpha_code>/}
#' \code{_log.txt}
#'
#' If available, the function will call
#' \code{save_trend_plot()} to generate a PNG visualization
#' of the trend raster from the GeoTIFF output.
#'
#' The function writes clear progress messages to the console.
#'
#' \strong{Trend estimation}:
#' The Theil-Sen slope (median slope) is computed per raster cell using
#' \code{mblm::mblm()}. This provides a robust estimate of temporal
#' change in suitability (units per year).
#'
#' \strong{Significance testing}:
#' Mann-Kendall trend statistics (Z, p, and Kendall's tau) are computed
#' internally using a dependency-free implementation.
#'
#' \strong{Data requirements}:
#' Cells with fewer than three non-NA time points are assigned NA in all
#' outputs.
#'
#' @param alpha_code Character. Four-letter banding code identifying the
#'   species (e.g., "CASP").
#'
#' @return Invisibly returns a List with the following structure:
#' \itemize{
#'   \item paths:
#'     \itemize{
#'       \item asc: List of file paths to AAIGrid outputs
#'       \item tif: List of file paths to GeoTIFF outputs
#'       \item png: Character path to PNG output (if generated)
#'       \item log: Character path to log file
#'     }
#'   \item n_years:
#'     Integer number of years used in the analysis
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes raster outputs to the Trends/suitability directory
#'   \item Appends a processing summary to the log file
#'   \item Optionally generates a PNG visualization
#' }
#'
#' @importFrom terra rast writeRaster app ncell values
#' @importFrom mblm mblm
#' @importFrom stats coef median pnorm
#'
#' @examples
#' \dontrun{
#'   find_suitability_trend("CASP")
#' }
#'
#' @seealso
#' [save_trend_plot()], [plot_trend()]
#'
#' @export
find_suitability_trend <- function(alpha_code) {

  ## --------------------------- Input validation --------------------------- ##
  if (!is.character(alpha_code) || length(alpha_code) != 1L || nchar(alpha_code) != 4L) {
    stop("`alpha_code` must be a single four-letter string.", call. = FALSE)
  }
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required.", call. = FALSE)
  }
  if (!requireNamespace("mblm", quietly = TRUE)) {
    stop("Package 'mblm' is required.", call. = FALSE)
  }

  code  <- toupper(alpha_code)
  years <- seq(1980, 2020, by = 5L)

  # start timer
  t_start <- Sys.time()

  ## ------------------------------- Paths ---------------------------------- ##
  root_dir <- file.path(rENM_project_dir(), "runs", code)

  ts_path <- function(yr) file.path(root_dir, "TimeSeries", yr, "model")
  in_file <- function(yr) file.path(ts_path(yr), sprintf("%s-%d-Prediction.asc", code, yr))

  trend_dir <- file.path(root_dir, "Trends", "suitability")
  if (!dir.exists(trend_dir)) dir.create(trend_dir, recursive = TRUE, showWarnings = FALSE)

  log_path <- file.path(root_dir, "_log.txt")

  cat("------------------------------------------------------------------------\n")
  cat(sprintf("find_suitability_trend(): starting for %s\n", code))

  ## ------------------------ Discover available years ---------------------- ##
  files <- vapply(
    years,
    function(yr) {
      f <- in_file(yr)
      if (file.exists(f)) f else NA_character_
    },
    FUN.VALUE = character(1)
  )

  ok_years <- years[!is.na(files)]
  files    <- files[!is.na(files)]

  if (length(ok_years) < 3L) {
    stop("Need at least 3 available years.", call. = FALSE)
  }

  years <- ok_years

  ## ---------------------------- Load rasters ------------------------------ ##
  r_stack <- terra::rast(files)

  ## -------------------- Per-cell Theil-Sen & MK -------------------------- ##
  theilsen_mk <- function(vals) {
    idx <- which(!is.na(vals))
    n   <- length(idx)
    if (n < 3L) return(c(NA_real_, NA_real_, NA_real_, NA_real_, n))

    y <- vals[idx]
    x <- years[idx]

    slope <- NA_real_
    fit <- try(mblm::mblm(y ~ x, repeated = FALSE), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      cf <- try(stats::coef(fit), silent = TRUE)
      if (!inherits(cf, "try-error") && length(cf) >= 2L) {
        slope <- unname(cf[2L])
      }
    }

    S <- 0
    for (k in 1:(n - 1)) {
      S <- S + sum(sign(y[(k + 1):n] - y[k]))
    }

    varS <- (n * (n - 1) * (2 * n + 5)) / 18

    if (S > 0) {
      Z <- (S - 1) / sqrt(varS)
    } else if (S < 0) {
      Z <- (S + 1) / sqrt(varS)
    } else {
      Z <- 0
    }

    p <- 2 * (1 - stats::pnorm(abs(Z)))
    tau <- S / (0.5 * n * (n - 1))

    c(slope, Z, p, tau, n)
  }

  stats_r <- terra::app(r_stack, fun = theilsen_mk)
  names(stats_r) <- c("slope", "Z", "p", "tau", "n")

  slope_r <- stats_r[["slope"]]
  z_r     <- stats_r[["Z"]]
  p_r     <- stats_r[["p"]]

  ## ------------------------------- Write outputs -------------------------- ##
  asc_out <- list(
    trend = file.path(trend_dir, sprintf("%s-Suitability-Trend.asc", code)),
    p     = file.path(trend_dir, sprintf("%s-Suitability-Trend-p.asc", code)),
    z     = file.path(trend_dir, sprintf("%s-Suitability-Trend-z.asc", code))
  )

  tif_out <- list(
    trend = file.path(trend_dir, sprintf("%s-Suitability-Trend.tif", code)),
    p     = file.path(trend_dir, sprintf("%s-Suitability-Trend-p.tif", code)),
    z     = file.path(trend_dir, sprintf("%s-Suitability-Trend-z.tif", code))
  )

  terra::writeRaster(slope_r, tif_out$trend, overwrite = TRUE)
  terra::writeRaster(p_r,     tif_out$p,     overwrite = TRUE)
  terra::writeRaster(z_r,     tif_out$z,     overwrite = TRUE)

  ## -------------------------- PNG output --------------------------------- ##
  png_path <- NA_character_

  if (exists("save_trend_plot", where = asNamespace("rENM.analysis"), mode = "function")) {

    cat("Calling save_trend_plot() to generate PNG...\n")

    png_path <- tryCatch(
      {
        rENM.analysis::save_trend_plot(code, tif_out$trend)
      },
      error = function(e) {
        warning(
          sprintf("save_trend_plot() failed: %s", conditionMessage(e)),
          call. = FALSE
        )
        NA_character_
      }
    )

  } else {
    cat("save_trend_plot() not found in namespace; skipping PNG generation.\n")
  }

  invisible(list(
    paths   = list(asc = asc_out, tif = tif_out, png = png_path, log = log_path),
    n_years = length(years)
  ))

  # ## ------------------------------ Logging -------------------------------- ##
  # log_block <- c(
  #   "",
  #   "------------------------------------------------------------------------",
  #   "Processing summary (find_suitability_trend)",
  #   sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  #   sprintf("Alpha code: %s", code),
  #   sprintf("Years used: %s", paste(years, collapse = ", "))
  # )
  #
  # cat(log_block, file = log_path, sep = "\n", append = TRUE)

  ## ------------------------------ Logging -------------------------------- ##

  # --- build summary stats ---
  n_years <- length(years)
  n_cells <- terra::ncell(r_stack)
  n_na_slope <- sum(is.na(terra::values(slope_r)))
  n_valid_slope <- n_cells - n_na_slope

  # --- elapsed time ---
  t_end <- Sys.time()
  elapsed <- round(as.numeric(difftime(t_end, t_start, units = "secs")), 2)

  log_block <- c(
    "------------------------------------------------------------------------",
    "Processing summary (find_suitability_trend)",
    sprintf("Start time:   %s", format(t_start, "%Y-%m-%d %H:%M:%S")),
    sprintf("End time:     %s", format(t_end,   "%Y-%m-%d %H:%M:%S")),
    sprintf("Elapsed time: %s seconds", elapsed),
    sprintf("Alpha code: %s", code),
    sprintf("Years used (%d): %s", n_years, paste(years, collapse = ", ")),
    "Input:",
    sprintf("  TimeSeries root: %s", file.path(root_dir, "TimeSeries")),
    sprintf("  Files loaded:    %d", length(files)),
    "Processing steps:",
    "  [1/4] Discover available years",
    "  [2/4] Load raster stack",
    "  [3/4] Compute Theil-Sen slope and Mann-Kendall statistics",
    "  [4/4] Write output rasters",
    "Outputs (GeoTIFF):",
    sprintf("  Trend: %s", tif_out$trend),
    sprintf("  p:     %s", tif_out$p),
    sprintf("  z:     %s", tif_out$z),
    "Outputs (AAIGrid):",
    sprintf("  Trend: %s", asc_out$trend),
    sprintf("  p:     %s", asc_out$p),
    sprintf("  z:     %s", asc_out$z),
    sprintf("PNG output: %s", ifelse(is.na(png_path), "Not generated", png_path)),
    "Diagnostics:",
    sprintf("  Total cells:        %d", n_cells),
    sprintf("  Valid slope cells:  %d", n_valid_slope),
    sprintf("  NA slope cells:     %d", n_na_slope)
  )

  cat(log_block, file = log_path, sep = "\n", append = TRUE)

}
