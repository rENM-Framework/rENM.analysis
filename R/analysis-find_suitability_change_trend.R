#' Find climatic suitability change trend
#'
#' Builds a multi-year prediction stack for a species (by alpha code)
#' and computes robust Theil-Sen trend layers representing suitability
#' change over time (i.e., acceleration or deceleration).
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item GeoTIFF prediction layers located at:
#'   \code{<project_dir>/runs/<alpha_code>/TimeSeries/<year>/model/}
#'   \code{<alpha_code>-<year>-Prediction.tif}
#'   \item Years included: 1980-2020 in 5-year steps
#' }
#'
#' \strong{Processing}
#' \itemize{
#'   \item Creates a multi-layer prediction stack across all available years
#'   \item Creates a difference stack of consecutive layers:
#'   \code{year[i + 1] - year[i]}
#'   \item Computes robust Theil-Sen per-pixel trend layers:
#'   \itemize{
#'     \item \code{ts_pred}: trend across prediction stack
#'     \item \code{ts_diff}: trend across difference stack
#'   }
#'   \item \code{ts_diff} represents the suitability change trend
#'   (rate-of-change in suitability over time)
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item The Theil-Sen slope is computed per pixel as the median of all
#'   pairwise slopes across time:
#'   \itemize{
#'     \item \eqn{median((y_j - y_i) / (t_j - t_i))} for \eqn{i < j}
#'   }
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Writes \code{ts_diff} (suitability change trend) to:
#'   \itemize{
#'     \item \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'     \code{<alpha_code>-Suitability-Change-Trend.tif}
#'     \item \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'     \code{<alpha_code>-Suitability-Change-Trend.asc}
#'   }
#'   \item Removes GDAL sidecar files (.aux.xml, .prj) if created
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item At least two prediction layers must be present
#'   \item Missing yearly files are tolerated with warnings
#' }
#'
#' @param alpha_code Character. Four-letter species code used to resolve
#' input and output paths.
#'
#' @return List with the following components:
#' \itemize{
#'   \item \code{stack}: SpatRaster of prediction layers
#'   \item \code{diff}: SpatRaster of consecutive differences
#'   \item \code{ts_pred}: SpatRaster of Theil-Sen trend on predictions
#'   \item \code{ts_diff}: SpatRaster of Theil-Sen trend on differences
#' }
#' Side effects:
#' \itemize{
#'   \item Writes GeoTIFF and ASCII raster outputs to disk
#'   \item Removes GDAL sidecar files if present
#' }
#'
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#'   find_suitability_change_trend("CASP")
#' }
#'
#' @export
find_suitability_change_trend <- function(alpha_code) {
  # ---- Dependencies and input checks ----
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required. Install it with install.packages('terra').", call. = FALSE)
  }
  stopifnot(is.character(alpha_code), length(alpha_code) == 1, nzchar(alpha_code))

  # ---- Project directory (CRAN-compliant) ----
  project_dir <- rENM_project_dir()

  # ---- Build list of expected files (1980-2020 by 5) ----
  years <- seq(1980, 2020, by = 5)
  base_dir <- project_dir
  files <- file.path(
    base_dir, "runs", alpha_code, "TimeSeries", years, "model",
    sprintf("%s-%d-Prediction.tif", alpha_code, years)
  )

  # Determine which files exist; warn if any are missing
  keep <- file.exists(files)
  if (!any(keep)) {
    stop("No prediction files found. Check paths/alpha_code and file naming.", call. = FALSE)
  }
  if (any(!keep)) {
    warning(
      sprintf("Missing files for years: %s", paste(years[!keep], collapse = ", ")),
      call. = FALSE
    )
  }

  years_found <- years[keep]

  # ---- Read available prediction rasters into a single SpatRaster stack ----
  stk <- terra::rast(files[keep])
  names(stk) <- sprintf("%s_%d", alpha_code, years_found)

  if (terra::nlyr(stk) < 2) {
    stop("Need at least two layers to compute differences and trends.", call. = FALSE)
  }

  # ---- Build difference stack ----
  diffs <- lapply(2:terra::nlyr(stk), function(i) stk[[i]] - stk[[i - 1]])
  ds <- do.call(c, diffs)
  names(ds) <- sprintf(
    "%s_%d_minus_%d",
    alpha_code, years_found[-1], years_found[-length(years_found)]
  )

  # ---- Theil-Sen slope helper ----
  theil_sen_slope <- function(y, t) {
    ok <- is.finite(y) & is.finite(t)
    y <- y[ok]; t <- t[ok]
    if (length(y) < 2) return(NA_real_)
    dy <- outer(y, y, "-")
    dt <- outer(t, t, "-")
    stats::median(dy[upper.tri(dy)] / dt[upper.tri(dt)], na.rm = TRUE)
  }

  # ---- Compute Theil-Sen trend layers ----
  t_pred <- years_found
  ts_pred <- terra::app(stk, fun = function(v) theil_sen_slope(v, t_pred))
  names(ts_pred) <- "ts_pred_slope"

  t_diff <- years_found[-1]
  ts_diff <- terra::app(ds, fun = function(v) theil_sen_slope(v, t_diff))
  names(ts_diff) <- "ts_diff_slope"

  # ---- Save outputs ----
  out_dir <- file.path(base_dir, "runs", alpha_code, "Trends", "suitability")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  tif_path <- file.path(out_dir, sprintf("%s-Suitability-Change-Trend.tif", alpha_code))
  asc_path <- file.path(out_dir, sprintf("%s-Suitability-Change-Trend.asc", alpha_code))

  terra::writeRaster(ts_diff, tif_path, overwrite = TRUE)

  asc_ok <- TRUE
  tryCatch({
    terra::writeRaster(ts_diff, asc_path, overwrite = TRUE)
  }, error = function(e) {
    asc_ok <<- FALSE
    tryCatch({
      terra::writeRaster(ts_diff, asc_path, overwrite = TRUE, gdal = "AAIGrid")
      asc_ok <<- TRUE
    }, error = function(e2) {
      stop(
        paste0(
          "Failed to write ASCII Grid (.asc). Your GDAL/terra build may lack the AAIGrid driver.\n",
          "GeoTIFF was written to: ", tif_path
        ),
        call. = FALSE
      )
    })
  })

  # ---- Remove GDAL sidecar files ----
  sidecars <- list.files(
    out_dir,
    pattern = paste0("^", alpha_code, ".*\\.(aux\\.xml|prj)$"),
    full.names = TRUE
  )
  if (length(sidecars) > 0) invisible(file.remove(sidecars))

  # ---- Final messages ----
  message("Saved suitability change trend raster (ts_diff):")
  message("  ", tif_path)
  if (asc_ok) message("  ", asc_path)

  # ---- Return ----
  list(
    stack   = stk,
    diff    = ds,
    ts_pred = ts_pred,
    ts_diff = ts_diff
  )
}
