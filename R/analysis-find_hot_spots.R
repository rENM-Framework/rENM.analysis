#' Find areas where declining suitability is accelerating (hot spots)
#'
#' Create a binary hot-spot raster identifying cells where baseline
#' suitability trend is negative and suitability change trend is
#' positive. Hot-spot cells are coded as 1 and all others as 0.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' Creates a binary hot-spot raster identifying cells where:
#' A (<alpha_code>-Suitability-Trend.tif) < 0 AND
#' B (<alpha_code>-Suitability-Difference-Trend.tif) > 0.
#' Hot-spot cells are coded 1; all other cells are coded 0.
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Baseline suitability trend raster (A).
#'   \item Suitability change trend raster (B).
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Binary hot-spot raster written to disk.
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Input rasters are read from:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#'   \code{<alpha_code>-Suitability-Difference-Trend.tif}
#'   \item If grids differ in extent, resolution, or CRS, raster B is
#'   resampled or projected to match raster A.
#'   \item Binary mask is computed where (A < 0) AND (B > 0).
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Input rasters located at:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#'   \code{<alpha_code>-Suitability-Difference-Trend.tif}
#'   \item Output written to:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-CCEI-Baseline.tif}
#' }
#'
#' @param alpha_code Character. Species alpha code (e.g., "CASP").
#'
#' @return
#' A SpatRaster object representing the binary (0, 1) hot-spot mask.
#' Side effects:
#' \itemize{
#'   \item Output raster is written to disk as a compressed GeoTIFF.
#' }
#'
#' @importFrom terra rast same.crs project compareGeom resample ifel writeRaster
#'
#' @examples
#' \dontrun{
#' # Example:
#' # hs <- find_hot_spots("CASP")
#' # terra::plot(hs)
#' }
#'
#' @export
find_hot_spots <- function(alpha_code) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required but not installed.")
  }

  if (!exists("rENM_project_dir", mode = "function")) {
    stop("Function `rENM_project_dir()` must be available (e.g., exported by this package).")
  }

  project_dir <- rENM_project_dir()

  base_dir <- file.path(project_dir, "runs", alpha_code, "Trends", "suitability")
  a_path   <- file.path(base_dir, sprintf("%s-Suitability-Trend.tif", alpha_code))
  b_path   <- file.path(base_dir, sprintf("%s-Suitability-Change-Trend.tif", alpha_code))
  out_path <- file.path(base_dir, sprintf("%s-Hot-Spot-Mask.tif", alpha_code))

  dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(a_path)) stop("Missing input A: ", a_path)
  if (!file.exists(b_path)) stop("Missing input B: ", b_path)

  A <- terra::rast(a_path)
  B <- terra::rast(b_path)

  # Align CRS if different
  if (!terra::same.crs(A, B)) {
    B <- terra::project(B, A, method = "near")
  }

  # Ensure same geometry (extent/resolution/origin)
  if (!terra::compareGeom(A, B, stopOnError = FALSE)) {
    B <- terra::resample(B, A, method = "near")
  }

  # Compute binary mask (A < 0 and B > 0)
  cond <- (A < 0) & (B > 0)
  C <- terra::ifel(cond, 1, 0)
  C[is.na(C)] <- 0

  # Save output as compressed 8-bit integer TIFF
  C <- terra::writeRaster(
    C,
    filename = out_path,
    overwrite = TRUE,
    datatype = "INT1U",
    gdal = c("COMPRESS=LZW", "PREDICTOR=2")
  )

  message("Hot-spot raster written: ", out_path)
  return(C)
}
