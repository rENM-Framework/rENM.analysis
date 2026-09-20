#' Summarize suitability trends inside the range and in the buffer ring
#'
#' Compares trend behavior within the species USGS GAP range against the
#' buffer ring surrounding it: the zone between the raw GAP polygon and the
#' buffered polygon written by \code{find_range_extent()}.
#'
#' @details
#' \strong{Why the ring}
#' A static range polygon cannot show a range edge that is moving. Strong
#' positive trends concentrated just outside the historic boundary are the
#' signature of an advancing leading edge, and a range-based statistic hides
#' that entirely by never looking beyond the boundary. Comparing the ring
#' against the interior makes the pattern visible.
#'
#' A ring that is markedly more positive than the interior suggests
#' conditions improving where the species would expand into. A ring that is
#' markedly less positive suggests the opposite: a range with nowhere
#' favorable to go.
#'
#' \strong{Why this is reported range-wide rather than per state}
#' The ring extends well beyond the range, so it reaches into states holding
#' no range at all, for which a per-state ring figure has no counterpart in
#' the range-based table. More importantly, a statistic computed over few
#' raster cells is unstable between model realizations, whereas the
#' range-wide ring and interior figures reproduce closely. The comparison is
#' only meaningful at a scale where both sides are well estimated.
#'
#' \strong{Inputs}
#' \itemize{
#'   \item \code{data/shapefiles/<GAP.RANGE>/<GAP.RANGE>.shp}, the raw range.
#'   \item \code{runs/<alpha_code>/_occs/range_buffered.gpkg}, written by
#'     \code{find_range_extent()}. The ring needs this polygon itself, not
#'     the bounding box in \code{extent.txt}: a rectangle would pull in
#'     large areas near its corners that are nowhere near the range
#'     boundary, diluting the signal the statistic exists to detect.
#'   \item \code{runs/<alpha_code>/Trends/suitability/<alpha_code>-Suitability-Trend.tif}
#'   \item \code{runs/<alpha_code>/Trends/suitability/<alpha_code>-Hot-Spot-Mask.tif}
#' }
#'
#' \strong{How areas are measured}
#' Cells straddling a zone boundary contribute only the fraction of their
#' area lying inside it, matching \code{create_hot_spot_map()}. Percentages
#' are taken over the area that carries trend data, which can be smaller
#' than the zone where the model produced no prediction.
#'
#' @param alpha_code Character scalar. Four-letter species alpha code.
#'
#' @return Invisibly returns a data.frame with one row per zone
#'   (\code{"interior"} and \code{"ring"}) and columns \code{zone},
#'   \code{area_km2}, \code{data_area_km2}, \code{pos_pct}, \code{neg_pct},
#'   \code{hotspot_area_km2}, \code{hotspot_pct}.
#'
#' @seealso \code{\link{create_hot_spot_map}},
#'   \code{\link[rENM.data]{find_range_extent}}
#'
#' @importFrom sf st_read st_make_valid st_union st_difference st_transform
#' @importFrom sf st_geometry st_crs
#' @importFrom terra rast crs vect rasterize cellSize global project resample
#' @importFrom terra same.crs compareGeom
#' @importFrom utils read.csv
#' @importFrom readr write_csv
#'
#' @examples
#' \dontrun{
#'   find_boundary_trend_statistics("CASP")
#' }
#'
#' @export
find_boundary_trend_statistics <- function(alpha_code) {

  t_start <- Sys.time()
  step_i  <- 0L
  bump <- function(label) {
    step_i <<- step_i + 1L
    message(sprintf("[Step %02d] %s", step_i, label))
  }

  for (p in c("sf", "terra", "readr")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required package '", p, "' is not installed.", call. = FALSE)
    }
  }
  if (!is.character(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  code        <- toupper(alpha_code)
  project_dir <- rENM_project_dir()
  runs_root   <- file.path(project_dir, "runs", code)
  out_dir     <- file.path(runs_root, "Trends", "suitability")

  trend_path  <- file.path(out_dir, sprintf("%s-Suitability-Trend.tif", code))
  hs_path     <- file.path(out_dir, sprintf("%s-Hot-Spot-Mask.tif", code))
  buffer_path <- file.path(runs_root, "_occs", "range_buffered.gpkg")
  out_csv     <- file.path(out_dir,
                   sprintf("%s-Suitability-Trend-Boundary-Statistics.csv", code))
  log_file    <- file.path(runs_root, "_log.txt")

  # ---- Resolve the raw GAP shapefile --------------------------------------
  bump("Resolving GAP range shapefile")
  species_csv <- file.path(project_dir, "data", "_species.csv")
  if (!file.exists(species_csv)) {
    stop("Species table not found at: ", species_csv, call. = FALSE)
  }
  sp       <- utils::read.csv(species_csv, stringsAsFactors = FALSE, check.names = FALSE)
  norm     <- function(x) gsub("[^A-Z0-9]", "", toupper(x))
  col_norm <- norm(names(sp))
  a_idx    <- match("ALPHACODE", col_norm, nomatch = 0L)
  g_idx    <- match("GAPRANGE",  col_norm, nomatch = 0L)
  if (a_idx == 0L || g_idx == 0L) {
    stop("Species table must contain ALPHA.CODE and GAP.RANGE columns.", call. = FALSE)
  }
  row_idx <- which(toupper(trimws(sp[[names(sp)[a_idx]]])) == code)
  if (!length(row_idx)) {
    stop("No row found with alpha code '", code, "' in: ", species_csv, call. = FALSE)
  }
  gap_range <- trimws(sp[[names(sp)[g_idx]]][row_idx[1L]])
  gap_path  <- file.path(project_dir, "data", "shapefiles",
                         gap_range, paste0(gap_range, ".shp"))

  for (f in c(gap_path, buffer_path, trend_path, hs_path)) {
    if (!file.exists(f)) stop("Required input not found: ", f, call. = FALSE)
  }

  # ---- Build the interior and ring polygons -------------------------------
  bump("Differencing buffered and raw range polygons")
  raw_geom <- sf::st_union(sf::st_make_valid(sf::st_read(gap_path, quiet = TRUE)))
  buffered <- sf::st_geometry(sf::st_read(buffer_path, quiet = TRUE))

  raw_eq  <- sf::st_transform(raw_geom, sf::st_crs(buffered))
  ring_eq <- suppressWarnings(sf::st_difference(buffered, raw_eq))

  # ---- Load rasters -------------------------------------------------------
  bump("Reading trend and hot-spot rasters")
  trend <- terra::rast(trend_path)
  hs    <- terra::rast(hs_path)
  if (!terra::same.crs(trend, hs)) {
    hs <- terra::project(hs, trend, method = "near")
  }
  if (!terra::compareGeom(trend, hs, stopOnError = FALSE)) {
    hs <- terra::resample(hs, trend, method = "near")
  }
  cell_km2 <- terra::cellSize(trend, unit = "km")
  rast_crs <- terra::crs(trend, proj = TRUE)

  # ---- Summarize a zone ---------------------------------------------------
  # Coverage weighting, as in create_hot_spot_map(): a cell straddling the
  # zone boundary contributes only the fraction of its area inside it.
  # Percentages are taken over the area carrying trend data, which is
  # smaller than the zone wherever the model produced no prediction.
  .zone_stats <- function(geom, label) {
    v   <- terra::vect(sf::st_transform(geom, rast_crs))
    cov <- terra::rasterize(v, trend, cover = TRUE)

    wt        <- cell_km2 * cov
    area      <- terra::global(wt, "sum", na.rm = TRUE)[1, 1]
    data_area <- terra::global(wt * !is.na(trend), "sum", na.rm = TRUE)[1, 1]
    pos_area  <- terra::global(wt * (trend > 0), "sum", na.rm = TRUE)[1, 1]
    neg_area  <- terra::global(wt * (trend < 0), "sum", na.rm = TRUE)[1, 1]
    hot_area  <- terra::global(wt * (hs == 1),   "sum", na.rm = TRUE)[1, 1]

    zero_if_na <- function(x) if (is.null(x) || is.na(x)) 0 else x
    area      <- zero_if_na(area)
    data_area <- zero_if_na(data_area)
    pos_area  <- zero_if_na(pos_area)
    neg_area  <- zero_if_na(neg_area)
    hot_area  <- zero_if_na(hot_area)

    pct <- function(x) if (data_area > 0) 100 * x / data_area else NA_real_

    data.frame(
      zone             = label,
      area_km2         = area,
      data_area_km2    = data_area,
      pos_pct          = pct(pos_area),
      neg_pct          = pct(neg_area),
      hotspot_area_km2 = hot_area,
      hotspot_pct      = pct(hot_area),
      stringsAsFactors = FALSE
    )
  }

  bump("Summarizing range interior")
  interior <- .zone_stats(raw_eq,  "interior")
  bump("Summarizing buffer ring")
  ring     <- .zone_stats(ring_eq, "ring")

  out <- rbind(interior, ring)
  readr::write_csv(out, out_csv)
  message("Boundary statistics written: ", out_csv)

  # ---- Log ----------------------------------------------------------------
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  sep72   <- paste(rep("-", 72), collapse = "")
  f <- function(k, v) sprintf("%-18s %s", k, v)
  try({
    cat(paste0(paste(c(
      "",
      sep72,
      "Processing summary (find_boundary_trend_statistics)",
      f("Timestamp:",     format(Sys.time(), tz = "UTC", usetz = TRUE)),
      f("Alpha code:",    code),
      f("GAP range:",     gap_range),
      f("Buffered poly:", buffer_path),
      f("Interior pos %:", sprintf("%.3f", interior$pos_pct)),
      f("Ring pos %:",     sprintf("%.3f", ring$pos_pct)),
      f("Ring - interior:", sprintf("%+.3f points", ring$pos_pct - interior$pos_pct)),
      f("Output file:",   out_csv),
      f("Total elapsed:", sprintf("%.1f s", elapsed)),
      ""
    ), collapse = "\n")), file = log_file, append = TRUE)
  }, silent = TRUE)

  invisible(out)
}
