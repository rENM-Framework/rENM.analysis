#' Plot suitability trend map with centroids
#'
#' Create a publication-ready map of climatic suitability trend values
#' with overlaid geographic context and centroid displacement. The map
#' highlights increasing and decreasing suitability patterns spatially.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Climatic suitability trend raster:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#'   or
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.asc}
#'   \item U.S. state boundaries shapefile:
#'   \code{<project_dir>/data/shapefiles/tl_2012_us_state/}
#'   \code{tl_2012_us_state.shp}
#'   \item Species range shapefile derived from
#'   \code{get_species_info(alpha_code)} using the \code{GAP.RANGE}
#'   field in
#'   \code{<project_dir>/data/}
#'   \code{_species.csv}
#'   and loaded from
#'   \code{<project_dir>/data/shapefiles/<gap_range>/}
#'   \code{<gap_range>.shp}
#'   \item Bioclimatic centroid velocity CSV:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/centroids/}
#'   \code{<alpha_code>-Bioclimatic-Velocity.csv}
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item A \code{ggplot2} map object with:
#'   \itemize{
#'     \item Raster layer showing climatic suitability trends
#'     \item U.S. state boundaries (thin gray outlines, no fill)
#'     \item Species range outline (thin black outline)
#'     \item Bioclimatic centroid shift (yellow): a solid dot at the
#'     start position and a closed arrow head pointing to the end
#'     position
#'   }
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Raster values are visualized using a red-green diverging
#'   palette centered at zero
#'   \item Color limits are chosen symmetrically using the 2-98 percent
#'   quantile range to reduce the influence of outliers
#'   \item A widened neutral band around zero (white) sharpens the
#'   distinction between positive and negative values
#'   \item Default neutral band width is +/- 80 percent
#'   (\code{zero_band_frac = 0.80})
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Centroid CSV must contain columns:
#'   \code{start_lon}, \code{start_lat}, \code{end_lon}, \code{end_lat}
#'   (WGS84 longitude and latitude in degrees)
#'   \item A common typo \code{end_on} is tolerated and treated as
#'   \code{end_lon}
#'   \item If multiple rows are present in the centroid CSV, all
#'   segments are drawn
#'   \item All vector layers (including centroid points) are
#'   reprojected to the raster CRS (if defined) and cropped or fit to
#'   the raster extent for performance and correct overlay
#' }
#'
#' @param alpha_code Character. Four-letter species code (e.g., "CASP").
#' @param raster_file Character, SpatRaster. NULL (default) to
#' auto-discover, a file path to a raster file, or a
#' \code{terra::SpatRaster}. If NULL, the function searches:
#' \itemize{
#'   \item \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#'   \item \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.asc}
#' }
#' and uses the first one found.
#' @param zero_band_frac Numeric. Fraction of the symmetric limit to
#' allocate as a neutral band around zero in the diverging palette
#' (default 0.80 = +/- 80 percent; clamped to \code{(0, 0.90)}).
#'
#' @return Invisibly returns a List with:
#' \itemize{
#'   \item plot: ggplot2 object representing the rendered map
#'   \item lims: Numeric vector of color scale limits used
#'   \item n_cells: Integer number of raster cells plotted
#'   \item centroids: Data frame of centroid segments used in the plot
#'   (in raster CRS)
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Reads raster and vector data from the project directory
#'   \item Performs CRS transformations and spatial cropping
#'   \item Generates a publication-ready visualization
#' }
#'
#' @importFrom paletteer paletteer_c
#' @importFrom terra rast vect crop as.data.frame crs ext
#' @importFrom sf read_sf st_transform st_as_sf st_coordinates
#' @importFrom ggplot2 ggplot geom_raster geom_sf geom_segment geom_point
#' @importFrom ggplot2 scale_fill_gradientn coord_sf labs theme_minimal aes
#' @importFrom scales rescale squish
#' @importFrom stats quantile
#' @importFrom utils read.csv
#'
#' @examples
#' \dontrun{
#'   res <- plot_trend_with_centroids("CASP")
#'   print(res$plot)
#' }
#'
#' @export
plot_trend_with_centroids <- function(alpha_code,
                                      raster_file = NULL,
                                      zero_band_frac = 0.80) {

  project_dir <- rENM_project_dir()

  # R CMD check appeasement for NSE columns used in ggplot2::aes()
  start_lon <- start_lat <- end_lon <- end_lat <- NULL
  x <- y <- trend <- xend <- yend <- NULL

  # ---- Dependencies ----
  if (!requireNamespace("terra", quietly = TRUE))   stop("Package 'terra' is required.",   call. = FALSE)
  if (!requireNamespace("sf", quietly = TRUE))      stop("Package 'sf' is required.",      call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
  if (!requireNamespace("scales", quietly = TRUE))  stop("Package 'scales' is required.",  call. = FALSE)
  if (!requireNamespace("grid", quietly = TRUE))    stop("Package 'grid' is required.",    call. = FALSE)

  if (!requireNamespace("paletteer", quietly = TRUE)) {
    stop("Package 'paletteer' is required.", call. = FALSE)
  }

  # ---- Species info ----
  code <- toupper(alpha_code %||% "")
  if (!nzchar(code)) {
    stop("`alpha_code` must be a non-empty character string.", call. = FALSE)
  }

  if (!exists("get_species_info", mode = "function")) {
    stop("get_species_info() not found; is the package loaded?", call. = FALSE)
  }

  si <- try(get_species_info(code), silent = TRUE)
  if (inherits(si, "try-error") || is.null(si)) {
    stop("get_species_info() failed for alpha_code: ", code, call. = FALSE)
  }
  if (!("GAP.RANGE" %in% names(si))) {
    stop("`GAP.RANGE` column not found for: ", code, call. = FALSE)
  }

  gap_range <- as.character(si[["GAP.RANGE"]])[1]
  if (!nzchar(gap_range) || is.na(gap_range)) {
    stop("Missing GAP.RANGE for alpha_code: ", code, call. = FALSE)
  }

  # ---- Paths ----
  cat("------------------------------------------------------------------------\n",
      "plot_trend_with_centroids(): starting (divergent colors + centroid arrow)\n",
      sep = "")

  states_path <- file.path(project_dir, "data", "shapefiles", "tl_2012_us_state", "tl_2012_us_state.shp")
  if (!file.exists(states_path)) {
    stop(sprintf("States shapefile not found: %s", states_path), call. = FALSE)
  }

  range_dir  <- file.path(project_dir, "data", "shapefiles", gap_range)
  range_path <- file.path(range_dir, paste0(gap_range, ".shp"))
  if (!file.exists(range_path)) {
    stop("Species range shapefile not found at: ", range_path, call. = FALSE)
  }

  cent_path <- file.path(
    project_dir, "runs",
    code, "Trends", "centroids",
    paste0(code, "-Bioclimatic-Velocity.csv")
  )
  if (!file.exists(cent_path)) {
    stop("Centroid velocity CSV not found at: ", cent_path, call. = FALSE)
  }

  # ---- Resolve raster_file default ----
  if (is.null(raster_file)) {
    base_dir <- file.path(project_dir, "runs", code, "Trends", "suitability")
    cand_tif <- file.path(base_dir, paste0(code, "-Suitability-Trend.tif"))
    cand_asc <- file.path(base_dir, paste0(code, "-Suitability-Trend.asc"))
    if (file.exists(cand_tif)) {
      raster_file <- cand_tif
    } else if (file.exists(cand_asc)) {
      raster_file <- cand_asc
    } else {
      stop("No raster found. Looked for:\n  ", cand_tif, "\n  ", cand_asc,
           "\nProvide `raster_file` explicitly or generate the trend raster first.",
           call. = FALSE)
    }
  }

  # ---- Raster input ----
  if (inherits(raster_file, "SpatRaster")) {
    r <- raster_file
  } else if (is.character(raster_file) && length(raster_file) == 1L) {
    if (!file.exists(raster_file)) {
      stop(sprintf("Raster not found: %s", raster_file), call. = FALSE)
    }
    r <- terra::rast(raster_file)
  } else {
    stop("`raster_file` must be NULL (to auto-discover), a file path, or a terra::SpatRaster.",
         call. = FALSE)
  }

  # ---- Read vectors ----
  states <- sf::read_sf(states_path, quiet = TRUE)
  range  <- sf::read_sf(range_path,  quiet = TRUE)

  # ---- Align CRS ----
  r_crs <- try(terra::crs(r, proj = TRUE), silent = TRUE)
  if (!inherits(r_crs, "try-error") && !is.na(r_crs) && nchar(r_crs) > 0) {
    states <- sf::st_transform(states, r_crs)
    range  <- sf::st_transform(range,  r_crs)
  }

  # ---- Crop vectors ----
  states_v <- try(terra::vect(states), silent = TRUE)
  range_v  <- try(terra::vect(range),  silent = TRUE)
  if (!inherits(states_v, "try-error")) states <- sf::st_as_sf(terra::crop(states_v, r))
  if (!inherits(range_v,  "try-error")) range  <- sf::st_as_sf(terra::crop(range_v,  r))

  # ---- Raster -> data.frame ----
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  if (!nrow(df)) stop("Raster has no finite values to plot.", call. = FALSE)
  names(df)[ncol(df)] <- "trend"
  n_cells <- nrow(df)

  # ---- Symmetric limits ----
  qs <- stats::quantile(df$trend, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
  max_abs <- max(abs(qs), na.rm = TRUE)
  lims <- c(-max_abs, max_abs)
  if (!all(is.finite(lims)) || diff(lims) == 0) {
    rng <- range(df$trend, na.rm = TRUE)
    max_abs <- max(abs(rng), na.rm = TRUE)
    lims <- max_abs * c(-1, 1)
  }

  # ---- Palette ----
  pal_full <- as.character(
    paletteer::paletteer_c("ggthemes::Classic Area Red-Green", n = 401)
  )
  pick_col <- function(p) pal_full[
    pmax(1, pmin(length(pal_full), round(p * (length(pal_full) - 1)) + 1))
  ]

  zero_band_frac <- max(0, min(0.90, zero_band_frac))
  z0 <- max_abs * zero_band_frac

  grad_values <- scales::rescale(
    c(-max_abs, -z0, 0, +z0, +max_abs),
    to = c(0, 1),
    from = c(-max_abs, +max_abs)
  )
  grad_colors <- c(
    pick_col(0.00),
    pick_col(0.35),
    "white",
    pick_col(0.65),
    pick_col(1.00)
  )

  # ---- Read centroid CSV ----
  cent_raw <- try(utils::read.csv(cent_path, stringsAsFactors = FALSE), silent = TRUE)
  if (inherits(cent_raw, "try-error")) {
    stop("Failed to read centroid CSV at: ", cent_path, call. = FALSE)
  }
  if (!("end_lon" %in% names(cent_raw)) && ("end_on" %in% names(cent_raw))) {
    cent_raw$end_lon <- cent_raw$end_on
  }

  needed <- c("start_lon", "start_lat", "end_lon", "end_lat")
  if (!all(needed %in% names(cent_raw))) {
    stop("Centroid CSV must contain columns: ", paste(needed, collapse = ", "), call. = FALSE)
  }

  cent_raw <- subset(
    cent_raw,
    is.finite(start_lon) & is.finite(start_lat) &
      is.finite(end_lon) & is.finite(end_lat)
  )

  cent_starts_wgs <- sf::st_as_sf(
    cent_raw, coords = c("start_lon", "start_lat"),
    crs = 4326, remove = FALSE
  )
  cent_ends_wgs <- sf::st_as_sf(
    cent_raw, coords = c("end_lon", "end_lat"),
    crs = 4326, remove = FALSE
  )

  if (!inherits(r_crs, "try-error") && !is.na(r_crs) && nchar(r_crs) > 0) {
    cent_starts <- sf::st_transform(cent_starts_wgs, r_crs)
    cent_ends   <- sf::st_transform(cent_ends_wgs,   r_crs)
  } else {
    cent_starts <- cent_starts_wgs
    cent_ends   <- cent_ends_wgs
  }

  starts_xy <- sf::st_coordinates(cent_starts)
  ends_xy   <- sf::st_coordinates(cent_ends)
  cent_segments <- data.frame(
    x    = starts_xy[, 1],
    y    = starts_xy[, 2],
    xend = ends_xy[, 1],
    yend = ends_xy[, 2]
  )

  # ---- Plot ----
  ex <- terra::ext(r)
  dx <- ex$xmax - ex$xmin
  dy <- ex$ymax - ex$ymin
  pad <- 0.02
  xlim <- c(ex$xmin - pad * dx, ex$xmax + pad * dx)
  ylim <- c(ex$ymin - pad * dy, ex$ymax + pad * dy)

  gp <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes(x = x, y = y, fill = trend)) +
    ggplot2::geom_sf(data = states, fill = NA, color = "grey40", linewidth = 0.3) +
    ggplot2::geom_sf(data = range, fill = NA, color = "black", linewidth = 0.7) +
    ggplot2::scale_fill_gradientn(
      colors = grad_colors,
      values = grad_values,
      limits = lims,
      oob    = scales::squish,
      name   = "Trend"
    ) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(
      title    = paste(code, "Climatic Suitability Trend"),
      subtitle = "Suitability centroid shift (Circle = 1980, Diamond = 2020)",
      x = "Longitude", y = "Latitude"
    ) +
    ggplot2::theme_minimal(base_size = 13)

  if (nrow(cent_segments)) {
    gp <- gp +
      ggplot2::geom_segment(data = cent_segments,
                            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
                            color = "black", linewidth = 1.0) +
      ggplot2::geom_point(data = cent_segments,
                          ggplot2::aes(x = x, y = y),
                          shape = 21, fill = "orange", color = "black", size = 6.0) +
      ggplot2::geom_point(data = cent_segments,
                          ggplot2::aes(x = xend, y = yend),
                          shape = 23, fill = "orange", color = "black", size = 6.0)
  }

  invisible(list(plot = gp, lims = lims, n_cells = n_cells, centroids = cent_segments))
}

`%||%` <- function(a, b) if (is.null(a)) b else a
