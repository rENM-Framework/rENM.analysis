#' Plot climatic suitability change trend
#'
#' Create a publication-ready map of climatic suitability change trend
#' values using a diverging color scheme. Positive values indicate
#' acceleration, while negative values indicate deceleration.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' This function visualizes raster-based climatic suitability change
#' trends derived from retrospective ecological niche modeling
#' workflows. It integrates raster outputs with species range and
#' administrative boundary overlays to produce publication-quality
#' figures.
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Species identifier used to retrieve GAP range information
#'   via \code{get_species_info()}.
#'   \item Raster input representing suitability change trend values,
#'   either as a file path or an in-memory SpatRaster.
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item A ggplot object containing raster visualization with vector
#'   overlays.
#'   \item Diagnostic console output reporting number of cells plotted.
#'   \item Invisible list containing plot, symmetric limits, and cell
#'   count.
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Raster values are converted to a data frame for plotting.
#'   \item Symmetric color limits are derived from the (0.02, 0.98)
#'   quantiles of the trend distribution.
#'   \item A diverging red-white-green palette is applied using
#'   \code{scale_fill_gradient2()}.
#'   \item Spatial alignment is enforced by projecting vector data to
#'   match the raster coordinate reference system.
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item U.S. state boundaries located at:
#'   \code{<project_dir>/data/shapefiles/}
#'   \code{tl_2012_us_state/tl_2012_us_state.shp}
#'   \item Species range shapefile located at:
#'   \code{<project_dir>/data/shapefiles/<gap_range>/}
#'   \code{<gap_range>.shp}
#' }
#'
#' @param alpha_code Character. Species alpha code used to retrieve
#' GAP range information via \code{get_species_info()}.
#' @param raster_file Character, SpatRaster. Raster representing
#' climatic suitability change trend values, provided either as a file
#' path or an in-memory object.
#' @param zero_band_frac Numeric. Fraction (0, 0.90) defining the width
#' of the near-zero band used in diverging color scaling.
#'
#' @return
#' A List with the following elements:
#' \itemize{
#'   \item \code{plot}: ggplot object representing the rendered map.
#'   \item \code{lims}: Numeric vector of symmetric color scale limits.
#'   \item \code{n_cells}: Integer count of raster cells plotted.
#' }
#' Side effects:
#' \itemize{
#'   \item Console output reporting processing status and cell counts.
#' }
#'
#' @importFrom terra rast vect crop as.data.frame crs ext
#' @importFrom sf read_sf st_transform st_as_sf
#' @importFrom ggplot2 ggplot geom_raster geom_sf aes scale_fill_gradient2
#' @importFrom ggplot2 coord_sf labs theme_minimal theme element_rect
#' @importFrom ggplot2 element_line element_blank element_text unit
#' @importFrom paletteer paletteer_c
#' @importFrom scales rescale
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun{
#' result <- plot_suitability_change_trend(
#'   alpha_code = "CASSIN",
#'   raster_file = "path/to/trend.tif"
#' )
#' print(result$plot)
#' }
#'
#' @export
plot_suitability_change_trend <- function(alpha_code, raster_file, zero_band_frac = 0.80) {
  # ---- Dependencies ----
  if (!requireNamespace("terra", quietly = TRUE))     stop("Package 'terra' is required.")
  if (!requireNamespace("sf", quietly = TRUE))        stop("Package 'sf' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE))   stop("Package 'ggplot2' is required.")
  if (!requireNamespace("paletteer", quietly = TRUE)) stop("Package 'paletteer' is required.")
  if (!requireNamespace("scales", quietly = TRUE))    stop("Package 'scales' is required.")

  # ---- Species info ----
  code <- toupper(alpha_code %||% "")
  if (!nzchar(code)) stop("`alpha_code` must be a non-empty character string.")

  # get_species_info() is assumed to be available in the package or search path
  if (!exists("get_species_info", mode = "function")) {
    stop("Function `get_species_info()` must be available (e.g., exported by this package).")
  }

  si <- try(get_species_info(code), silent = TRUE)
  if (inherits(si, "try-error") || is.null(si)) {
    stop("get_species_info() failed for alpha_code: ", code)
  }
  if (!("GAP.RANGE" %in% names(si))) {
    stop("`GAP.RANGE` column not found for: ", code)
  }

  gap_range <- as.character(si[["GAP.RANGE"]])[1]
  if (!nzchar(gap_range) || is.na(gap_range)) {
    stop("Missing GAP.RANGE for alpha_code: ", code)
  }

  # ---- Paths ----
  cat("------------------------------------------------------------------------\nplot_suitability_change_trend(): starting (divergent colors)\n")

  # Use project directory helper
  if (!exists("rENM_project_dir", mode = "function")) {
    stop("Function `rENM_project_dir()` must be available (e.g., exported by this package).")
  }
  project_dir <- rENM_project_dir()

  states_path <- file.path(project_dir, "data", "shapefiles",
                           "tl_2012_us_state", "tl_2012_us_state.shp")
  if (!file.exists(states_path)) stop(sprintf("States shapefile not found: %s", states_path))

  range_dir  <- file.path(project_dir, "data", "shapefiles", gap_range)
  range_path <- file.path(range_dir, paste0(gap_range, ".shp"))
  if (!file.exists(range_path)) {
    stop("Species range shapefile not found at: ", range_path)
  }

  # ---- Raster input ----
  if (inherits(raster_file, "SpatRaster")) {
    r <- raster_file
  } else if (is.character(raster_file) && length(raster_file) == 1L) {
    if (!file.exists(raster_file)) stop(sprintf("Raster not found: %s", raster_file))
    r <- terra::rast(raster_file)
  } else {
    stop("`raster_file` must be a file path or a terra::SpatRaster.")
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
  if (!nrow(df)) stop("Raster has no finite values to plot.")
  names(df)[ncol(df)] <- "trend"
  n_cells <- nrow(df)

  # ---- Symmetric limits ----
  qs <- stats::quantile(df$trend, probs = c(0.02, 0.98), na.rm = TRUE, names = FALSE)
  max_abs <- max(abs(qs), na.rm = TRUE)
  lims <- c(-max_abs, max_abs)
  if (!all(is.finite(lims)) || diff(lims) == 0) {
    rng <- range(df$trend, na.rm = TRUE)
    lims <- max(abs(rng)) * c(-1, 1)
    max_abs <- max(abs(rng), na.rm = TRUE)
  }

  # ---- Palette ----
  pal_full <- as.character(
    paletteer::paletteer_c("ggthemes::Classic Area Red-Green", n = 401)
  )
  pick_col <- function(p) {
    pal_full[pmax(1, pmin(length(pal_full), round(p * (length(pal_full) - 1)) + 1))]
  }

  zero_band_frac <- max(0, min(0.90, zero_band_frac))
  z0 <- max_abs * zero_band_frac

  grad_values <- scales::rescale(
    c(-max_abs, -z0, 0, +z0, +max_abs),
    to   = c(0, 1),
    from = c(-max_abs, +max_abs)
  )
  grad_colors <- c(
    pick_col(0.00),
    pick_col(0.35),
    "white",
    pick_col(0.65),
    pick_col(1.00)
  )

  # ---- Plot ----
  ex <- terra::ext(r)
  dx <- ex$xmax - ex$xmin; dy <- ex$ymax - ex$ymin
  pad <- 0.02
  xlim <- c(ex$xmin - pad * dx, ex$xmax + pad * dx)
  ylim <- c(ex$ymin - pad * dy, ex$ymax + pad * dy)

  gp <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = df,
      ggplot2::aes(x = x, y = y, fill = trend)
    ) +
    ggplot2::geom_sf(
      data  = states,
      fill  = NA,
      color = "grey40",
      linewidth = 0.3
    ) +
    ggplot2::geom_sf(
      data  = range,
      fill  = NA,
      color = "black",
      linewidth = 0.7
    ) +
    ggplot2::scale_fill_gradient2(
      name     = "Change Trend",
      low      = "#D64241",
      mid      = "white",
      high     = "#58B440",
      midpoint = 0
    ) +
    ggplot2::coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    ggplot2::labs(
      title    = paste(code, "Climatic Suitability Change Trend (1980-2020)"),
      subtitle = "Suitability rate-of-change trend (acceleration/deceleration)",
      x = "Longitude",
      y = "Latitude"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.background  = ggplot2::element_rect(fill = "grey98", color = NA),
      panel.border      = ggplot2::element_rect(
        color = "grey30",
        fill  = NA,
        linewidth = 0.6
      ),
      panel.grid.major  = ggplot2::element_line(
        linewidth = 0.3,
        color     = "grey88"
      ),
      panel.grid.minor  = ggplot2::element_blank(),
      axis.title        = ggplot2::element_text(face = "plain"),
      axis.text         = ggplot2::element_text(color = "grey20"),
      legend.title      = ggplot2::element_text(size = 9),
      legend.text       = ggplot2::element_text(size = 8),
      legend.key.height = ggplot2::unit(1.0, "lines"),
      legend.key.width  = ggplot2::unit(0.4, "lines"),
      legend.position   = "right",
      plot.title        = ggplot2::element_text(
        face = "bold",
        size = 12,
        hjust = 0
      )
    )

  cat(sprintf("Cells plotted: %s\n",
              format(n_cells, big.mark = ",", scientific = FALSE)))
  cat("plot_suitability_change_trend(): done. Returning ggplot object.\n")
  invisible(list(plot = gp, lims = lims, n_cells = n_cells))
}

# Small infix helper
`%||%` <- function(a, b) if (is.null(a)) b else a
