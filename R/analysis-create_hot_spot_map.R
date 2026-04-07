#' Create and save a "hot spot" map for a species
#'
#' Build a publication-ready hot spot map for a species. In this experimental
#' metric we define hot spots as areas where declining trends are accelerating.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context}
#' Builds a publication-ready hot-spot map for a species. Hot spots are
#' where the baseline suitability trend is negative and the suitability
#' change trend (acceleration) is positive. Automatically generates
#' missing inputs and logs an eBird-standard processing summary.
#'
#' In addition, for U.S. states intersecting the species GAP range, this
#' function computes:
#' \itemize{
#'   \item hot-spot area (km^2) within each intersecting state, and
#'   \item percent of the state covered by hot spots,
#' }
#' and writes a per-state table to:
#' \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#' \code{<alpha_code>-Suitability-Trend-State-Analysis-Hotspots-Stats.csv}
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Species alpha code used to locate input rasters and shapefiles.
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Hot-spot map saved as PNG.
#'   \item Per-state statistics written as CSV.
#'   \item ggplot object returned invisibly.
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Hot-spot mask is created by \code{find_hot_spots()} if needed,
#'   using (A < 0) and (B > 0) with A = baseline trend and
#'   B = suitability change trend.
#'   \item Area and percentage method:
#'   \itemize{
#'     \item Hot-spot area is computed by summing cell areas (km^2) for
#'     cells with value 1 in the hot-spot mask within each state polygon.
#'     \item To avoid inflating percentages when rasters do not span an
#'     entire state, the state denominator is the area of the state
#'     clipped to the raster footprint (state intersect raster extent),
#'     in km^2.
#'     \item Percent coverage is 100 * hotspot_area_km2 / state_area_km2.
#'   }
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Non-CONUS areas (AK, HI, PR, GU, VI, AS, MP, UM) are excluded
#'   before intersecting with the GAP range.
#'   \item Input rasters and shapefiles must exist within:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#' }
#'
#' @param alpha_code Character. Four-letter species code (e.g., "CASP",
#' "LOSH").
#'
#' @return
#' A List containing:
#' \itemize{
#'   \item \code{plot}: ggplot object representing the rendered map.
#'   \item \code{png}: Character file path to the saved PNG map.
#'   \item \code{stats}: Character file path to the saved CSV file
#'   containing per-state hot-spot statistics.
#' }
#' Side effects:
#' \itemize{
#'   \item PNG map is written to disk.
#'   \item CSV file of per-state statistics is written to disk.
#'   \item Log entries are appended to:
#'   \code{<project_dir>/runs/<alpha_code>/_log.txt}
#' }
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_tile geom_sf geom_sf_text aes scale_fill_gradient2
#' @importFrom ggplot2 coord_sf labs theme_minimal theme element_rect element_text ggsave
#' @importFrom sf st_read st_make_valid st_crs st_transform st_union st_intersection
#' @importFrom sf st_intersects st_geometry st_point_on_surface st_bbox st_as_sfc st_sf
#' @importFrom terra rast same.crs project compareGeom resample vect mask crop ext crs
#' @importFrom terra cellSize expanse global
#' @importFrom dplyr filter bind_rows arrange desc
#' @importFrom readr write_csv
#'
#' @examples
#' \dontrun{
#' create_hot_spot_map("CASP")
#' }
#'
#' @export
create_hot_spot_map <- function(alpha_code) {
  t_start <- Sys.time()
  created_change_trend <- FALSE
  mask_recreated <- FALSE

  for (p in c("sf", "dplyr", "ggplot2", "terra", "readr")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package not installed: ", p, call. = FALSE)
    }
  }

  if (!exists("rENM_project_dir", mode = "function")) {
    stop("Function `rENM_project_dir()` must be available (e.g., exported by this package).")
  }
  project_dir <- rENM_project_dir()

  exists_any <- function(paths) any(file.exists(paths))
  code <- toupper(alpha_code)

  # Paths
  base_dir <- file.path(project_dir, "runs", code, "Trends", "suitability")
  states_path <- file.path(
    project_dir, "data", "shapefiles",
    "tl_2012_us_state", "tl_2012_us_state.shp"
  )
  gap_path <- file.path(
    project_dir, "data", "shapefiles",
    sprintf("b%sx_CONUS_Range_2001v1", code),
    sprintf("b%sx_CONUS_Range_2001v1.shp", code)
  )
  trend_tif   <- file.path(base_dir, sprintf("%s-Suitability-Trend.tif", code))
  trend_asc   <- file.path(base_dir, sprintf("%s-Suitability-Trend.asc", code))
  change_tif  <- file.path(base_dir, sprintf("%s-Suitability-Change-Trend.tif", code))
  change_asc  <- file.path(base_dir, sprintf("%s-Suitability-Change-Trend.asc", code))

  out_dir  <- base_dir
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  out_png  <- file.path(
    out_dir,
    sprintf("%s-Suitability-Trend-State-Analysis-Hotspots.png", code)
  )
  out_csv  <- file.path(
    out_dir,
    sprintf("%s-Suitability-Trend-State-Analysis-Hotspots-Stats.csv", code)
  )

  # ---- Check baseline trend ----
  if (!exists_any(c(trend_tif, trend_asc))) {
    stop(
      "Baseline trend raster not found:\n  ",
      trend_tif, "\n  ",
      trend_asc,
      call. = FALSE
    )
  }

  # ---- Ensure suitability-change-trend exists (use package function) ----
  if (!exists_any(c(change_tif, change_asc))) {
    message(
      "Suitability-Change-Trend raster not found - generating via ",
      "find_suitability_change_trend()..."
    )
    if (!exists("find_suitability_change_trend", mode = "function")) {
      stop(
        "Function find_suitability_change_trend() must be available ",
        "(e.g., exported by this package).",
        call. = FALSE
      )
    }
    find_suitability_change_trend(code)
    if (!exists_any(c(change_tif, change_asc))) {
      stop(
        "After running find_suitability_change_trend(), raster still missing:\n  ",
        change_tif, "\n  ",
        change_asc,
        call. = FALSE
      )
    }
    message("[OK] Created suitability change trend raster.")
    created_change_trend <- TRUE
  }

  # ---- Ensure hot-spot mask (use package function) ----
  if (!exists("find_hot_spots", mode = "function")) {
    stop(
      "Function find_hot_spots() must be available ",
      "(e.g., exported by this package).",
      call. = FALSE
    )
  }
  hs <- find_hot_spots(code)
  mask_recreated <- TRUE

  # ---- Load rasters ----
  trend_path <- if (file.exists(trend_tif)) trend_tif else trend_asc
  trend_raster <- terra::rast(trend_path)
  if (!terra::same.crs(trend_raster, hs)) {
    hs <- terra::project(hs, trend_raster, method = "near")
  }
  if (!terra::compareGeom(trend_raster, hs, stopOnError = FALSE)) {
    hs <- terra::resample(hs, trend_raster, method = "near")
  }

  # ---- Read vectors ----
  if (!file.exists(states_path)) stop("States shapefile not found: ", states_path)
  if (!file.exists(gap_path))    stop("GAP range shapefile not found: ", gap_path)
  states <- sf::st_read(states_path, quiet = TRUE) |>
    sf::st_make_valid()
  gap    <- sf::st_read(gap_path,    quiet = TRUE) |>
    sf::st_make_valid()
  if (sf::st_crs(states) != sf::st_crs(gap)) {
    gap <- sf::st_transform(gap, sf::st_crs(states))
  }

  # ---- Crop, plot, and save ------------------------------------
  non_conus <- c("AK", "HI", "PR", "GU", "VI", "AS", "MP", "UM")
  states_conus <- dplyr::filter(states, !.data$STUSPS %in% non_conus)
  gap_conus <- suppressWarnings(
    sf::st_intersection(sf::st_union(gap), sf::st_union(states_conus))
  )
  idx <- sf::st_intersects(states_conus, gap_conus, sparse = TRUE)
  states_in_gap <- states_conus[lengths(idx) > 0, , drop = FALSE]
  raster_wkt <- terra::crs(trend_raster, proj = TRUE)
  states_conus_tr  <- sf::st_transform(states_conus,  raster_wkt)
  states_in_gap_tr <- sf::st_transform(states_in_gap, raster_wkt)
  gap_conus_tr     <- sf::st_transform(gap_conus,     raster_wkt)

  states_conus_v <- terra::vect(states_conus_tr)
  trend_crop <- terra::mask(
    terra::crop(trend_raster, states_conus_v),
    states_conus_v
  )
  hs_crop    <- terra::mask(
    terra::crop(hs,           states_conus_v),
    states_conus_v
  )

  trend_df <- as.data.frame(trend_crop, xy = TRUE, na.rm = TRUE)
  names(trend_df)[3] <- "trend"
  hs_df <- as.data.frame(hs_crop, xy = TRUE, na.rm = TRUE)
  names(hs_df)[3] <- "hotspot"
  hs_df <- hs_df[hs_df$hotspot == 1, , drop = FALSE]

  lab_vals <- as.character(states_in_gap_tr$STUSPS)
  lab_pts_geom <- suppressWarnings(
    sf::st_point_on_surface(sf::st_geometry(states_in_gap_tr))
  )
  label_points <- sf::st_sf(
    label    = lab_vals,
    geometry = lab_pts_geom,
    crs      = sf::st_crs(states_in_gap_tr)
  )
  conus_bbox_tr <- sf::st_bbox(
    sf::st_union(
      sf::st_geometry(states_in_gap_tr),
      sf::st_geometry(gap_conus_tr)
    )
  )

  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = trend_df,
      ggplot2::aes(x = .data$x, y = .data$y, fill = .data$trend)
    ) +
    ggplot2::scale_fill_gradient2(
      name     = "Trend",
      low      = "darkgray",
      mid      = "white",
      high     = "darkgray",
      midpoint = 0
    ) +
    ggplot2::geom_tile(
      data        = hs_df,
      ggplot2::aes(x = .data$x, y = .data$y),
      fill        = "orange",
      alpha       = 0.85,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_sf(
      data      = states_conus_tr,
      fill      = NA,
      color     = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_sf(
      data      = states_in_gap_tr,
      fill      = NA,
      color     = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_sf(
      data      = gap_conus_tr,
      fill      = NA,
      color     = "black",
      linewidth = 0.7,
      linetype  = "dashed"
    ) +
    ggplot2::geom_sf_text(
      data = label_points,
      ggplot2::aes(label = .data$label),
      size     = 6,
      color    = "black",
      fontface = "bold"
    ) +
    ggplot2::coord_sf(
      xlim   = c(conus_bbox_tr["xmin"], conus_bbox_tr["xmax"]),
      ylim   = c(conus_bbox_tr["ymin"], conus_bbox_tr["ymax"]),
      expand = FALSE,
      crs    = raster_wkt
    ) +
    ggplot2::labs(
      title    = sprintf("%s Hot Spots within GAP Range States", code),
      subtitle = "Orange = Areas of accelerating declines in climatic suitability",
      x        = "Longitude",
      y        = "Latitude"
    ) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.border     = ggplot2::element_rect(fill = NA, color = "darkgray", linewidth = 0.8),
      plot.title       = ggplot2::element_text(size = 28, face = "bold", hjust = 0.0),
      plot.subtitle    = ggplot2::element_text(size = 20),
      axis.title       = ggplot2::element_text(size = 20),
      axis.text        = ggplot2::element_text(size = 14)
    )

  ggplot2::ggsave(
    out_png,
    plot   = p,
    width  = 1000,
    height = 750,
    units  = "px",
    dpi    = 72,
    bg     = "white"
  )

  {
    cell_km2 <- terra::cellSize(hs, unit = "km")
    rst_ext_poly_sf <- sf::st_as_sfc(
      sf::st_bbox(terra::ext(hs), crs = raster_wkt)
    )

    .pick_col <- function(x, candidates) {
      nm <- names(x)
      k <- candidates[candidates %in% nm][1]
      ifelse(length(k) == 0, NA_character_, k)
    }
    nm_col <- .pick_col(states_in_gap_tr, c("NAME","name","STATE_NAME","State","state","STATE"))
    ab_col <- .pick_col(states_in_gap_tr, c("STUSPS","stusps","postal","abbr","STATE_ABBR","STATEFP"))
    get_vals <- function(df, i, col, fallback) {
      if (!is.na(col)) as.character(df[[col]][i]) else fallback
    }

    out_rows <- vector("list", nrow(states_in_gap_tr))
    for (i in seq_len(nrow(states_in_gap_tr))) {
      st_sf <- states_in_gap_tr[i, , drop = FALSE]
      st_name <- get_vals(states_in_gap_tr, i, nm_col, as.character(i))
      st_abbr <- get_vals(states_in_gap_tr, i, ab_col, NA_character_)

      st_clip <- suppressWarnings(sf::st_intersection(st_sf, rst_ext_poly_sf))
      if (nrow(st_clip) == 0) {
        out_rows[[i]] <- data.frame(
          state = st_name, abbr = st_abbr,
          state_area_km2 = 0, hotspot_area_km2 = 0,
          hotspot_pct_of_state = NA_real_, stringsAsFactors = FALSE
        )
        next
      }
      st_clip_sv <- terra::vect(st_clip)
      state_area_km2 <- terra::expanse(st_clip_sv, unit = "km")
      if (length(state_area_km2) == 0 || is.na(state_area_km2)) state_area_km2 <- 0

      st_sv <- terra::vect(st_sf)
      hs_masked   <- terra::mask(hs, st_sv)
      area_masked <- terra::mask(cell_km2, st_sv)
      h_area <- terra::global(area_masked * (hs_masked == 1), "sum", na.rm = TRUE)[1,1]
      if (is.na(h_area)) h_area <- 0

      pct <- if (state_area_km2 > 0) 100 * (h_area / state_area_km2) else NA_real_

      out_rows[[i]] <- data.frame(
        state = st_name, abbr = st_abbr,
        state_area_km2 = state_area_km2,
        hotspot_area_km2 = h_area,
        hotspot_pct_of_state = pct,
        stringsAsFactors = FALSE
      )
    }

    stats_df <- dplyr::bind_rows(out_rows) |>
      dplyr::arrange(dplyr::desc(.data$hotspot_area_km2))

    readr::write_csv(stats_df, out_csv)
    message("Per-state hot-spot stats written: ", out_csv)
  }

  # ---- Logging ----
  t_end <- Sys.time()
  elapsed <- difftime(t_end, t_start, units = "secs")
  s <- as.numeric(elapsed)
  hh <- sprintf("%02d", s %/% 3600); s <- s %% 3600
  mm <- sprintf("%02d", s %/% 60);   ss <- sprintf("%02d", round(s %% 60))
  fmt_elapsed <- paste0(hh, ":", mm, ":", ss)

  log_file <- file.path(project_dir, "runs", code, "_log.txt")

  sep_line <- paste(rep("-", 72), collapse = "")
  header   <- "Processing summary (create_hot_spot_map)"
  ts_str   <- format(t_end, "%Y-%m-%d %H:%M:%S %Z")
  f <- function(k, v) sprintf("%-18s : %s", k, v)

  outputs <- c(out_png, out_csv)
  log_block <- c(
    "",
    sep_line,
    header,
    f("Timestamp", ts_str),
    f("Alpha code", code),
    f("Change trend made", if (created_change_trend) "TRUE" else "FALSE"),
    f("Hot-spot mask made", if (mask_recreated) "TRUE" else "FALSE"),
    f("States shapefile", states_path),
    f("GAP range file", gap_path),
    f("Outputs saved", length(outputs)),
    vapply(outputs, function(x) paste0(" - ", x), character(1)),
    f("Total elapsed", fmt_elapsed),
    ""
  )

  try({
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
    cat(paste0(log_block, collapse = "\n"), file = log_file, append = TRUE)
  }, silent = TRUE)

  invisible(list(plot = p, png = out_png, stats = out_csv))
}
