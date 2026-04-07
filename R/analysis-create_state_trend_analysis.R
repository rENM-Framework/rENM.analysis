#' Create per-state suitability trend analysis and map for a species
#'
#' Create a per-state analysis of climatic suitability trends and a
#' publication-ready map showing spatial patterns and state-level metrics.
#' Results include tabular summaries, visualization, and logging.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item U.S. state boundaries:
#'   \code{<project_dir>/data/shapefiles/tl_2012_us_state/}
#'   \code{tl_2012_us_state.shp}
#'   \item GAP range shapefile:
#'   \code{<project_dir>/data/shapefiles/}
#'   \code{b<alpha_code>x_CONUS_Range_2001v1/}
#'   \code{b<alpha_code>x_CONUS_Range_2001v1.shp}
#'   \item Suitability trend raster:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item CSV summary:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend-State-Analysis.csv}
#'   \item PNG map:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend-State-Analysis.png}
#'   \item Processing summary appended to:
#'   \code{<project_dir>/runs/<alpha_code>/}
#'   \code{_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Filter to CONUS-only states and clip the GAP range to CONUS
#'   \item Identify states intersecting the GAP range
#'   \item Compute per-state metrics:
#'   \itemize{
#'     \item \code{GAP.RANGE.AREA} (km^2): GAP area within the state
#'     (vector area in EPSG:5070)
#'     \item \code{GAP.RANGE.PCT} (percent): share of total CONUS GAP
#'     area
#'     \item \code{GAP.RANGE.POS.PCT} (percent): share of the state's GAP
#'     area with positive trend (raster > 0)
#'     \item \code{GAP.RANGE.NEG.PCT} (percent): share with negative
#'     trend (raster < 0)
#'   }
#'   \item Map construction includes:
#'   \itemize{
#'     \item Underlay of trend raster
#'     \item Dashed GAP outline
#'     \item Highlighting intersecting states
#'     \item State labeling
#'     \item Cropping to intersection extent
#'   }
#'   \item Log: append a 72-dash processing summary block with aligned
#'   fields including Timestamp, Alpha code, Raster source,
#'   Total/Valid/Positive/Negative/Zero cells, Outputs saved,
#'   Total elapsed, and Output file
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item State shapefile must include the \code{STUSPS} column
#'   \item Raster and vector data must share compatible CRS or be
#'   transformable
#'   \item GAP range geometry must be valid and non-empty
#' }
#'
#' @param alpha_code Character. Scalar species code (e.g., "CASP") used
#' in file paths, titles, and outputs.
#'
#' @return Invisibly returns a List with:
#' \itemize{
#'   \item plot: ggplot object of the final map
#'   \item table: data.frame of per-state metrics
#'   \item csv: Character path to the saved CSV
#'   \item png: Character path to the saved PNG
#'   \item log: Character path to the appended log file
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes CSV and PNG outputs to disk
#'   \item Appends a processing summary to \code{_log.txt}
#' }
#'
#' @importFrom sf st_read st_make_valid st_crs st_transform st_union st_bbox
#'   st_geometry st_crop st_point_on_surface st_sf st_area st_is_empty
#'   st_collection_extract
#' @importFrom terra rast same.crs vect project ext mask crop cellSize global
#'   ncell ifel as.data.frame
#' @importFrom dplyr bind_rows
#' @importFrom ggplot2 ggplot geom_raster aes scale_fill_gradient2 geom_sf
#'   geom_sf_text coord_sf labs theme_void theme element_rect element_text
#'   ggsave theme_minimal
#' @importFrom grid unit
#' @importFrom rlang .data
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#'   create_state_trend_analysis("CASP")
#' }
#'
#' @export
create_state_trend_analysis <- function(alpha_code) {
  stopifnot(is.character(alpha_code), length(alpha_code) == 1L, nzchar(alpha_code))

  # Ensure required namespaces are present
  for (p in c("sf", "dplyr", "ggplot2", "terra", "grid", "rlang")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required package '", p, "' is not installed.", call. = FALSE)
    }
  }

  start_time <- Sys.time()

  # ---- Paths ---------------------------------------------------------------
  project_dir <- rENM_project_dir()

  states_path <- file.path(
    project_dir, "data", "shapefiles", "tl_2012_us_state", "tl_2012_us_state.shp"
  )

  gap_path <- file.path(
    project_dir, "data", "shapefiles",
    sprintf("b%sx_CONUS_Range_2001v1", alpha_code),
    sprintf("b%sx_CONUS_Range_2001v1.shp", alpha_code)
  )

  trend_path <- file.path(
    project_dir, "runs", alpha_code, "Trends", "suitability",
    sprintf("%s-Suitability-Trend.tif", alpha_code)
  )

  out_dir <- file.path(
    project_dir, "runs", alpha_code, "Trends", "suitability"
  )
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  out_csv <- file.path(
    out_dir,
    sprintf("%s-Suitability-Trend-State-Analysis.csv", alpha_code)
  )

  out_png <- file.path(
    out_dir,
    sprintf("%s-Suitability-Trend-State-Analysis.png", alpha_code)
  )

  out_log <- file.path(
    project_dir, "runs", alpha_code, "_log.txt"
  )

  # ---- Read vectors --------------------------------------------------------
  message("==> [", alpha_code, "] Reading vector data...")
  states <- sf::st_read(states_path, quiet = TRUE)
  states <- sf::st_make_valid(states)

  gap <- sf::st_read(gap_path, quiet = TRUE)
  gap <- sf::st_make_valid(gap)

  if (!identical(sf::st_crs(states), sf::st_crs(gap))) {
    gap <- sf::st_transform(gap, sf::st_crs(states))
  }

  trySuppressWarnings <- function(expr) {
    suppressWarnings(
      tryCatch(expr, error = function(e) NULL)
    )
  }

  safe_intersection <- function(a, b) {
    out <- trySuppressWarnings(sf::st_intersection(a, b))
    if (inherits(out, "sfc")) {
      out <- sf::st_sf(geometry = out, crs = sf::st_crs(a))
    }
    if (is.null(out)) {
      return(sf::st_sf(
        geometry = sf::st_sfc(),
        crs      = sf::st_crs(a)
      ))
    }
    out <- suppressWarnings(
      tryCatch(
        sf::st_collection_extract(out, "POLYGON"),
        error = function(e) out
      )
    )
    if (inherits(out, "sf")) {
      out <- out[!sf::st_is_empty(out), , drop = FALSE]
    }
    out
  }

  is_empty_sf <- function(x) {
    if (is.null(x)) {
      TRUE
    } else if (inherits(x, "sf")) {
      nrow(x) == 0L || all(sf::st_is_empty(x))
    } else {
      FALSE
    }
  }

  # ---- Filter to CONUS -----------------------------------------------------
  message("==> Filtering to CONUS states...")
  non_conus <- c("AK", "HI", "PR", "GU", "VI", "AS", "MP", "UM")
  state_id_col <- "STUSPS"

  if (!state_id_col %in% names(states)) {
    stop("Expected 'STUSPS' column not found in states shapefile.", call. = FALSE)
  }

  states_conus <- states[!(states[[state_id_col]] %in% non_conus), , drop = FALSE]

  # ---- Clip GAP to CONUS ---------------------------------------------------
  message("==> Clipping GAP to CONUS...")
  gap_conus <- safe_intersection(
    sf::st_union(gap),
    sf::st_union(states_conus)
  )

  idx <- sf::st_intersects(states_conus, gap_conus, sparse = FALSE)[, 1]
  states_in_gap <- states_conus[idx, , drop = FALSE]

  conus_bbox <- sf::st_bbox(
    sf::st_union(sf::st_geometry(states_in_gap), sf::st_geometry(gap_conus))
  )
  states_crop <- sf::st_crop(states_conus, conus_bbox)

  # ---- Read trend raster ---------------------------------------------------
  message("==> Reading trend raster...")
  trend_raster <- terra::rast(trend_path)

  if (!terra::same.crs(trend_raster, terra::vect(states_conus))) {
    trend_raster <- terra::project(
      trend_raster,
      terra::crs(terra::vect(states_conus))
    )
  }

  r_ext <- terra::ext(
    conus_bbox["xmin"],
    conus_bbox["xmax"],
    conus_bbox["ymin"],
    conus_bbox["ymax"]
  )
  trend_crop <- terra::mask(
    terra::crop(trend_raster, r_ext),
    terra::vect(states_conus)
  )

  # ---- Labels --------------------------------------------------------------
  lab_vals <- as.character(states_in_gap[[state_id_col]])
  lab_pts_geom <- suppressWarnings(
    sf::st_point_on_surface(sf::st_geometry(states_in_gap))
  )
  label_points <- sf::st_sf(
    label   = lab_vals,
    geometry = lab_pts_geom,
    crs     = sf::st_crs(states_in_gap)
  )

  # ---- Area + trend by state ----------------------------------------------
  message("==> Computing per-state metrics...")
  eq_crs <- sf::st_crs("EPSG:5070")
  states_eq <- sf::st_transform(states_in_gap, eq_crs)
  gap_eq    <- sf::st_transform(gap_conus, eq_crs)

  total_gap_m2 <- as.numeric(sf::st_area(gap_eq))

  trend_eq <- terra::project(
    trend_crop,
    as.character(eq_crs$wkt)
  )
  cell_area <- terra::cellSize(trend_eq, unit = "m")

  results_list <- lapply(seq_len(nrow(states_eq)), function(i) {
    st_row_eq <- states_eq[i, , drop = FALSE]
    state_id  <- lab_vals[i]

    gap_state_eq <- safe_intersection(gap_eq, st_row_eq)

    if (is_empty_sf(gap_state_eq)) {
      return(data.frame(
        STATE             = state_id,
        GAP.RANGE.AREA    = 0,
        GAP.RANGE.PCT     = 0,
        GAP.RANGE.POS.PCT = 0,
        GAP.RANGE.NEG.PCT = 0,
        stringsAsFactors  = FALSE
      ))
    }

    area_gap_m2 <- sum(as.numeric(sf::st_area(gap_state_eq)), na.rm = TRUE)
    pct_total   <- if (total_gap_m2 > 0) (area_gap_m2 / total_gap_m2) * 100 else 0

    gap_spat <- terra::vect(gap_state_eq)
    r_state  <- terra::mask(terra::crop(trend_eq, gap_spat), gap_spat)
    ca_state <- terra::mask(terra::crop(cell_area, gap_spat), gap_spat)

    area_total_m2 <- trySuppressWarnings(
      terra::global(ca_state, "sum", na.rm = TRUE)[1, 1]
    )

    if (is.null(area_total_m2) || is.na(area_total_m2) || area_total_m2 == 0) {
      pos_pct <- 0
      neg_pct <- 0
    } else {
      pos_mask <- r_state > 0
      neg_mask <- r_state < 0

      area_pos <- trySuppressWarnings(
        terra::global(ca_state * pos_mask, "sum", na.rm = TRUE)[1, 1]
      )
      if (is.na(area_pos)) area_pos <- 0

      area_neg <- trySuppressWarnings(
        terra::global(ca_state * neg_mask, "sum", na.rm = TRUE)[1, 1]
      )
      if (is.na(area_neg)) area_neg <- 0

      pos_pct <- (area_pos / area_total_m2) * 100
      neg_pct <- (area_neg / area_total_m2) * 100
    }

    data.frame(
      STATE             = state_id,
      GAP.RANGE.AREA    = round(area_gap_m2 / 1e6, 3),
      GAP.RANGE.PCT     = round(pct_total, 3),
      GAP.RANGE.POS.PCT = round(pos_pct, 3),
      GAP.RANGE.NEG.PCT = round(neg_pct, 3),
      stringsAsFactors  = FALSE
    )
  })

  results <- dplyr::bind_rows(results_list)
  utils::write.csv(results, out_csv, row.names = FALSE)

  # ---- Map -----------------------------------------------------------------
  message("==> Building map and saving PNG...")
  trend_df <- terra::as.data.frame(trend_crop, xy = TRUE, na.rm = TRUE)
  names(trend_df)[3] <- "value"

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data    = trend_df,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value)
    ) +
    ggplot2::scale_fill_gradient2(
      name      = "Trend",
      low       = "#D64241",
      mid       = "white",
      high      = "#58B440",
      midpoint  = 0
    ) +
    ggplot2::geom_sf(
      data   = states_crop,
      fill   = NA,
      color  = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_sf(
      data   = states_in_gap,
      fill   = NA,
      color  = "black",
      linewidth = 0.4
    ) +
    ggplot2::geom_sf(
      data   = gap_conus,
      fill   = NA,
      color  = "black",
      linewidth = 0.7,
      linetype = "dashed"
    ) +
    ggplot2::geom_sf_text(
      data   = label_points,
      mapping = ggplot2::aes(label = .data$label),
      size   = 6,
      color  = "black",
      fontface = "bold"
    ) +
    ggplot2::coord_sf(
      xlim   = c(conus_bbox["xmin"], conus_bbox["xmax"]),
      ylim   = c(conus_bbox["ymin"], conus_bbox["ymax"]),
      expand = FALSE
    ) +
    ggplot2::labs(
      title    = sprintf("%s Climatic Suitability Trend", alpha_code),
      subtitle = sprintf("States encompassing %s GAP Range", alpha_code),
      x        = "Longitude",
      y        = "Latitude"
    ) +
    ggplot2::theme_void() +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5)) +
    ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border     = ggplot2::element_rect(
        fill   = NA,
        color  = "darkgray",
        linewidth = 0.8
      ),
      plot.title       = ggplot2::element_text(size = 28, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = 18, face = "plain"),
      legend.text      = ggplot2::element_text(size = 16),
      legend.title     = ggplot2::element_text(size = 16, face = "plain"),
      legend.key.width  = grid::unit(0.3, "cm"),
      legend.key.height = grid::unit(0.8, "cm"),
      axis.title       = ggplot2::element_text(size = 24, face = "plain"),
      axis.text        = ggplot2::element_text(size = 18),
      legend.position  = "right"
    )

  ggplot2::ggsave(
    filename  = out_png,
    plot      = p,
    width     = 1000,
    height    = 750,
    units     = "px",
    dpi       = 72,
    pointsize = 15,
    bg        = "white"
  )

  # ---- Log -----------------------------------------------------------------
  message("==> Appending processing summary to log...")
  dir.create(dirname(out_log), recursive = TRUE, showWarnings = FALSE)

  gap_for_counts <- sf::st_transform(gap_conus, sf::st_crs(trend_crop))
  r_crop_to_gap  <- terra::crop(trend_crop, terra::vect(gap_for_counts))
  r_mask         <- terra::mask(r_crop_to_gap, terra::vect(gap_for_counts))
  lyr            <- r_mask[[1]]

  total_cells <- terra::ncell(r_crop_to_gap[[1]])
  valid_cells <- as.numeric(
    terra::global(terra::ifel(!is.na(lyr), 1, 0), "sum", na.rm = TRUE)[1, 1]
  )
  pos_cells <- as.numeric(
    terra::global(terra::ifel(lyr > 0, 1, 0), "sum", na.rm = TRUE)[1, 1]
  )
  neg_cells <- as.numeric(
    terra::global(terra::ifel(lyr < 0, 1, 0), "sum", na.rm = TRUE)[1, 1]
  )
  zer_cells <- as.numeric(
    terra::global(terra::ifel(lyr == 0, 1, 0), "sum", na.rm = TRUE)[1, 1]
  )

  elapsed <- difftime(Sys.time(), start_time, units = "secs")
  sep     <- paste(rep("-", 72L), collapse = "")

  log_block <- paste0(
    "\n", sep, "\nProcessing summary (create_state_trend_analysis)\n",
    sprintf("%-12s %s\n", "Timestamp:",    format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    sprintf("%-12s %s\n", "Alpha code:",   alpha_code),
    sprintf("%-12s %s\n", "Raster source:", trend_path),
    sprintf("%-12s %s\n", "Total cells:",   format(total_cells, big.mark = ",")),
    sprintf("%-12s %s\n", "Valid cells:",   format(valid_cells, big.mark = ",")),
    sprintf("%-12s %s\n", "Positive cells:", format(pos_cells, big.mark = ",")),
    sprintf("%-12s %s\n", "Negative cells:", format(neg_cells, big.mark = ",")),
    sprintf("%-12s %s\n", "Zero cells:",     format(zer_cells, big.mark = ",")),
    sprintf(
      "%-12s %s\n",
      "Outputs saved:",
      paste(
        basename(out_csv),
        basename(out_png),
        sep = " | "
      )
    ),
    sprintf(
      "%-12s %s seconds\n",
      "Total elapsed:",
      round(as.numeric(elapsed), 2)
    ),
    sprintf("%-12s %s", "Output file:", out_csv)
  )

  try({
    con <- file(out_log, open = "a", encoding = "UTF-8")
    writeLines(log_block, con = con)
    close(con)
    message("==> Log updated: ", out_log)
  }, silent = TRUE)

  message("==> Done.")

  invisible(list(
    plot = p,
    table = results,
    csv = out_csv,
    png = out_png,
    log = out_log
  ))
}
