#' Compute GAP range change percentages for a species trend raster
#'
#' Compute percentages of positive, negative, and zero suitability trend
#' values within a species GAP range using a masked raster. Results are
#' summarized, written to disk, and logged.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Species table:
#'   \code{<project_dir>/data/}
#'   \code{_species.csv}
#'   \item GAP.RANGE shapefile resolved from entries such as:
#'   \code{<project_dir>/data/shapefiles/<name>/}
#'   \code{<name>.shp}
#'   \item Suitability trend raster:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item CSV summary written to:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend-GAP-Range-Percentages.csv}
#'   \item Processing summary appended to:
#'   \code{<project_dir>/runs/<alpha_code>/}
#'   \code{_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Raster is masked by the GAP range polygon
#'   \item Percentages of positive (>0), negative (<0), and zero (=0)
#'   values are computed among non-NA cells
#'   \item Denominator for percentages is the number of non-NA cells
#'   inside the GAP range
#'   \item Zeros are reported separately
#'   \item Area is computed using EPSG:6933 (World Cylindrical Equal
#'   Area)
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item GAP.RANGE must resolve to a valid shapefile path
#'   \item Raster and vector data must be readable by
#'   \code{terra::rast} and \code{sf::read_sf}
#' }
#'
#' For a given species \code{alpha_code}, this function:
#' \itemize{
#'   \item Reads \code{<project_dir>/data/_species.csv} to locate the
#'   GAP.RANGE shapefile reference
#'   \item Resolves the GAP.RANGE shapefile path, including common
#'   layouts such as
#'   \code{<project_dir>/data/shapefiles/<name>/<name>.shp}
#'   \item Reads the species suitability trend raster at
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend.tif}
#'   \item Masks the raster by the GAP range polygon and computes the
#'   percentage of positive (>0), negative (<0), and zero (=0) cells
#'   among non-NA cells
#'   \item Computes total GAP range area (km^2) in an equal-area CRS
#'   (EPSG:6933)
#'   \item Writes a CSV summary to
#'   \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'   \code{<alpha_code>-Suitability-Trend-GAP-Range-Percentages.csv}
#'   \item Appends a processing summary to
#'   \code{<project_dir>/runs/<alpha_code>/}
#'   \code{_log.txt}
#'   using the standard format (72-dash separator; labeled header;
#'   aligned fields; no "Notes")
#' }
#'
#' @param alpha_code Character. Species ALPHA.CODE used to resolve
#' inputs and outputs.
#'
#' @return Invisibly returns a tibble with:
#' \itemize{
#'   \item counts: positive, negative, zero, and NA cells
#'   \item percentages: percent positive, negative, and zero
#'   \item area: total GAP range area in km^2
#'   \item paths: raster and GAP shapefile paths
#' }
#'
#' Side effects:
#' \itemize{
#'   \item CSV written to disk
#'   \item Processing summary appended to \code{_log.txt}
#' }
#'
#' @importFrom readr read_csv write_csv
#' @importFrom dplyr filter slice_head
#' @importFrom sf read_sf st_make_valid st_union st_transform st_area
#' @importFrom terra rast crop mask values crs vect
#' @importFrom tibble tibble
#' @importFrom tools file_path_sans_ext
#'
#' @examples
#' \dontrun{
#'   find_range_change_percentages("CASP")
#' }
#'
#' @export
find_range_change_percentages <- function(alpha_code) {
  # ---------------------------- Setup & Timing ----------------------------
  t_start <- Sys.time()
  step_i  <- 0L
  bump <- function(label) {
    step_i <<- step_i + 1L
    message(sprintf("[Step %02d] %s", step_i, label))
  }

  for (p in c("readr", "dplyr", "sf", "terra", "tibble")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required package '", p, "' is not installed.", call. = FALSE)
    }
  }

  # ---------------------------- Paths & Inputs ----------------------------
  bump("Validating alpha_code and core paths")
  if (missing(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }

  project_dir <- rENM_project_dir()

  species_csv <- file.path(project_dir, "data", "_species.csv")
  if (!file.exists(species_csv)) {
    stop("Species CSV not found at: ", species_csv, call. = FALSE)
  }

  runs_root <- file.path(project_dir, "runs", alpha_code)
  out_dir   <- file.path(runs_root, "Trends", "suitability")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  out_csv <- file.path(
    out_dir,
    paste0(alpha_code, "-Suitability-Trend-GAP-Range-Percentages.csv")
  )
  out_log <- file.path(runs_root, "_log.txt")

  trend_raster_path <- file.path(
    runs_root, "Trends", "suitability",
    paste0(alpha_code, "-Suitability-Trend.tif")
  )
  if (!file.exists(trend_raster_path)) {
    stop("Trend raster not found at: ", trend_raster_path, call. = FALSE)
  }

  # ---------------------------- Read species table -------------------------
  bump("Reading _species.csv and locating GAP.RANGE entry")
  sp <- readr::read_csv(species_csv, show_col_types = FALSE)

  alpha_col <- intersect(names(sp), c("ALPHA.CODE", "alpha_code", "Alpha.code"))[1]
  gap_col   <- intersect(names(sp), c("GAP.RANGE", "gap_range", "GAP_RANGE", "Gap.Range"))[1]
  if (is.na(alpha_col)) {
    stop("Could not find an ALPHA.CODE column in _species.csv", call. = FALSE)
  }
  if (is.na(gap_col)) {
    stop("Could not find a GAP.RANGE column in _species.csv", call. = FALSE)
  }

  row <- sp |>
    dplyr::filter(.data[[alpha_col]] == alpha_code) |>
    dplyr::slice_head(n = 1)
  if (nrow(row) == 0L) {
    stop("ALPHA.CODE not found in _species.csv: ", alpha_code, call. = FALSE)
  }

  gap_entry <- row[[gap_col]][1]
  if (is.na(gap_entry) || !nzchar(gap_entry)) {
    stop("GAP.RANGE entry is missing for ALPHA.CODE: ", alpha_code, call. = FALSE)
  }

  # ---------------------------- Resolve GAP shapefile path -----------------
  bump("Resolving GAP.RANGE shapefile path")
  expand <- function(p) p
  entry  <- gap_entry
  entry_no_ext <- sub("\\.shp$", "", entry, ignore.case = TRUE)
  base_data   <- file.path(project_dir, "data")
  base_shapes <- file.path(project_dir, "data", "shapefiles")

  candidates <- c(
    expand(entry),
    file.path(expand(entry), paste0(basename(entry_no_ext), ".shp")),
    file.path(base_shapes, entry_no_ext, paste0(basename(entry_no_ext), ".shp")),
    file.path(base_shapes, entry_no_ext),
    file.path(base_data, paste0(entry_no_ext, ".shp")),
    file.path(base_data, entry_no_ext, paste0(basename(entry_no_ext), ".shp"))
  )

  pick_shp_in_dir <- function(d) {
    if (!dir.exists(d)) return(NA_character_)
    shp <- list.files(d, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
    if (length(shp) == 0L) return(NA_character_)
    preferred <- shp[basename(tools::file_path_sans_ext(shp)) == basename(d)]
    if (length(preferred)) return(preferred[1])
    shp[1]
  }

  resolved <- NA_character_
  tried    <- character(0)
  for (cand in candidates) {
    if (is.na(cand)) next
    if (dir.exists(cand)) {
      shp_inside <- pick_shp_in_dir(cand)
      tried <- c(tried, cand)
      if (!is.na(shp_inside) && file.exists(shp_inside)) {
        resolved <- shp_inside
        break
      }
    } else {
      tried <- c(tried, cand)
      if (file.exists(cand) && grepl("\\.shp$", cand, ignore.case = TRUE)) {
        resolved <- cand
        break
      }
    }
  }
  if (is.na(resolved)) {
    stop(
      "GAP.RANGE shapefile not found for entry '", gap_entry, "'.\n",
      "Tried:\n - ", paste(tried, collapse = "\n - "), "\n",
      "Hint: If your layout is like <project_dir>/data/shapefiles/<name>/<name>.shp,\n",
      "store '<name>' (no extension) in GAP.RANGE; this function searches there.",
      call. = FALSE
    )
  }
  gap_path <- resolved
  message("Resolved GAP shapefile: ", gap_path)

  # ---------------------------- Read data ----------------------------------
  bump("Reading GAP range and suitability trend raster")
  rng_sf <- tryCatch(
    sf::read_sf(gap_path, quiet = TRUE),
    error = function(e) {
      stop("Failed to read GAP.RANGE shapefile at ", gap_path, ": ",
           conditionMessage(e), call. = FALSE)
    }
  )
  if (nrow(rng_sf) == 0L) {
    stop("GAP.RANGE shapefile has no features: ", gap_path, call. = FALSE)
  }

  r <- tryCatch(
    terra::rast(trend_raster_path),
    error = function(e) {
      stop("Failed to read trend raster: ", conditionMessage(e), call. = FALSE)
    }
  )

  # ---------------------------- Prepare geometry ---------------------------
  bump("Validating geometry and aligning CRS")
  rng_sf <- sf::st_make_valid(rng_sf)  # repair self-intersections, etc.
  rng_sf <- sf::st_union(rng_sf)       # dissolve to single multipart polygon
  if (!is.na(terra::crs(r))) {
    rng_to_rast <- sf::st_transform(rng_sf, terra::crs(r))
  } else {
    warning(
      "Raster has no CRS. Proceeding without reprojection; results may be inaccurate.",
      call. = FALSE
    )
    rng_to_rast <- rng_sf
  }
  rng_v <- terra::vect(rng_to_rast)

  # ---------------------------- Masking & statistics -----------------------
  bump("Cropping and masking raster to GAP range")
  r_crop <- terra::crop(r, rng_v, snap = "out")
  r_mask <- terra::mask(r_crop, rng_v)

  bump("Computing cell statistics inside GAP range")
  vals <- terra::values(r_mask, na.rm = FALSE)

  n_total <- sum(!is.na(vals))                # denominator for percentages
  if (n_total == 0L) {
    stop("No non-NA raster cells found within the GAP.RANGE.", call. = FALSE)
  }

  n_pos  <- sum(vals > 0, na.rm = TRUE)       # positive change cells
  n_neg  <- sum(vals < 0, na.rm = TRUE)       # negative change cells
  n_zero <- sum(vals == 0, na.rm = TRUE)      # zero change cells
  n_na   <- sum(is.na(vals))                  # NA within mask (for transparency)

  pct_pos  <- 100 * n_pos  / n_total
  pct_neg  <- 100 * n_neg  / n_total
  pct_zero <- 100 * n_zero / n_total
  pct_sum  <- pct_pos + pct_neg + pct_zero    # should be ~100 (rounding)

  # ---------------------------- Area computation --------------------------
  bump("Computing GAP range area (km^2) in equal-area CRS (EPSG:6933)")
  rng_eq   <- sf::st_transform(rng_sf, 6933)  # World Cylindrical Equal Area
  area_m2  <- as.numeric(sf::st_area(rng_eq))
  total_area_km2 <- sum(area_m2, na.rm = TRUE) / 1e6

  # ---------------------------- Assemble & write CSV ----------------------
  bump("Assembling results and writing CSV")
  res <- tibble::tibble(
    alpha_code            = alpha_code,
    total_cells           = n_total,
    positive_cells        = n_pos,
    negative_cells        = n_neg,
    zero_cells            = n_zero,
    na_cells              = n_na,
    percent_positive      = round(pct_pos,  6),
    percent_negative      = round(pct_neg,  6),
    percent_zero          = round(pct_zero, 6),
    gap_range_area_km2    = round(total_area_km2, 6),
    raster_path           = trend_raster_path,
    gap_range_path        = gap_path,
    computed_at_utc       = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

  readr::write_csv(res, out_csv)
  message("Wrote CSV: ", out_csv)

  # ---------------------------- Append processing summary log -------------
  bump("Appending processing summary to runs root _log.txt")
  t_end   <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))

  # 72-dash separator per standard
  sep72  <- paste(rep("-", 72), collapse = "")
  header <- sprintf("Processing summary (%s)", "find_range_change_percentages")

  # Ensure the runs_root exists (should already)
  if (!dir.exists(runs_root)) {
    dir.create(runs_root, recursive = TRUE, showWarnings = FALSE)
  }

  # If the log already has content, prefix a blank line before the dashed separator
  size_before <- if (file.exists(out_log)) {
    suppressWarnings(file.info(out_log)$size)
  } else {
    0L
  }
  pre_blank <- if (is.finite(size_before) && size_before > 0) "" else character(0)

  # Format numeric fields for the log (stable width, readable)
  fmt_pct  <- function(x) sprintf("%.6f %%", x)
  fmt_area <- function(x) sprintf("%.6f", x)

  log_lines <- c(
    pre_blank,  # blank line before dashes if appending to existing log
    sep72,
    header,
    sprintf("%-18s %s", "Timestamp:",   format(Sys.time(), tz = "UTC", usetz = TRUE)),
    sprintf("%-18s %s", "Alpha code:",  alpha_code),
    sprintf("%-18s %s", "Raster source:", trend_raster_path),
    sprintf("%-18s %d", "Total cells:", n_total),
    sprintf("%-18s %d", "Valid cells:", n_total),  # denominator == non-NA
    sprintf("%-18s %d", "Positive cells:", n_pos),
    sprintf("%-18s %d", "Negative cells:", n_neg),
    sprintf("%-18s %d", "Zero cells:",    n_zero),
    sprintf("%-18s %s", "Percent positive:", fmt_pct(pct_pos)),
    sprintf("%-18s %s", "Percent negative:", fmt_pct(pct_neg)),
    sprintf("%-18s %s", "Percent zero:",     fmt_pct(pct_zero)),
    sprintf("%-18s %s", "Percent sum:",      fmt_pct(pct_sum)),
    sprintf("%-18s %s", "GAP area (km^2):",  fmt_area(total_area_km2)),
    sprintf("%-18s %s", "Outputs saved:", "1 CSV"),
    sprintf("%-18s %.3f sec", "Total elapsed:", elapsed),
    sprintf("%-18s %s", "Output file:", out_csv),
    ""  # final newline
  )

  # Robust append with explicit connection; soft verification
  wrote_ok <- TRUE
  tryCatch({
    con <- file(out_log, open = "a", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(log_lines, con, sep = "\n")
    flush(con)
  }, error = function(e) {
    wrote_ok <<- FALSE
    stop("Failed to append to log file at ", out_log, " : ",
         conditionMessage(e), call. = FALSE)
  })

  if (wrote_ok) {
    size_after <- if (file.exists(out_log)) {
      suppressWarnings(file.info(out_log)$size)
    } else {
      NA_integer_
    }
    if (is.finite(size_before) && is.finite(size_after)) {
      delta <- as.integer(size_after - size_before)
      if (delta > 0L) {
        message("Appended log: ", out_log, " (", delta, " bytes)")
      } else {
        warning(
          "Log append completed but file size did not increase; keeping as success. Path: ",
          out_log,
          call. = FALSE
        )
      }
    } else {
      message("Appended log: ", out_log)
    }
  }

  bump("Done")
  invisible(res)
}
