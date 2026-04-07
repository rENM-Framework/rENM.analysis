#' Find the weighted centroid for suitability in time series predictions
#'
#' Computes the longitude/latitude weighted centroid of a species prediction
#' raster for a given 5-year time step. If year is omitted, all standard
#' time steps are scanned and processed where inputs exist.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Input file:
#'   \code{<rENM_project_dir()>/runs/<alpha_code>/TimeSeries/<year>/model/}
#'   \code{<alpha_code>-<year>-Prediction.(tif|asc)}
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Output file:
#'   \code{<rENM_project_dir()>/runs/<alpha_code>/TimeSeries/<year>/model/}
#'   \code{<alpha_code>-<year>-Centroid.csv}
#'   \item Log file (appended):
#'   \code{<rENM_project_dir()>/runs/<alpha_code>/}
#'   \code{_log.txt}
#' }
#'
#' \strong{Processing}
#' The raster is reprojected to WGS84 (EPSG:4326) if necessary so that
#' returned and saved coordinates are in lon/lat degrees. The centroid is
#' computed as a value-weighted average of cell-center coordinates.
#'
#' \strong{Weights}
#' All finite raster cell values are used as weights. If you prefer to
#' ignore negative values, uncomment the indicated line to clamp negatives
#' to zero.
#'
#' \strong{Logging format}
#' The log block mirrors the eBird style, e.g.:
#' \preformatted{
#' ------------------------------------------------------------------------
#' Processing summary (find_weighted_centroid)
#' Timestamp:          2025-10-03 10:26:25 MDT
#' Alpha code:         CASP
#' Years processed:    9
#' Completed:          9
#' Failed:             0
#' Total elapsed:      00:00:04.85 (4.85 s)
#' Outputs saved:      centroid_csv
#' Per-year status:
#'   1980  status=ok file=tif cells=12345 wsum=456.78   lon=-100.123456 lat=40.654321 notes=reprojected
#'   ...
#' }
#'
#' @param alpha_code Character. Four-letter banding code for the species
#'   (e.g., "AMRO"). Must be length 4.
#' @param year Integer, NULL. A year within \{1980, 1985, ..., 2020\}.
#'   If NULL (default), the function attempts all of those years and
#'   processes any that have an available input raster.
#'
#' @return
#' Data frame with columns:
#' \itemize{
#'   \item alpha_code: Species code passed in.
#'   \item year: Year processed.
#'   \item lon: Weighted centroid longitude (WGS84) or NA on failure.
#'   \item lat: Weighted centroid latitude (WGS84) or NA on failure.
#'   \item file: Path to the input raster used (if any).
#'   \item status: Processing status ("ok" on success, otherwise a tag).
#'   \item file_type: "tif" or "asc" when available.
#'   \item n_cells_used: Number of non-NA cells used.
#'   \item weight_sum: Sum of weights used after cleaning.
#'   \item notes: Short note (e.g., "reprojected to EPSG:4326").
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes one CSV per successful year with columns lon, lat.
#'   \item Appends a multi-line summary block to _log.txt.
#' }
#'
#' @importFrom utils write.csv
#'
#' @examples
#' \dontrun{
#' # Single year
#' find_weighted_centroid("AMRO", 2005)
#'
#' # All standard time steps with available inputs
#' find_weighted_centroid("AMRO")
#' }
#'
#' @seealso terra
#'
#' @export
find_weighted_centroid <- function(alpha_code, year = NULL) {
  ## ------------------------------ Input validation ------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || nchar(alpha_code) != 4L) {
    stop("alpha_code must be a 4-letter character string (e.g., 'AMRO').")
  }
  valid_years <- seq(1980, 2020, by = 5)

  do_all <- is.null(year) || (missing(year))
  if (!do_all) {
    if (!is.numeric(year) || length(year) != 1L || !(year %in% valid_years)) {
      stop("year must be one of: ", paste(valid_years, collapse = ", "), ".")
    }
  }

  ## ------------------------------ Dependencies ------------------------------
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is required. Please install it: install.packages('terra')")
  }

  ## ------------------------------ Helpers ------------------------------
  .fmt_hms <- function(sec) {
    h <- floor(sec / 3600); m <- floor((sec %% 3600) / 60); s <- sec - h*3600 - m*60
    sprintf("%02d:%02d:%05.2f", h, m, s)
  }

  ## ------------------------------ Paths & Logging ------------------------------
  project_dir <- rENM_project_dir()
  species_root <- file.path(project_dir, "runs", alpha_code)
  if (!dir.exists(species_root)) {
    dir.create(species_root, recursive = TRUE, showWarnings = FALSE)
  }
  log_path <- file.path(species_root, "_log.txt")

  ## ------------------------------ Per-year worker ------------------------------
  .one_year <- function(y) {
    project_dir <- rENM_project_dir()
    base_dir <- file.path(project_dir, "runs", alpha_code, "TimeSeries", y, "model")
    stem_pred <- sprintf("%s-%s-Prediction", alpha_code, y)
    tif_path  <- file.path(base_dir, paste0(stem_pred, ".tif"))
    asc_path  <- file.path(base_dir, paste0(stem_pred, ".asc"))

    has_tif <- file.exists(tif_path)
    has_asc <- file.exists(asc_path)
    if (!has_tif && !has_asc) {
      return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                        file = NA_character_, status = "missing_input",
                        file_type = NA_character_, n_cells_used = NA_integer_,
                        weight_sum = NA_real_, notes = "-", stringsAsFactors = FALSE))
    }
    cand <- if (has_tif) tif_path else asc_path
    ftype <- if (has_tif) "tif" else "asc"

    r <- try(terra::rast(cand), silent = TRUE)
    if (inherits(r, "try-error")) {
      return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                        file = cand, status = "read_error", file_type = ftype,
                        n_cells_used = NA_integer_, weight_sum = NA_real_,
                        notes = "-", stringsAsFactors = FALSE))
    }

    reproj_note <- "-"
    if (!terra::is.lonlat(r)) {
      r <- try(terra::project(r, "EPSG:4326"), silent = TRUE)
      if (inherits(r, "try-error")) {
        return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                          file = cand, status = "reproject_error", file_type = ftype,
                          n_cells_used = NA_integer_, weight_sum = NA_real_,
                          notes = "-", stringsAsFactors = FALSE))
      }
      reproj_note <- "reprojected"
    }

    df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
    if (!nrow(df)) {
      return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                        file = cand, status = "no_data", file_type = ftype,
                        n_cells_used = 0L, weight_sum = NA_real_,
                        notes = reproj_note, stringsAsFactors = FALSE))
    }

    val_col <- setdiff(names(df), c("x", "y"))
    if (length(val_col) != 1L) {
      return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                        file = cand, status = "unexpected_layers", file_type = ftype,
                        n_cells_used = NA_integer_, weight_sum = NA_real_,
                        notes = reproj_note, stringsAsFactors = FALSE))
    }

    w <- df[[val_col]]
    w[!is.finite(w)] <- NA

    # w[w < 0] <- 0

    keep <- !is.na(w)
    if (!any(keep)) {
      return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                        file = cand, status = "all_weights_na", file_type = ftype,
                        n_cells_used = 0L, weight_sum = NA_real_,
                        notes = reproj_note, stringsAsFactors = FALSE))
    }
    w  <- w[keep]
    xs <- df$x[keep]
    ys <- df$y[keep]

    wsum <- sum(w, na.rm = TRUE)
    if (!is.finite(wsum) || wsum == 0) {
      return(data.frame(alpha_code = alpha_code, year = y, lon = NA_real_, lat = NA_real_,
                        file = cand, status = "zero_weight_sum", file_type = ftype,
                        n_cells_used = length(w), weight_sum = wsum,
                        notes = reproj_note, stringsAsFactors = FALSE))
    }

    lon <- sum(xs * w) / wsum
    lat <- sum(ys * w) / wsum

    out_dir <- base_dir
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    out_csv <- file.path(out_dir, sprintf("%s-%s-Centroid.csv", alpha_code, y))
    utils::write.csv(data.frame(lon = lon, lat = lat), out_csv, row.names = FALSE)

    data.frame(alpha_code = alpha_code, year = y, lon = lon, lat = lat,
               file = cand, status = "ok", file_type = ftype,
               n_cells_used = length(w), weight_sum = wsum,
               notes = reproj_note, stringsAsFactors = FALSE)
  }

  years <- if (do_all) valid_years else year
  t0 <- proc.time()[["elapsed"]]

  results <- do.call(rbind, lapply(years, .one_year))

  elapsed <- proc.time()[["elapsed"]] - t0
  elapsed_str <- .fmt_hms(elapsed)

  n_processed <- length(years)
  n_ok        <- sum(results$status == "ok", na.rm = TRUE)
  n_fail      <- sum(results$status != "ok", na.rm = TRUE)

  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  cat("\n", file = log_path, append = TRUE)

  header_line <- "------------------------------------------------------------------------"
  block <- c(
    header_line,
    "Processing summary (find_weighted_centroid)",
    sprintf("Timestamp:          %s", ts),
    sprintf("Alpha code:         %s", alpha_code),
    sprintf("Years processed:    %d", n_processed),
    sprintf("Completed:          %d", n_ok),
    sprintf("Failed:             %d", n_fail),
    sprintf("Total elapsed:      %s (%.2f s)", elapsed_str, elapsed),
    "Outputs saved:      centroid_csv",
    "Per-year status:"
  )

  results <- results[order(results$year), , drop = FALSE]
  per_year <- vapply(seq_len(nrow(results)), function(i) {
    row <- results[i, ]
    yr   <- sprintf("%4d", row$year)
    st   <- paste0("status=", row$status)
    ft   <- paste0("file=", ifelse(is.na(row$file_type), "-", row$file_type))
    cells<- paste0("cells=", ifelse(is.na(row$n_cells_used), "-", format(row$n_cells_used, big.mark = "")))
    wsum <- paste0("wsum=", ifelse(is.na(row$weight_sum), "-", sprintf("%.6g", row$weight_sum)))
    lon  <- paste0("lon=",  ifelse(is.na(row$lon), "-", sprintf("%.6f", row$lon)))
    lat  <- paste0("lat=",  ifelse(is.na(row$lat), "-", sprintf("%.6f", row$lat)))
    note <- paste0("notes=", ifelse(is.na(row$notes) || row$notes == "", "-", row$notes))
    sprintf("  %s  %s %s %s %s   %s %s %s", yr, st, ft, cells, wsum, lon, lat, note)
  }, character(1))

  block <- c(block, per_year)

  cat(paste0(paste(block, collapse = "\n"), "\n"), file = log_path, append = TRUE)

  if (n_ok > 0) {
    message("Centroids written for years: ",
            paste(results$year[results$status == "ok"], collapse = ", "))
  }
  if (n_fail > 0) {
    warn_rows <- results[results$status != "ok", , drop = FALSE]
    message("Skipped/problem years:\n",
            paste(sprintf("  %s: %s", warn_rows$year, warn_rows$status), collapse = "\n"))
  }

  results
}
