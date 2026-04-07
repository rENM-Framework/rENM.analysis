#' Find positive and negative trend percentages
#'
#' Computes cell-wise trend sign statistics from a suitability trend raster
#' and summarizes proportions of positive, negative, and zero values across
#' the full modeled extent.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Species identifier provided as \code{alpha_code}.
#'   \item A suitability trend raster stored in the project directory.
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item CSV summary file written to:
#'         \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'         \code{<alpha_code>-Suitability-Trend-Percentages.csv}
#'   \item Processing summary appended to:
#'         \code{<project_dir>/runs/<alpha_code>/}
#'         \code{_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Reads the suitability trend raster, preferring GeoTIFF and
#'         falling back to ASCII grid:
#'         \code{<project_dir>/runs/<alpha_code>/Trends/suitability/}
#'         \code{<alpha_code>-Suitability-Trend.(tif|asc)}
#'   \item Computes counts and percentages of positive (>0), negative (<0),
#'         and zero (=0) cells over the full raster extent.
#'   \item Excludes NA cells from percentage denominators.
#'   \item Computes extent area (km^2) using \code{terra::cellSize}
#'         (in m^2), converted to km^2.
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Raster must contain numeric trend values.
#'   \item Raster must exist in either GeoTIFF or ASCII grid format.
#' }
#'
#' Timestamp in CSV and log is UTC. All file paths are constructed
#' relative to the rENM project directory returned by
#' \code{rENM_project_dir()}, ensuring CRAN-compliant path handling.
#'
#' @param alpha_code Character. Species ALPHA.CODE (for example, "CASP").
#'
#' @return
#' Tibble. A single-row table with the following fields:
#' \itemize{
#'   \item \code{alpha_code}: Species identifier.
#'   \item \code{total_cells}: Total raster cell count.
#'   \item \code{valid_cells}: Non-NA cell count.
#'   \item \code{positive_cells}: Count of cells > 0.
#'   \item \code{negative_cells}: Count of cells < 0.
#'   \item \code{zero_cells}: Count of cells == 0.
#'   \item \code{percent_positive}: Percentage of positive cells.
#'   \item \code{percent_negative}: Percentage of negative cells.
#'   \item \code{percent_zero}: Percentage of zero cells.
#'   \item \code{percent_sum}: Sum of percentages.
#'   \item \code{extent_area_km2}: Total raster extent area in km^2.
#'   \item \code{raster_path}: Source raster file path.
#'   \item \code{computed_at_utc}: UTC timestamp of computation.
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes CSV summary file to disk.
#'   \item Appends processing summary to \code{_log.txt}.
#' }
#'
#' @importFrom readr write_csv
#' @importFrom terra rast ncell global cellSize values
#' @importFrom tibble tibble
#'
#' @examples
#' \dontrun{
#'   find_trend_percentages("CASP")
#' }
#'
#' @export
find_trend_percentages <- function(alpha_code) {
  # ---------------------------- Setup & Timing ----------------------------
  t_start <- Sys.time()
  step_i <- 0L
  bump <- function(label) {
    step_i <<- step_i + 1L
    message(sprintf("[Step %02d] %s", step_i, label))
  }

  # Ensure required packages are available (no library() in package code)
  for (p in c("terra", "readr", "tibble")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Required package '", p, "' is not installed.", call. = FALSE)
    }
  }

  # ---------------------------- Validate inputs ----------------------------
  bump("Validating alpha_code and core paths")
  if (missing(alpha_code) || length(alpha_code) != 1L || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar.", call. = FALSE)
  }
  alpha_code <- toupper(alpha_code)

  # ---------------------------- Paths ----------------------------
  project_dir <- rENM_project_dir()
  runs_root <- file.path(project_dir, "runs", alpha_code)
  trend_dir <- file.path(runs_root, "Trends", "suitability")
  if (!dir.exists(trend_dir)) {
    dir.create(trend_dir, recursive = TRUE, showWarnings = FALSE)
  }

  tif_path <- file.path(trend_dir, sprintf("%s-Suitability-Trend.tif", alpha_code))
  asc_path <- file.path(trend_dir, sprintf("%s-Suitability-Trend.asc", alpha_code))

  raster_path <- if (file.exists(tif_path)) {
    tif_path
  } else if (file.exists(asc_path)) {
    asc_path
  } else {
    NULL
  }
  if (is.null(raster_path)) {
    stop("Trend raster not found (.tif or .asc) at: ", trend_dir, call. = FALSE)
  }

  out_csv <- file.path(
    trend_dir,
    sprintf("%s-Suitability-Trend-Percentages.csv", alpha_code)
  )
  out_log <- file.path(runs_root, "_log.txt")

  # ---------------------------- Read raster ----------------------------
  bump("Reading suitability trend raster (full extent)")
  r <- tryCatch(
    terra::rast(raster_path),
    error = function(e) {
      stop("Failed to read trend raster: ", conditionMessage(e), call. = FALSE)
    }
  )
  message("Raster: ", raster_path)

  # ---------------------------- Compute cell stats ----------------------------
  bump("Computing positive, negative, zero, and NA counts over full extent")
  total_cells <- terra::ncell(r)

  # Count NA cells (fast global op)
  na_cells <- as.numeric(terra::global(is.na(r), fun = "sum", na.rm = TRUE)[1, 1])
  valid_cells <- total_cells - na_cells
  if (valid_cells <= 0) {
    stop("All cells are NA; cannot compute percentages.", call. = FALSE)
  }

  # Count signs (terra::global on logical rasters gives cell counts)
  positive_cells <- as.numeric(terra::global(r > 0, fun = "sum", na.rm = TRUE)[1, 1])
  negative_cells <- as.numeric(terra::global(r < 0, fun = "sum", na.rm = TRUE)[1, 1])
  zero_cells     <- as.numeric(terra::global(r == 0, fun = "sum", na.rm = TRUE)[1, 1])

  pct <- function(x, denom) 100 * x / denom
  percent_positive <- pct(positive_cells, valid_cells)
  percent_negative <- pct(negative_cells, valid_cells)
  percent_zero     <- pct(zero_cells,     valid_cells)
  percent_sum      <- percent_positive + percent_negative + percent_zero

  # ---------------------------- Extent area (km^2) ----------------------------
  bump("Computing extent area (km^2) via cell areas")
  extent_area_km2 <- tryCatch({
    cs <- terra::cellSize(r, unit = "m")
    sum(terra::values(cs), na.rm = TRUE) / 1e6
  }, error = function(e) {
    warning(
      "Failed to compute cellSize; extent area set to NA. Error: ",
      conditionMessage(e),
      call. = FALSE
    )
    NA_real_
  })

  # ---------------------------- Assemble & write CSV ------------------------
  bump("Assembling results and writing CSV")
  res <- tibble::tibble(
    alpha_code         = alpha_code,
    total_cells        = total_cells,
    valid_cells        = valid_cells,
    positive_cells     = positive_cells,
    negative_cells     = negative_cells,
    zero_cells         = zero_cells,
    percent_positive   = round(percent_positive, 6),
    percent_negative   = round(percent_negative, 6),
    percent_zero       = round(percent_zero, 6),
    percent_sum        = round(percent_sum, 6),
    extent_area_km2    = if (is.na(extent_area_km2)) NA_real_ else round(extent_area_km2, 6),
    raster_path        = raster_path,
    computed_at_utc    = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

  readr::write_csv(res, out_csv)
  message("Wrote CSV: ", out_csv)

  # ---------------------------- Append processing summary log ---------------
  bump("Appending processing summary to runs root _log.txt")
  t_end   <- Sys.time()
  elapsed <- as.numeric(difftime(t_end, t_start, units = "secs"))

  sep72  <- paste(rep("-", 72), collapse = "")
  header <- sprintf("Processing summary (%s)", "find_trend_percentages")
  src_type <- if (grepl("\\.tif$", raster_path, ignore.case = TRUE)) "tif" else "asc"

  if (!dir.exists(runs_root)) {
    dir.create(runs_root, recursive = TRUE, showWarnings = FALSE)
  }

  size_before <- if (file.exists(out_log)) {
    suppressWarnings(file.info(out_log)$size)
  } else {
    0L
  }
  pre_blank <- if (is.finite(size_before) && size_before > 0) "" else character(0)

  fmt_pct  <- function(x) sprintf("%.6f %%", x)
  fmt_area <- function(x) if (is.na(x)) "NA" else sprintf("%.6f", x)

  log_lines <- c(
    pre_blank,
    sep72,
    header,
    sprintf("%-18s %s", "Timestamp:",     format(Sys.time(), tz = "UTC", usetz = TRUE)),
    sprintf("%-18s %s", "Alpha code:",    alpha_code),
    sprintf("%-18s %s", "Raster source:", src_type),
    sprintf("%-18s %d", "Total cells:",   total_cells),
    sprintf("%-18s %d", "Valid cells:",   valid_cells),
    sprintf("%-18s %d", "Positive cells:", positive_cells),
    sprintf("%-18s %d", "Negative cells:", negative_cells),
    sprintf("%-18s %d", "Zero cells:",     zero_cells),
    sprintf("%-18s %s", "Percent positive:", fmt_pct(percent_positive)),
    sprintf("%-18s %s", "Percent negative:", fmt_pct(percent_negative)),
    sprintf("%-18s %s", "Percent zero:",     fmt_pct(percent_zero)),
    sprintf("%-18s %s", "Percent sum:",      fmt_pct(percent_sum)),
    sprintf("%-18s %s", "Extent area (km^2):", fmt_area(extent_area_km2)),
    sprintf("%-18s %s", "Outputs saved:", "1 CSV"),
    sprintf("%-18s %.3f sec", "Total elapsed:", elapsed),
    sprintf("%-18s %s", "Output file:", out_csv)
  )

  wrote_ok <- TRUE
  tryCatch({
    con <- file(out_log, open = "a", encoding = "UTF-8")
    on.exit(close(con), add = TRUE)
    writeLines(log_lines, con, sep = "\n")
    flush(con)
  }, error = function(e) {
    wrote_ok <<- FALSE
    stop("Failed to append to log file at ", out_log, " : ", conditionMessage(e),
         call. = FALSE)
  })

  if (wrote_ok) {
    size_after <- if (file.exists(out_log)) {
      suppressWarnings(file.info(out_log)$size)
    } else {
      NA_integer_
    }
    if (is.finite(size_before) && is.finite(size_after)) {
      delta <- as.integer(size_after - size_before)
      if (delta > 0) {
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

  # ---------------------------- Console summary -----------------------------
  bump("Done")
  message(sprintf(
    "Valid=%d | +%d (%.2f%%) | -%d (%.2f%%) | 0=%d (%.2f%%) | Extent %.6f km^2",
    valid_cells, positive_cells, percent_positive,
    negative_cells, percent_negative,
    zero_cells,     percent_zero,
    if (is.na(extent_area_km2)) NA_real_ else extent_area_km2
  ))

  invisible(res)
}
