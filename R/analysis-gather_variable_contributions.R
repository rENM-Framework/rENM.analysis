#' Gather variable contributions across all modeled years
#'
#' Consolidate ranked variable importance outputs across time-series runs
#' into a single summary dataset for a species. Produces a unified CSV and
#' appends a detailed processing log with per-year results.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Inputs}
#' \itemize{
#'   \item Ranked variable importance files (PI-Ranked):
#'   \code{<project_dir>/runs/<alpha_code>/TimeSeries/<year>/model/}
#'   \code{<alpha_code>-<year>-PI-Ranked.txt}
#'   \item Years processed: 1980, 1985, ..., 2020
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Consolidated CSV:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/variables/}
#'   \code{<alpha_code>-Variables-AllYears.csv}
#'   \item Processing log appended to:
#'   \code{<project_dir>/runs/<alpha_code>/}
#'   \code{_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Reads per-span ranked variable importance files across years
#'   \item Tolerates tab, comma, and whitespace-delimited formats
#'   \item Uses heuristics to identify Variable and Percent columns
#'   \item Cleans and converts percentage values
#'   \item Aggregates results into long and wide formats
#'   \item Computes summary statistics:
#'   \itemize{
#'     \item n_years
#'     \item mean_pct
#'     \item median_pct
#'     \item sd_pct
#'     \item min_pct
#'     \item max_pct
#'   }
#'   \item Writes a consolidated CSV and logs per-year processing status
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Input files must exist for one or more modeled years
#'   \item Files must contain variable and percent contribution columns
#'   \item Missing or malformed files are tolerated and logged
#' }
#'
#' @param alpha_code Character. Four-letter banding code (e.g., "CASP"),
#' coerced to upper case.
#'
#' @return Invisibly returns a List with:
#' \itemize{
#'   \item long: Data frame of all variable contributions by year
#'   \item wide_summary: Data frame of summary statistics and per-year columns
#'   \item saved_csv: Character path to the output CSV
#'   \item log_path: Character path to the appended log file
#' }
#'
#' Side effects:
#' \itemize{
#'   \item Writes a consolidated CSV file to disk
#'   \item Appends a detailed processing summary (with per-year lines) to
#'   \code{_log.txt}
#' }
#'
#' @importFrom utils read.table write.csv
#' @importFrom stats aggregate median sd
#'
#' @examples
#' \dontrun{
#'   gather_variable_contributions("CASP")
#' }
#'
#' @export
gather_variable_contributions <- function(alpha_code) {
  ## -------- Input validation -------------------------------------------------
  if (missing(alpha_code) || length(alpha_code) != 1 || !nzchar(alpha_code)) {
    stop("`alpha_code` must be a non-empty character scalar (e.g., 'CASP').")
  }
  alpha_code <- toupper(trimws(alpha_code))
  if (!grepl("^[A-Z]{4}$", alpha_code)) {
    stop("`alpha_code` must be exactly four letters A-Z (e.g., 'CASP').")
  }

  ## -------- Paths ------------------------------------------------------------
  project_dir <- rENM_project_dir()

  years <- seq(1980, 2020, by = 5)
  base_species_dir <- file.path(project_dir, "runs", alpha_code)
  time_series_dir  <- file.path(base_species_dir, "TimeSeries")
  out_dir          <- file.path(base_species_dir, "Trends", "variables")
  out_csv          <- file.path(out_dir, sprintf("%s-Variables-AllYears.csv", alpha_code))
  log_path         <- file.path(base_species_dir, "_log.txt")

  if (!dir.exists(base_species_dir)) dir.create(base_species_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(out_dir))         dir.create(out_dir,         recursive = TRUE, showWarnings = FALSE)

  ## -------- Helpers ----------------------------------------------------------
  .msg <- function(...) cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] "),
                            sprintf(...), "\n")

  .read_ranked <- function(path) {
    try_read <- function(sep, header_guess = TRUE) {
      utils::read.table(path, header = header_guess, sep = sep,
                        stringsAsFactors = FALSE, quote = "", comment.char = "",
                        fill = TRUE, check.names = FALSE)
    }
    for (sep in c("\t", ",", "")) {
      for (hdr in c(TRUE, FALSE)) {
        res <- try(suppressWarnings(try_read(sep, hdr)), silent = TRUE)
        if (!inherits(res, "try-error") && is.data.frame(res) && nrow(res) > 0) return(res)
      }
    }
    stop("Unable to read table: ", path)
  }

  .find_cols <- function(df) {
    nms <- tolower(gsub("\\s+", "", names(df)))
    var_candidates <- grepl("^(variable|var|varname|predictor|feature)$", nms) |
      grepl("variable|varname|predictor|feature", nms)
    pct_candidates <- grepl("^(percent|perc|percentage|%|contribution|percentcontribution|%contribution)$", nms) |
      grepl("percent|perc|percentage|contribution", nms) | grepl("%", nms, fixed = TRUE)
    var_col <- which(var_candidates); if (!length(var_col)) var_col <- 1
    pct_col <- which(pct_candidates); if (!length(pct_col)) pct_col <- if (length(nms) >= 2) 2 else 1
    list(var = var_col[1], pct = pct_col[1])
  }

  .as_pct <- function(x) {
    x <- trimws(as.character(x))
    x <- gsub("%", "", x, fixed = TRUE)
    x <- gsub(",", "", x, fixed = TRUE)
    suppressWarnings(as.numeric(x))
  }

  .fmt_hms <- function(seconds_num) {
    s <- max(0, as.numeric(seconds_num)); h <- s %/% 3600; m <- (s %% 3600) %/% 60; sec <- floor(s %% 60)
    sprintf("%02d:%02d:%02d", as.integer(h), as.integer(m), as.integer(sec))
  }

  ## -------- Log writer -------------------------------------------------------
  .write_log_gvc <- function(alpha_code, years, per_year, workers,
                             out_csv, total_elapsed, log_path) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
    years_processed <- length(years)
    completed_count <- sum(vapply(per_year, function(x) identical(x$status, "completed"), logical(1)))
    failed_count    <- years_processed - completed_count

    have_years <- integer(0)
    if (length(per_year)) {
      have_years <- vapply(per_year, function(x) as.integer(x$year), integer(1))
      per_year <- per_year[order(have_years)]
    }

    missing_yrs <- setdiff(as.integer(years),
                           if (length(per_year)) vapply(per_year, function(x) as.integer(x$year), integer(1)) else integer(0))
    if (length(missing_yrs)) {
      for (yr in missing_yrs) {
        per_year[[length(per_year) + 1L]] <- list(
          year = as.integer(yr), status = "missing", elapsed = 0, outputs = "-", notes = "not found"
        )
      }
      per_year <- per_year[order(vapply(per_year, function(x) as.integer(x$year), integer(1)))]
    }

    per_year_lines <- vapply(per_year, function(x) {
      sprintf("  %4d  status=%-9s elapsed=%s   outputs=%-12s   notes=%s",
              as.integer(x$year),
              x$status,
              .fmt_hms(x$elapsed),
              if (is.null(x$outputs) || !nzchar(x$outputs)) "pi-ranked" else x$outputs,
              if (is.null(x$notes) || is.na(x$notes) || !nzchar(x$notes)) "-" else x$notes)
    }, character(1))

    header <- c(
      "\n", strrep("-", 72),
      "Processing summary (gather_variable_contributions)",
      sprintf("%-20s %s", "Timestamp:", ts),
      sprintf("%-20s %s", "Alpha code:", alpha_code),
      sprintf("%-20s %d", "Years processed:", years_processed),
      sprintf("%-20s %d", "Parallel workers:", as.integer(workers)),
      sprintf("%-20s %d", "Completed:", completed_count),
      sprintf("%-20s %d", "Failed:", failed_count),
      sprintf("%-20s %s (%.2f s)", "Total elapsed:", .fmt_hms(total_elapsed), total_elapsed),
      "Per-year status:"
    )
    footer <- c(
      "Outputs:",
      paste0("  ", out_csv),
      ""
    )

    block <- c(header, per_year_lines, footer)

    if (!dir.exists(dirname(log_path))) dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
    con <- file(log_path, open = "a", encoding = "native.enc")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    writeLines(block, con = con, sep = "\n", useBytes = TRUE)
    flush(con)
  }

  ## -------- Processing -------------------------------------------------------
  .msg("Starting gather_variable_contributions for '%s'.", alpha_code)
  start_time <- Sys.time()

  per_year <- list()
  rows <- list()

  for (yr in years) {
    yr_t0 <- Sys.time()
    in_file <- file.path(time_series_dir, as.character(yr), "model",
                         sprintf("%s-%s-PI-Ranked.txt", alpha_code, yr))

    if (!file.exists(in_file)) {
      .msg("... %d: file not found -> %s", yr, in_file)
      per_year[[length(per_year) + 1L]] <- list(
        year = as.integer(yr), status = "missing",
        elapsed = as.numeric(difftime(Sys.time(), yr_t0, units = "secs")),
        outputs = "-", notes = "file not found"
      )
      next
    }

    .msg("... %d: reading %s", yr, in_file)
    df <- try(.read_ranked(in_file), silent = TRUE)
    if (inherits(df, "try-error")) {
      .msg("WARNING: Failed to parse %s (skipping).", in_file)
      per_year[[length(per_year) + 1L]] <- list(
        year = as.integer(yr), status = "failed",
        elapsed = as.numeric(difftime(Sys.time(), yr_t0, units = "secs")),
        outputs = "pi-ranked", notes = "parse error"
      )
      next
    }

    cols <- .find_cols(df)
    v <- df[[cols$var]]
    p <- .as_pct(df[[cols$pct]])
    keep <- !(is.na(v) | v == "" | is.na(p))
    v <- v[keep]; p <- p[keep]

    elapsed_sec <- as.numeric(difftime(Sys.time(), yr_t0, units = "secs"))
    if (!length(v)) {
      .msg("WARNING: No valid rows after cleaning for %d (skipping).", yr)
      per_year[[length(per_year) + 1L]] <- list(
        year = as.integer(yr), status = "no-data",
        elapsed = elapsed_sec, outputs = "pi-ranked", notes = "no valid rows"
      )
      next
    }

    rows[[length(rows) + 1L]] <- data.frame(
      Variable = as.character(v),
      Year     = as.integer(yr),
      Percent  = as.numeric(p),
      stringsAsFactors = FALSE
    )

    .msg("...... %d: kept %d rows.", yr, length(v))

    per_year[[length(per_year) + 1L]] <- list(
      year = as.integer(yr), status = "completed",
      elapsed = elapsed_sec, outputs = "pi-ranked", notes = paste0("kept=", length(v))
    )
  }

  end_time <- Sys.time()
  total_elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (!length(rows)) {
    .write_log_gvc(alpha_code, years, per_year, workers = 1L,
                   out_csv = "(no CSV written)", total_elapsed = total_elapsed, log_path = log_path)
    stop("No input files parsed successfully for alpha_code='", alpha_code, "'.")
  }

  long <- do.call(rbind, rows)

  ## -------- Summaries -------------------------------------------------------
  .msg("Computing summaries and shaping wide table...")

  agg_n      <- aggregate(Percent ~ Variable, data = long, FUN = function(x) sum(!is.na(x)))
  names(agg_n)[2] <- "n_years"

  agg_mean   <- aggregate(Percent ~ Variable, data = long, FUN = function(x) mean(x, na.rm = TRUE))
  names(agg_mean)[2] <- "mean_pct"

  agg_median <- aggregate(Percent ~ Variable, data = long, FUN = function(x) stats::median(x, na.rm = TRUE))
  names(agg_median)[2] <- "median_pct"

  agg_sd     <- aggregate(Percent ~ Variable, data = long, FUN = function(x) stats::sd(x, na.rm = TRUE))
  names(agg_sd)[2] <- "sd_pct"

  agg_min    <- aggregate(Percent ~ Variable, data = long, FUN = function(x) min(x, na.rm = TRUE))
  names(agg_min)[2] <- "min_pct"

  agg_max    <- aggregate(Percent ~ Variable, data = long, FUN = function(x) max(x, na.rm = TRUE))
  names(agg_max)[2] <- "max_pct"

  summary_df <- Reduce(function(a, b) merge(a, b, by = "Variable", all = TRUE),
                       list(agg_n, agg_mean, agg_median, agg_sd, agg_min, agg_max))

  all_vars <- sort(unique(long$Variable))
  wide <- data.frame(Variable = all_vars, stringsAsFactors = FALSE)

  for (yr in years) {
    col_name <- paste0("Y", yr)
    sub <- long[long$Year == yr, c("Variable", "Percent")]
    m <- match(wide$Variable, sub$Variable)
    wide[[col_name]] <- ifelse(is.na(m), NA_real_, sub$Percent[m])
  }

  wide_summary <- merge(summary_df, wide, by = "Variable", all.x = TRUE, sort = TRUE)
  wide_summary <- wide_summary[order(-wide_summary$mean_pct, wide_summary$Variable), ]

  ## -------- Write CSV --------------------------------------------------------
  wrote_csv <- TRUE
  tryCatch(
    utils::write.csv(wide_summary, file = out_csv, row.names = FALSE),
    error = function(e) {
      wrote_csv <<- FALSE
      .msg("ERROR writing CSV (%s): %s", out_csv, conditionMessage(e))
    }
  )
  if (wrote_csv) .msg("Saved summary CSV: %s", out_csv)

  ## -------- Log (always) -----------------------------------------------------
  .write_log_gvc(alpha_code, years, per_year, workers = 1L,
                 out_csv = if (wrote_csv) out_csv else "(CSV write failed)",
                 total_elapsed = total_elapsed, log_path = log_path)

  invisible(list(
    long = long,
    wide_summary = wide_summary,
    saved_csv = if (wrote_csv) out_csv else NA_character_,
    log_path = log_path
  ))
}
