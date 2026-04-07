#' Analyze temporal trends in suitability centroids
#'
#' Fits Bayesian linear trends to species' annual weighted centroids
#' (latitude and longitude) returned by [find_weighted_centroid()].
#' Exports plots, CSV summaries, and model objects to the species run
#' directory.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Pipeline context.}
#' The function:
#' \itemize{
#'   \item Loads annual centroid data
#'   \item Performs input validation and dependency checks
#'   \item Fits separate \code{lat ~ year} and \code{lon ~ year}
#'         Bayesian GLMs via rstanarm
#'   \item Summarizes slope posteriors with bayestestR (95\% CI, PD)
#'         and ROPE for the slope with a user-set practical band
#'   \item Produces prediction ribbons and lines with ggplot2
#'   \item Writes PNG, CSV, and RDS outputs
#'   \item Appends a human-readable block to the species log file
#' }
#'
#' \strong{Inputs.}
#' Upstream requirement: A function
#' \code{find_weighted_centroid(alpha_code)} must be discoverable on the
#' R search path or in one of the standard rENM locations (see file
#' discovery below). It must return a data frame with columns:
#' \itemize{
#'   \item \code{year} (numeric or integer; analysis year)
#'   \item \code{lon}  (numeric; centroid longitude)
#'   \item \code{lat}  (numeric; centroid latitude)
#' }
#' and may optionally include \code{status} where rows with
#' \code{status == "ok"} are kept.
#'
#' \strong{Directory layout.}
#' \itemize{
#'   \item Species root:
#'     \code{<project_dir>/runs/}
#'     \code{<alpha_code>/}
#'   \item Outputs:
#'     \code{<project_dir>/runs/}
#'     \code{<alpha_code>/Trends/centroids/}
#'   \item Log:
#'     \code{<project_dir>/runs/}
#'     \code{<alpha_code>/_log.txt}
#' }
#'
#' \strong{File discovery for find_weighted_centroid.R.}
#' Searched in order:
#' \itemize{
#'   \item \code{find_weighted_centroid.R} (working directory)
#'   \item \code{<project_dir>/find_weighted_centroid.R}
#'   \item \code{<project_dir>/R/find_weighted_centroid.R}
#'   \item \code{<project_dir>/scripts/find_weighted_centroid.R}
#'   \item \code{<project_dir>/code/find_weighted_centroid.R}
#' }
#'
#' \strong{Methods.}
#' \itemize{
#'   \item Separate Gaussian models are fit with
#'         \code{rstanarm::stan_glm()}:
#'         \code{lat ~ (year - mean(year))} and
#'         \code{lon ~ (year - mean(year))}
#'   \item Slope posteriors are summarized with 95\% CI,
#'         \code{p_direction} (PD), and ROPE (percent of posterior
#'         within the practical band)
#'   \item Posterior expected values (\code{posterior_epred}) are
#'         computed on a 5-year grid for plotting ribbons and trend lines
#'   \item A fixed seed (1234) is used for reproducibility
#' }
#'
#' \strong{ROPE decision rule (slope).}
#' Let \code{w = rope_width_deg_per_year}. We compute
#' \code{rope_pct = 100 * mean(beta in [-w, +w])} from the slope
#' posterior \code{beta}. We report:
#' \itemize{
#'   \item \code{inside}   if \code{rope_pct == 100}
#'   \item \code{outside}  if \code{rope_pct == 0}
#'   \item \code{overlaps} otherwise
#' }
#'
#' \strong{Outputs.}
#' The following files are created in \code{Trends/centroids/}
#' (overwritten if present):
#' \itemize{
#'   \item \code{"<ALPHA>-Centroids-Latitude-Trend.png"}
#'   \item \code{"<ALPHA>-Centroids-Longitude-Trend.png"}
#'   \item \code{"<ALPHA>-Centroids-Latitude-Summary.csv"}
#'   \item \code{"<ALPHA>-Centroids-Longitude-Summary.csv"}
#'   \item \code{"<ALPHA>-Centroids-Latitude-Model.rds"}
#'   \item \code{"<ALPHA>-Centroids-Longitude-Model.rds"}
#' }
#' A summary block is also appended to \code{"_log.txt"} under the
#' species root.
#'
#' \strong{Data requirements and failure modes.}
#' \itemize{
#'   \item Missing packages: a clear error lists what to install
#'   \item Missing or invalid find_weighted_centroid.R or malformed return
#'   \item Fewer than 3 valid years (insufficient to estimate a trend)
#' }
#'
#' @param alpha_code Character. Four-letter species alpha code
#' (e.g., \code{"CASP"}). Must be exactly 4 uppercase letters; used to
#' locate and name outputs under
#' \code{<project_dir>/runs/}
#' \code{<alpha_code>/}.
#'
#' @param rope_width_deg_per_year Numeric. Half-width of the Region of
#' Practical Equivalence (ROPE) for the slope (trend per year), expressed
#' in degrees per year. The ROPE is taken as
#' \code{[-rope_width_deg_per_year, +rope_width_deg_per_year]}.
#' Defaults to \code{0.05}, i.e., +/- 0.05 degrees per year. Set based on
#' ecological relevance (e.g., what annual shift magnitude is practically
#' negligible for your application).
#'
#' @param rope_ci Numeric. Credible interval level (0-1) used for slope
#' CI summaries. Defaults to \code{0.95} for 95\% CI (and does not affect
#' ROPE itself).
#'
#' @return
#' Invisibly returns a named list with:
#' \itemize{
#'   \item \code{data}: Data frame containing filtered analysis data with
#'         columns \code{year}, \code{lon}, \code{lat}, \code{year_c}
#'   \item \code{fit_lat}, \code{fit_lon}: rstanarm \code{stanreg} model
#'         objects for latitude and longitude trends
#'   \item \code{lat_graph}, \code{lon_graph}: ggplot objects for the
#'         fitted trends
#'   \item \code{summaries}: List with data frames \code{lat} and
#'         \code{lon} containing:
#'         \code{term}, \code{mean}, \code{median}, \code{ci_low},
#'         \code{ci_high}, \code{pd}, \code{rope_low},
#'         \code{rope_high}, \code{rope_pct}, \code{rope_decision}
#'   \item \code{outputs}: List of absolute file paths to written PNG,
#'         CSV, and RDS files
#' }
#'
#' @importFrom bayestestR ci p_direction
#' @importFrom rstanarm stan_glm posterior_epred
#' @importFrom stats gaussian quantile
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_ribbon geom_line labs
#' @importFrom ggplot2 theme_bw theme element_blank ggsave
#'
#' @examples
#' \dontrun{
#' analyze_weighted_centroids("CASP")
#' }
#'
#' @seealso
#' [find_weighted_centroid()], [rstanarm::stan_glm()],
#' [rstanarm::posterior_epred()], [bayestestR::ci()],
#' [bayestestR::p_direction()], [ggplot2::ggplot()]
#'
#' @export
analyze_weighted_centroids <- function(alpha_code,
                                       rope_width_deg_per_year = 0.05,
                                       rope_ci = 0.95) {

  # ---- Input checks ---------------------------------------------------------
  if (!is.character(alpha_code) || length(alpha_code) != 1L || nchar(alpha_code) != 4L) {
    stop("alpha_code must be a 4-letter string, e.g., 'CASP'.")
  }
  if (!is.numeric(rope_width_deg_per_year) ||
      length(rope_width_deg_per_year) != 1L ||
      rope_width_deg_per_year < 0) {
    stop("rope_width_deg_per_year must be a single non-negative numeric (degrees/year).")
  }
  if (!is.numeric(rope_ci) || length(rope_ci) != 1L || rope_ci <= 0 || rope_ci >= 1) {
    stop("rope_ci must be a single numeric in (0, 1), e.g., 0.95.")
  }

  # ---- Project directory ----------------------------------------------------
  project_dir <- rENM_project_dir()

  # ---- Dependencies ---------------------------------------------------------
  required_pkgs <- c("ggplot2", "rstanarm", "bayestestR", "dplyr", "readr")
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs)) {
    stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), ".")
  }
  `%>%` <- dplyr::`%>%`

  # ---- Source find_weighted_centroid.R --------------------------------------
  candidate_paths <- c(
    "find_weighted_centroid.R",
    file.path(project_dir, "find_weighted_centroid.R"),
    file.path(project_dir, "R", "find_weighted_centroid.R"),
    file.path(project_dir, "scripts", "find_weighted_centroid.R"),
    file.path(project_dir, "code", "find_weighted_centroid.R")
  )
  src_ok <- FALSE
  for (p in candidate_paths) {
    if (file.exists(p)) {
      tryCatch({
        sys.source(p, envir = environment())
        src_ok <- TRUE
        break
      }, error = function(e) NULL)
    }
  }
  if (!src_ok && !exists("find_weighted_centroid", mode = "function")) {
    stop("Cannot find find_weighted_centroid.R. Place it in one of the searched locations.")
  }

  # ---- Paths ----------------------------------------------------------------
  species_root <- file.path(project_dir, "runs", alpha_code)
  out_dir      <- file.path(species_root, "Trends", "centroids")
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  log_path <- file.path(species_root, "_log.txt")

  # ---- Filenames ------------------------------------------------------------
  f_lat_png <- file.path(out_dir, sprintf("%s-Centroids-Latitude-Trend.png",  alpha_code))
  f_lon_png <- file.path(out_dir, sprintf("%s-Centroids-Longitude-Trend.png", alpha_code))
  f_lat_csv <- file.path(out_dir, sprintf("%s-Centroids-Latitude-Summary.csv",  alpha_code))
  f_lon_csv <- file.path(out_dir, sprintf("%s-Centroids-Longitude-Summary.csv", alpha_code))
  f_lat_rds <- file.path(out_dir, sprintf("%s-Centroids-Latitude-Model.rds",    alpha_code))
  f_lon_rds <- file.path(out_dir, sprintf("%s-Centroids-Longitude-Model.rds",   alpha_code))

  # Simple HH:MM:SS.ss formatter for log readability.
  .fmt_hms <- function(sec) {
    h <- floor(sec / 3600)
    m <- floor((sec %% 3600) / 60)
    s <- sec - h * 3600 - m * 60
    sprintf("%02d:%02d:%05.2f", h, m, s)
  }
  t_start <- proc.time()[["elapsed"]]

  # ---- Data -----------------------------------------------------------------
  p_all <- find_weighted_centroid(alpha_code)
  if (!all(c("year", "lon", "lat") %in% names(p_all))) {
    stop("find_weighted_centroid() must return columns: year, lon, lat.")
  }

  p <- p_all %>%
    dplyr::filter(
      is.finite(.data$year),
      is.finite(.data$lon),
      is.finite(.data$lat)
    )
  if ("status" %in% names(p)) {
    p <- p %>% dplyr::filter(.data$status == "ok")
  }
  if (nrow(p) < 3) {
    stop("Need at least 3 valid years to fit trends.")
  }

  # ---- Modeling -------------------------------------------------------------
  yr_mean   <- mean(p$year)
  p$year_c  <- p$year - yr_mean

  set.seed(1234)
  fit_lat <- rstanarm::stan_glm(
    lat ~ year_c,
    data   = p,
    family = stats::gaussian(),
    refresh = 0
  )
  fit_lon <- rstanarm::stan_glm(
    lon ~ year_c,
    data   = p,
    family = stats::gaussian(),
    refresh = 0
  )

  # Helper: summarize slope posterior with CI, PD, and ROPE -------------------
  .slope_summary <- function(fit, term = "year_c",
                             rope_w = 0.05,
                             ci_level = 0.95) {
    draws <- as.data.frame(as.matrix(fit))
    if (!term %in% names(draws)) {
      stop("Expected term '", term, "' in posterior draws but did not find it.")
    }
    s <- draws[[term]]

    ci <- bayestestR::ci(s, ci = ci_level)
    pd <- bayestestR::p_direction(s)[["pd"]]

    rope_low  <- -abs(rope_w)
    rope_high <-  abs(rope_w)
    rope_pct  <- mean(s >= rope_low & s <= rope_high) * 100

    rope_decision <- if (isTRUE(all(s >= rope_low & s <= rope_high))) {
      "inside"
    } else if (isTRUE(all(s < rope_low | s > rope_high))) {
      "outside"
    } else {
      "overlaps"
    }

    data.frame(
      term          = "year (per 1 yr)",
      mean          = mean(s),
      median        = median(s),
      ci_low        = ci$CI_low,
      ci_high       = ci$CI_high,
      pd            = pd,
      rope_low      = rope_low,
      rope_high     = rope_high,
      rope_pct      = rope_pct,
      rope_decision = rope_decision,
      stringsAsFactors = FALSE
    )
  }

  sum_lat <- .slope_summary(
    fit_lat,
    rope_w   = rope_width_deg_per_year,
    ci_level = rope_ci
  )
  sum_lon <- .slope_summary(
    fit_lon,
    rope_w   = rope_width_deg_per_year,
    ci_level = rope_ci
  )

  # ---- Predictions ----------------------------------------------------------
  make_pred_df <- function(fit, years_seq) {
    grid <- data.frame(year = years_seq)
    grid$year_c <- grid$year - yr_mean
    ep <- rstanarm::posterior_epred(fit, newdata = grid)
    data.frame(
      year    = grid$year,
      mean    = apply(ep, 2, mean),
      lower95 = apply(ep, 2, function(z) stats::quantile(z, 0.025)),
      upper95 = apply(ep, 2, function(z) stats::quantile(z, 0.975))
    )
  }

  year_seq <- seq(min(p$year), max(p$year), by = 5)
  pred_lat <- make_pred_df(fit_lat, year_seq)
  pred_lon <- make_pred_df(fit_lon, year_seq)

  # ---- Plotting -------------------------------------------------------------
  base_theme <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  lat_graph <- ggplot2::ggplot(
    p,
    ggplot2::aes(x = .data$year, y = .data$lat)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      data = pred_lat,
      mapping = ggplot2::aes(
        x    = .data$year,
        ymin = .data$lower95,
        ymax = .data$upper95
      ),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = pred_lat,
      mapping = ggplot2::aes(
        x = .data$year,
        y = .data$mean
      ),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = paste0(toupper(alpha_code), " Centroid Latitude Trend"),
      x     = "Year",
      y     = "Latitude"
    ) +
    base_theme

  lon_graph <- ggplot2::ggplot(
    p,
    ggplot2::aes(x = .data$year, y = .data$lon)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      data = pred_lon,
      mapping = ggplot2::aes(
        x    = .data$year,
        ymin = .data$lower95,
        ymax = .data$upper95
      ),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = pred_lon,
      mapping = ggplot2::aes(
        x = .data$year,
        y = .data$mean
      ),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = paste0(toupper(alpha_code), " Centroid Longitude Trend"),
      x     = "Year",
      y     = "Longitude"
    ) +
    base_theme

  # ---- Exports --------------------------------------------------------------
  ggplot2::ggsave(f_lat_png, lat_graph, width = 7.5, height = 5.0, dpi = 300)
  ggplot2::ggsave(f_lon_png, lon_graph, width = 7.5, height = 5.0, dpi = 300)
  readr::write_csv(sum_lat, f_lat_csv)
  readr::write_csv(sum_lon, f_lon_csv)
  saveRDS(fit_lat, f_lat_rds)
  saveRDS(fit_lon, f_lon_rds)

  # ---- Log block ------------------------------------------------------------
  elapsed     <- proc.time()[["elapsed"]] - t_start
  elapsed_str <- .fmt_hms(elapsed)
  ts          <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  cat("\n", file = log_path, append = TRUE)
  header_line <- "------------------------------------------------------------------------"

  fmt_slope <- function(s) {
    sprintf(
      paste0(
        "    slope: mean=%.6f  95%%CI=[%.6f, %.6f]  PD=%.1f%%  ",
        "ROPE=[%.3f, %.3f] deg/yr  inside=%.3f%%  decision=%s"
      ),
      s$mean, s$ci_low, s$ci_high, 100 * s$pd,
      s$rope_low, s$rope_high, s$rope_pct, s$rope_decision
    )
  }

  block <- c(
    header_line,
    "Processing summary (analyze_weighted_centroids)",
    sprintf("Timestamp:          %s", ts),
    sprintf("Alpha code:         %s", alpha_code),
    sprintf("Years processed:    %d", length(unique(p$year))),
    sprintf("Completed:          %d", nrow(p)),
    "Failed:             0",
    sprintf("Total elapsed:      %s (%.2f s)", elapsed_str, elapsed),
    sprintf(
      "Outputs saved:      %s-Centroids-Latitude-Trend.png+%s-Centroids-Longitude-Trend.png+%s-Centroids-Latitude-Summary.csv+%s-Centroids-Longitude-Summary.csv",
      alpha_code, alpha_code, alpha_code, alpha_code
    ),
    sprintf("ROPE width (slope): +/- %.3f deg/yr", rope_width_deg_per_year),
    "Model summaries:",
    paste0("  Latitude trend\n",  fmt_slope(sum_lat)),
    paste0("  Longitude trend\n", fmt_slope(sum_lon)),
    "Per-year status:"
  )

  per_year <- vapply(
    sort(unique(p$year)),
    function(y) {
      row <- p[p$year == y, , drop = FALSE][1, ]
      sprintf(
        "  %4d  status=plotted lon=%.6f lat=%.6f notes=-",
        y, row$lon, row$lat
      )
    },
    character(1)
  )

  cat(
    paste0(paste(c(block, per_year), collapse = "\n"), "\n"),
    file   = log_path,
    append = TRUE
  )

  # ---- Return ---------------------------------------------------------------
  invisible(list(
    data      = p,
    fit_lat   = fit_lat,
    fit_lon   = fit_lon,
    lat_graph = lat_graph,
    lon_graph = lon_graph,
    summaries = list(lat = sum_lat, lon = sum_lon),
    outputs   = list(
      lat_png = f_lat_png,
      lon_png = f_lon_png,
      lat_csv = f_lat_csv,
      lon_csv = f_lon_csv,
      lat_rds = f_lat_rds,
      lon_rds = f_lon_rds
    )
  ))
}
