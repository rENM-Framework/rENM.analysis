#' Analyze temporal trends in suitability centroids
#'
#' Fits Bayesian linear trends to a species' weighted centroids (latitude and
#' longitude) returned by \code{\link{find_weighted_centroid}} and exports
#' plots, CSV summaries, and model objects to the species run directory.
#'
#' @details
#' Fits separate Gaussian models via \code{rstanarm::stan_glm()}:
#' \code{lat ~ (year - mean(year))} and \code{lon ~ (year - mean(year))}.
#' Slope posteriors are summarized with 95\% CI, \code{p_direction} (PD),
#' and ROPE (percent of posterior within the practical equivalence band).
#' Posterior expected values are computed on a 5-year grid for plotting
#' ribbons and trend lines. A fixed seed (1234) is used for reproducibility.
#'
#' \strong{ROPE decision rule (slope).}
#' Let \code{w = rope_width_deg_per_year}. The posterior is evaluated against
#' \code{[-w, +w]}. The decision is \code{"inside"} if 100\% of draws fall
#' within ROPE, \code{"outside"} if 0\% do, and \code{"overlaps"} otherwise.
#'
#' Outputs written to \code{<project_dir>/runs/<alpha_code>/Trends/centroids/}:
#' \itemize{
#'   \item \code{<ALPHA>-Centroids-Latitude-Trend.png}
#'   \item \code{<ALPHA>-Centroids-Longitude-Trend.png}
#'   \item \code{<ALPHA>-Centroids-Latitude-Summary.csv}
#'   \item \code{<ALPHA>-Centroids-Longitude-Summary.csv}
#'   \item \code{<ALPHA>-Centroids-Latitude-Model.rds}
#'   \item \code{<ALPHA>-Centroids-Longitude-Model.rds}
#' }
#' A summary block is appended to
#' \code{<project_dir>/runs/<alpha_code>/_log.txt}.
#'
#' @param alpha_code Character. Four-letter species alpha code (e.g.,
#'   \code{"CASP"}). Must be exactly 4 uppercase letters.
#' @param rope_width_deg_per_year Numeric. Half-width of the ROPE for the
#'   slope (degrees per year). Default 0.05.
#' @param rope_ci Numeric. Credible interval level (0--1) for slope CI
#'   summaries. Default 0.95.
#'
#' @return Invisible list with:
#'   \code{data} (filtered data frame),
#'   \code{fit_lat}, \code{fit_lon} (stanreg objects),
#'   \code{lat_graph}, \code{lon_graph} (ggplot objects),
#'   \code{summaries} (list of per-dimension summary data frames),
#'   \code{outputs} (named list of file paths).
#'   Primary side effects are six files written to disk and a log entry.
#'
#' @seealso \code{\link{find_weighted_centroid}}, \code{\link[rENM.core]{rENM_project_dir}}
#' @importFrom bayestestR ci p_direction
#' @importFrom rstanarm stan_glm posterior_epred
#' @importFrom stats gaussian quantile
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_point geom_ribbon geom_line labs
#' @importFrom ggplot2 theme_bw theme element_blank ggsave
#'
#' @examples
#' \dontrun{
#'   analyze_weighted_centroids("CASP")
#' }
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

  # HH:MM:SS.ss formatter
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
    data    = p,
    family  = stats::gaussian(),
    refresh = 0
  )
  fit_lon <- rstanarm::stan_glm(
    lon ~ year_c,
    data    = p,
    family  = stats::gaussian(),
    refresh = 0
  )

  # Summarize slope posterior with CI, PD, and ROPE
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

  sum_lat <- .slope_summary(fit_lat, rope_w = rope_width_deg_per_year, ci_level = rope_ci)
  sum_lon <- .slope_summary(fit_lon, rope_w = rope_width_deg_per_year, ci_level = rope_ci)

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
      mapping = ggplot2::aes(x = .data$year, ymin = .data$lower95, ymax = .data$upper95),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = pred_lat,
      mapping = ggplot2::aes(x = .data$year, y = .data$mean),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = paste0(toupper(alpha_code), " Centroid Latitude Trend"),
      x = "Year", y = "Latitude"
    ) +
    base_theme

  lon_graph <- ggplot2::ggplot(
    p,
    ggplot2::aes(x = .data$year, y = .data$lon)
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_ribbon(
      data = pred_lon,
      mapping = ggplot2::aes(x = .data$year, ymin = .data$lower95, ymax = .data$upper95),
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = pred_lon,
      mapping = ggplot2::aes(x = .data$year, y = .data$mean),
      linewidth = 1,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(
      title = paste0(toupper(alpha_code), " Centroid Longitude Trend"),
      x = "Year", y = "Longitude"
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
