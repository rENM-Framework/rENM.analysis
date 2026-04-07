#' Summarize and plot variable contributions
#'
#' Reads a species-specific wide-format CSV of variable contributions
#' across years, selects variables by mean contribution, and produces
#' a complete suite of visualizations and statistics.
#'
#' @details
#' This function is part of the rENM framework's processing pipeline
#' and operates within the project directory structure defined by
#' rENM_project_dir().
#'
#' \strong{Input}
#' \itemize{
#'   \item File:
#'   \code{<project_dir>/runs/<alpha_code>/Trends/variables/}
#'   \code{<alpha_code>-Variables-AllYears.csv}
#' }
#'
#' \strong{Processing rules}
#' \itemize{
#'   \item Eligibility rule (2025-09): Only variables that appear
#'   (non-NA contribution) in 3 or more of the years 1980-2020 are
#'   considered for \code{top_n} ranking and plotting.
#'   \item Adds Region of Practical Equivalence (ROPE) calculations
#'   for the Bayesian slope (and optional intercept).
#' }
#'
#' \strong{Outputs}
#' \itemize{
#'   \item Multiple PNG visualizations:
#'   \itemize{
#'     \item Heatmap
#'     \item Line plots
#'     \item Frequentist linear regression (lines and points)
#'     \item Bayesian regression (lines and points, with and without ribbon)
#'   }
#'   \item CSV summaries:
#'   \itemize{
#'     \item Linear regression statistics
#'     \item Bayesian regression statistics (including ROPE metrics)
#'   }
#'   \item Processing summary appended to:
#'   \code{<project_dir>/runs/<alpha_code>/}
#'   \code{_log.txt}
#' }
#'
#' \strong{Methods}
#' \itemize{
#'   \item Variable selection based on mean contribution and eligibility rule
#'   \item Frequentist linear regression (LR) across years
#'   \item Bayesian regression (BR) with posterior summaries
#'   \item ROPE calculations:
#'   \itemize{
#'     \item Probability mass within ROPE
#'     \item Decision classification: inside, outside, overlaps
#'   }
#' }
#'
#' \strong{Data requirements}
#' \itemize{
#'   \item Input CSV must exist and contain per-year contribution columns
#'   \item Years must span 1980-2020 in 5-year increments
#'   \item Variables must have sufficient non-NA observations for inclusion
#' }
#'
#' @param alpha_code Character. Four-letter banding code (e.g., "CASP").
#' @param top_n Integer. Number of top variables to include in plots and
#' summaries; must be a positive integer or NULL. (Default = 10)
#' @param rope_slope Numeric. Two-element vector defining ROPE bounds for
#' Bayesian slope (default c(-0.05, 0.05)); a single value is expanded
#' symmetrically.
#' @param rope_intercept Numeric. Two-element vector defining ROPE bounds
#' for Bayesian intercept, or NULL (default) to disable.
#'
#' @return Invisibly returns NULL. Side effects:
#' \itemize{
#'   \item Writes multiple PNG plots to disk
#'   \item Writes CSV files containing LR and BR statistics
#'   \item Appends a processing summary to \code{_log.txt}
#' }
#'
#' @examples
#' \dontrun{
#'   summarize_variable_contributions("CASP")
#' }
#'
#' @export
summarize_variable_contributions <- function(alpha_code,
                                             top_n = 10,
                                             rope_slope = c(-0.05, 0.05),
                                             rope_intercept = NULL) {
  start_time <- Sys.time()

  # number of years to include in trend analyses
  included_years <- 3

  # ---- [0] helpers (ROPE) ---------------------------------------------------
  .norm_rope <- function(x) {
    if (is.null(x)) return(NULL)
    if (length(x) == 1L) x <- c(-abs(x), abs(x))
    stopifnot(is.numeric(x), length(x) == 2L)
    x <- sort(as.numeric(x))
    if (!is.finite(x[1]) || !is.finite(x[2]) || x[1] >= x[2]) stop("Invalid ROPE bounds")
    x
  }
  rope_slope <- .norm_rope(rope_slope)
  rope_intercept <- .norm_rope(rope_intercept)

  rope_prob <- function(draws, rope) mean(draws >= rope[1] & draws <= rope[2]) * 100
  rope_decision <- function(ci_lo, ci_hi, rope) {
    if (is.null(rope) || any(!is.finite(c(ci_lo, ci_hi)))) return(NA_character_)
    if (ci_hi < rope[1] || ci_lo > rope[2]) return("outside")
    if (ci_lo >= rope[1] && ci_hi <= rope[2]) return("inside")
    "overlaps"
  }

  # ---- [1/13] validate inputs ----------------------------------------------
  message(">>> [1/13] Validating inputs...")
  if (missing(alpha_code) || !is.character(alpha_code) || length(alpha_code) != 1) {
    stop("`alpha_code` must be a single character value (e.g., 'CASP').")
  }
  if (!is.null(top_n)) {
    if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) || top_n <= 0) {
      stop("`top_n` must be a single positive integer or NULL.")
    }
    top_n <- as.integer(top_n)
  }

  # ---- [2/13] paths ----------------------------------------------------------
  message(">>> [2/13] Resolving file paths...")
  project_dir <- rENM_project_dir()

  base_dir <- file.path(project_dir, "runs", alpha_code)
  vars_dir <- file.path(base_dir, "Trends", "variables")
  csv_path <- file.path(vars_dir, paste0(alpha_code, "-Variables-AllYears.csv"))
  log_path <- file.path(base_dir, "_log.txt")

  # PNGs
  png_heatmap       <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-Heatmap.png"))
  png_lines         <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-Lines.png"))
  png_lr            <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-LR.png"))
  png_lr_pts        <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-LR-Points.png"))
  png_br            <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-BR-Lines.png"))
  png_br_pts        <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-BR-Points.png"))
  png_br_norib      <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-BR-Lines-NoRibbon.png"))
  png_br_pts_norib  <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-BR-Points-NoRibbon.png"))

  # CSVs
  stats_lr_path <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-LR-Stats.csv"))
  stats_br_path <- file.path(vars_dir, paste0(alpha_code, "-Variable-Contributions-BR-Stats.csv"))

  if (!file.exists(csv_path)) stop(sprintf("Input not found: %s", csv_path))

  # ---- [3/13] packages -------------------------------------------------------
  message(">>> [3/13] Loading required packages...")
  req_pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "tidyselect")
  missing_pkgs <- req_pkgs[!vapply(req_pkgs, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(sprintf("Missing required packages: %s. Please install them.",
                 paste(missing_pkgs, collapse = ", ")))
  }
  if (!requireNamespace("rstanarm", quietly = TRUE)) {
    stop("Package 'rstanarm' is required for Bayesian regression. Please install it: install.packages('rstanarm')")
  }

  readr    <- asNamespace("readr")
  dplyr    <- asNamespace("dplyr")
  tidyr    <- asNamespace("tidyr")
  ggplot2  <- asNamespace("ggplot2")
  stringr  <- asNamespace("stringr")
  rstanarm <- asNamespace("rstanarm")

  # Optional: markdown-capable legend labels if ggtext is installed
  has_ggtext <- "ggtext" %in% rownames(utils::installed.packages())
  if (has_ggtext) ggtext <- asNamespace("ggtext")

  # ---- Reusable themes -------------------------------------------------------
  trend_theme_base <- function() {
    base <- ggplot2::theme_minimal(base_size = 16) +
      ggplot2::theme(
        plot.title          = ggplot2::element_text(face = "bold", size = 32, hjust = 0.5),
        plot.title.position = "plot",
        axis.title          = ggplot2::element_text(size = 18),
        axis.text           = ggplot2::element_text(size = 16),
        legend.title        = ggplot2::element_text(size = 18, face = "bold"),
        strip.text          = ggplot2::element_text(size = 16),
        legend.key.size     = grid::unit(22, "pt"),
        legend.key.width    = grid::unit(28, "pt"),
        legend.key.height   = grid::unit(18, "pt"),
        panel.background    = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background     = ggplot2::element_rect(fill = "white", colour = NA),
        panel.border        = ggplot2::element_rect(color = "darkgray", fill = NA, linewidth = 0.8)
      )
    if (has_ggtext) {
      base <- base + ggplot2::theme(legend.text = ggtext$element_markdown(size = 16))
    } else {
      base <- base + ggplot2::theme(legend.text = ggplot2::element_text(size = 16))
    }
    base
  }
  trend_theme_br <- function() {
    trend_theme_base()
  }

  # ---- [4/13] read -----------------------------------------------------------
  message(">>> [4/13] Reading input CSV: ", csv_path)
  df <- readr$read_csv(csv_path, show_col_types = FALSE)

  req_cols <- c("Variable", "mean_pct")
  if (!all(req_cols %in% names(df))) {
    stop("Input file must include columns: ", paste(req_cols, collapse = ", "))
  }

  # ---- [5/13] detect years ---------------------------------------------------
  year_cols <- names(df)[stringr::str_detect(names(df), "^Y\\d{4}$")]
  if (length(year_cols) == 0) {
    stop("No year columns found. Expected columns named like Y1980, Y1985, ..., Y2020.")
  }
  year_nums <- as.integer(stringr::str_replace(year_cols, "^Y", ""))
  ord <- order(year_nums); year_cols <- year_cols[ord]; year_nums <- year_nums[ord]
  message("    Detected years: ", paste(year_nums, collapse = ", "))

  # ---- [6/13] eligibility (>=3 years) + select top_n by mean_pct ------------
  message(">>> [6/13] Applying eligibility rule (present in >= 3 years) and ranking by mean_pct...")
  df$years_present <- rowSums(!is.na(df[, year_cols, drop = FALSE]))
  df_eligible <- df[df$years_present >= included_years, , drop = FALSE]
  if (nrow(df_eligible) == 0L) stop("No variables meet the eligibility rule (present in >= 3 years).")

  df_ranked <- df_eligible |>
    dplyr::arrange(dplyr::desc(mean_pct)) |>
    dplyr::mutate(rank = dplyr::row_number())

  n_available <- nrow(df_ranked)
  n_requested <- if (is.null(top_n)) n_available else top_n
  n_use <- min(n_available, n_requested)

  top_vars <- df_ranked |>
    dplyr::slice_head(n = n_use) |>
    dplyr::select(Variable, mean_pct, rank, years_present)

  message("    Eligible variables: ", n_available,
          " | Selected (top ", n_use, "): ",
          paste(top_vars$Variable, collapse = ", "))

  # ---- [7/13] reshape to long ------------------------------------------------
  message(">>> [7/13] Reshaping data to long format...")
  long <- df |>
    dplyr::filter(Variable %in% top_vars$Variable) |>
    tidyr::pivot_longer(
      cols = tidyselect::all_of(year_cols),
      names_to = "Year",
      values_to = "percent"
    ) |>
    dplyr::mutate(
      Variable = factor(
        Variable,
        levels = top_vars |> dplyr::arrange(dplyr::desc(mean_pct)) |> dplyr::pull(Variable)
      ),
      YearNum  = as.integer(stringr::str_replace(Year, "^Y", "")),
      Year     = factor(Year, levels = year_cols, labels = as.character(year_nums))
    ) |>
    dplyr::arrange(Variable, YearNum)

  plot_title <- paste0(alpha_code, " Variable Contributions (", min(year_nums), "-", max(year_nums), ")")

  # ---- [7.5] colorblind-safe palette + non-color cues ------------------------
  levels_var <- levels(long$Variable)
  K <- length(levels_var)
  pal_base <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#000000", "#F0E442",
    "#4477AA", "#66CCEE", "#228833", "#CCBB44",
    "#EE6677", "#AA3377", "#BBBBBB"
  )
  more_needed <- max(0, K - length(pal_base))
  if (more_needed > 0) {
    extra <- grDevices::hcl(
      h = seq(15, 375, length.out = more_needed + 1)[-1],
      c = 65, l = 45
    )
    pal_base <- c(pal_base, extra)
  }
  pal <- stats::setNames(pal_base[seq_len(K)], levels_var)
  pal_fill <- grDevices::adjustcolor(pal, alpha.f = 0.18)  # ribbons
  lt_vals  <- rep(c("solid","solid","dashed","dotdash"), length.out = K)
  pch_vals <- c(rep(21, min(K, 12)), rep(22, max(0, K - 12)))
  names(lt_vals)  <- levels_var
  names(pch_vals) <- levels_var

  # ---- [8/13] heatmap --------------------------------------------------------
  message(">>> [8/13] Creating heatmap...")
  p_heatmap <- ggplot2::ggplot(long, ggplot2::aes(x = Year, y = Variable, fill = percent)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_continuous(name = "% Contribution") +
    ggplot2::labs(title = plot_title, x = "Year", y = "Variable") +
    trend_theme_base() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  # print(p_heatmap)

  # ---------- Build Bayesian stats FIRST to create lab_map ----------
  message(">>> [9/13] Preparing data splits...")
  split_list <- split(long, long$Variable)

  # ---- [10/13] Bayesian regression (rstanarm) + ribbons + PD + CIs + ROPE ---
  message(">>> [10/13] Running Bayesian linear regression per variable (rstanarm)...")
  bayes_list <- lapply(names(split_list), function(v) {
    d <- split_list[[v]]
    d <- d[!is.na(d$percent) & !is.na(d$YearNum), ]
    if (nrow(d) >= 2 && length(unique(d$YearNum)) >= 2) {
      fit <- rstanarm$stan_glm(
        percent ~ YearNum, data = d, family = gaussian(),
        prior = rstanarm$normal(0, 10, autoscale = TRUE),
        prior_intercept = rstanarm$normal(0, 10, autoscale = TRUE),
        chains = 4, iter = 2000, seed = 1234, refresh = 0
      )
      post <- as.data.frame(as.matrix(fit))
      a_draw <- post[["(Intercept)"]]
      b_draw <- post[["YearNum"]]
      X <- stats::model.matrix(~ YearNum, data = d)
      fit_draws <- as.matrix(X %*% t(cbind(a_draw, b_draw)))
      var_fit <- apply(fit_draws, 2, stats::var)
      sigma2  <- post[["sigma"]]^2
      r2_draw <- var_fit / (var_fit + sigma2)

      summarize <- function(x) c(mean = mean(x), sd = stats::sd(x),
                                 q2.5 = as.numeric(stats::quantile(x, 0.025)),
                                 q97.5 = as.numeric(stats::quantile(x, 0.975)))
      PD <- function(x) 100 * max(mean(x > 0), mean(x < 0))

      sa <- summarize(a_draw); sb <- summarize(b_draw); sr2 <- summarize(r2_draw)
      pd_intercept <- PD(a_draw); pd_slope <- PD(b_draw)
      eq_mean <- paste0("y = ", sprintf("%.3f", sa["mean"]),
                        ifelse(sb["mean"] >= 0, " + ", " - "),
                        sprintf("%.3f", abs(sb["mean"])), "x; R^2 = ", sprintf("%.3f", sr2["mean"]))

      # ROPE metrics
      rope_slope_pct <- if (!is.null(rope_slope)) rope_prob(b_draw, rope_slope) else NA_real_
      rope_slope_dec <- rope_decision(sb["q2.5"], sb["q97.5"], rope_slope)
      rope_int_pct   <- if (!is.null(rope_intercept)) rope_prob(a_draw, rope_intercept) else NA_real_
      rope_int_dec   <- if (!is.null(rope_intercept)) rope_decision(sa["q2.5"], sa["q97.5"], rope_intercept) else NA_character_

      list(
        variable = v,
        n_points = nrow(d),
        a_draw = a_draw, b_draw = b_draw,
        summary = data.frame(
          Variable          = v,
          n_points          = nrow(d),
          intercept_mean    = sa["mean"], intercept_sd = sa["sd"],
          slope_mean        = sb["mean"], slope_sd     = sb["sd"],
          r2_mean           = sr2["mean"], r2_sd       = sr2["sd"],
          intercept_q2.5    = sa["q2.5"],  intercept_q97.5 = sa["q97.5"],
          slope_q2.5        = sb["q2.5"],  slope_q97.5     = sb["q97.5"],
          r2_q2.5           = sr2["q2.5"], r2_q97.5        = sr2["q97.5"],
          intercept_ci_lower= sa["q2.5"],  intercept_ci_upper = sa["q97.5"],
          slope_ci_lower    = sb["q2.5"],  slope_ci_upper     = sb["q97.5"],
          r2_ci_lower       = sr2["q2.5"], r2_ci_upper        = sr2["q97.5"],
          pd_intercept      = pd_intercept,
          pd_slope          = pd_slope,
          # ROPE (slope)
          rope_slope_lower  = if (!is.null(rope_slope)) rope_slope[1] else NA_real_,
          rope_slope_upper  = if (!is.null(rope_slope)) rope_slope[2] else NA_real_,
          rope_slope_pct    = rope_slope_pct,
          rope_slope_decision = rope_slope_dec,
          # ROPE (intercept; optional)
          rope_intercept_lower = if (!is.null(rope_intercept)) rope_intercept[1] else NA_real_,
          rope_intercept_upper = if (!is.null(rope_intercept)) rope_intercept[2] else NA_real_,
          rope_intercept_pct   = rope_int_pct,
          rope_intercept_decision = rope_int_dec,
          equation_mean     = eq_mean,
          stringsAsFactors  = FALSE
        )
      )
    } else {
      list(
        variable = v, n_points = nrow(d),
        a_draw = NULL, b_draw = NULL,
        summary = data.frame(
          Variable=v, n_points=nrow(d),
          intercept_mean=NA_real_, intercept_sd=NA_real_,
          slope_mean=NA_real_, slope_sd=NA_real_,
          r2_mean=NA_real_, r2_sd=NA_real_,
          intercept_q2.5=NA_real_, intercept_q97.5=NA_real_,
          slope_q2.5=NA_real_, slope_q97.5=NA_real_,
          r2_q2.5=NA_real_, r2_q97.5=NA_real_,
          intercept_ci_lower=NA_real_, intercept_ci_upper=NA_real_,
          slope_ci_lower=NA_real_,     slope_ci_upper=NA_real_,
          r2_ci_lower=NA_real_,        r2_ci_upper=NA_real_,
          pd_intercept=NA_real_, pd_slope=NA_real_,
          rope_slope_lower=NA_real_, rope_slope_upper=NA_real_, rope_slope_pct=NA_real_, rope_slope_decision=NA_character_,
          rope_intercept_lower=NA_real_, rope_intercept_upper=NA_real_, rope_intercept_pct=NA_real_, rope_intercept_decision=NA_character_,
          equation_mean=NA_character_,
          stringsAsFactors = FALSE
        )
      )
    }
  })
  stats_br_df <- do.call(rbind, lapply(bayes_list, `[[`, "summary"))

  # Posterior predictions per year
  message(">>> Preparing Bayesian mean lines and 95% ribbons...")
  year_grid <- sort(unique(long$YearNum))
  br_pred_df <- do.call(rbind, lapply(bayes_list, function(item) {
    if (is.null(item$a_draw)) return(NULL)
    a <- item$a_draw; b <- item$b_draw
    out <- lapply(year_grid, function(y) {
      draws <- a + b * y
      c(mean = mean(draws),
        q2.5 = as.numeric(stats::quantile(draws, 0.025)),
        q97.5 = as.numeric(stats::quantile(draws, 0.975)),
        YearNum = y)
    })
    out <- do.call(rbind, out)
    data.frame(
      Variable = item$variable,
      YearNum = as.numeric(out[, "YearNum"]),
      mean = as.numeric(out[, "mean"]),
      q2.5 = as.numeric(out[, "q2.5"]),
      q97.5 = as.numeric(out[, "q97.5"]),
      stringsAsFactors = FALSE
    )
  }))
  if (!is.null(br_pred_df) && nrow(br_pred_df) > 0) {
    br_pred_df$Variable <- factor(br_pred_df$Variable, levels = levels(long$Variable))
  }

  # ---- Legend linewidths and label stars/signs -------------------------------
  lw_df <- stats_br_df[, c("Variable", "pd_slope", "slope_mean")]
  lw_df$Variable <- factor(lw_df$Variable, levels = levels(long$Variable))

  # normalize PD to percent if needed (accepts 0-1 or 0-100)
  pdp <- suppressWarnings(as.numeric(lw_df$pd_slope))
  pdp <- ifelse(is.na(pdp), NA_real_, ifelse(pdp <= 1, 100 * pdp, pdp))

  # tiered stars: *** >=95, ** >=90 & <95, * >=85 & <90
  stars <- ifelse(!is.na(pdp) & pdp >= 95,  "***",
                  ifelse(!is.na(pdp) & pdp >= 90,  "**",
                         ifelse(!is.na(pdp) & pdp >= 85,  "*", "")))

  # (+)/(-) only for starred variables, based on slope_mean
  signs <- ifelse(stars != "" & !is.na(lw_df$slope_mean),
                  ifelse(lw_df$slope_mean > 0, "(+)",
                         ifelse(lw_df$slope_mean < 0, "(-)", "")),
                  "")
  star_sign <- paste0(stars, signs)

  # Label map for legends (markdown-capable if ggtext is available)
  base_names <- as.character(lw_df$Variable)
  if (has_ggtext) {
    lab_vals <- ifelse(star_sign != "",
                       paste0("<b>", base_names, "</b>", star_sign),
                       base_names)
  } else {
    lab_vals <- paste0(base_names, star_sign)
  }
  lab_map <- stats::setNames(lab_vals, as.character(lw_df$Variable))

  # Thicker lines for any starred tier (PD >=85) on BR plots
  lw_vals <- ifelse(!is.na(pdp) & pdp >= 85, 2.75, 0.50)
  lw_map  <- stats::setNames(lw_vals, as.character(lw_df$Variable))

  # ---- [11/13] Create plots (lab_map now defined) ----------------------------
  message(">>> [11/13] Creating line plot...")
  p_lines <- ggplot2::ggplot(
    long, ggplot2::aes(x = YearNum, y = percent, color = Variable,
                       group = Variable, linetype = Variable, shape = Variable)
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2.8, stroke = 0.7, fill = "white", show.legend = FALSE) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, labels = lab_map, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = lt_vals, labels = lab_map) +
    ggplot2::scale_shape_manual(values = pch_vals, guide = "none") +
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_base() +
    ggplot2::guides(
      color    = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA, alpha = 1)),
      linetype = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, alpha = 1))
    )
  # print(p_lines)

  message(">>> Creating frequentist LR plots...")
  base_lr_aes <- ggplot2::aes(x = YearNum, y = percent, color = Variable,
                              group = Variable, linetype = Variable)

  p_lr <- ggplot2::ggplot(long, base_lr_aes) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, labels = lab_map, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = lt_vals, labels = lab_map) +
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_base() +
    ggplot2::guides(
      color    = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA)),
      linetype = ggplot2::guide_legend(override.aes = list(linewidth = 2.4))
    )

  p_lr_pts <- ggplot2::ggplot(long, base_lr_aes) +
    ggplot2::geom_point(ggplot2::aes(shape = Variable), size = 2.8, stroke = 0.7, fill = "white", show.legend = FALSE) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, labels = lab_map, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = lt_vals, labels = lab_map) +
    ggplot2::scale_shape_manual(values = pch_vals, guide = "none") +
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_base() +
    ggplot2::guides(
      color    = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA)),
      linetype = ggplot2::guide_legend(override.aes = list(linewidth = 2.4))
    )
  # print(p_lr); print(p_lr_pts)

  message(">>> Creating Bayesian regression plots...")
  p_br <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = br_pred_df,
      ggplot2::aes(x = YearNum, ymin = q2.5, ymax = q97.5,
                   fill = Variable, group = Variable),
      alpha = 0.18, colour = NA, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = br_pred_df,
      ggplot2::aes(x = YearNum, y = mean, color = Variable,
                   group = Variable, linetype = Variable, linewidth = Variable)
    ) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE, labels = lab_map) +
    ggplot2::scale_fill_manual(values  = pal_fill, drop = FALSE, guide = "none") +
    ggplot2::scale_linetype_manual(values = lt_vals, guide = "none") +   # single legend (color)
    ggplot2::scale_linewidth_manual(values = lw_map, guide = "none") +   # linewidth not in legend
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_br() +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA))
    )

  p_br_pts <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = long,
      ggplot2::aes(x = YearNum, y = percent, color = Variable, shape = Variable),
      size = 2.6, stroke = 0.7, fill = "white", show.legend = FALSE
    ) +
    ggplot2::geom_ribbon(
      data = br_pred_df,
      ggplot2::aes(x = YearNum, ymin = q2.5, ymax = q97.5,
                   fill = Variable, group = Variable),
      alpha = 0.18, colour = NA, show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = br_pred_df,
      ggplot2::aes(x = YearNum, y = mean, color = Variable,
                   group = Variable, linetype = Variable, linewidth = Variable)
    ) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE, labels = lab_map) +
    ggplot2::scale_fill_manual(values  = pal_fill, drop = FALSE, guide = "none") +
    ggplot2::scale_linetype_manual(values = lt_vals, guide = "none") +
    ggplot2::scale_linewidth_manual(values = lw_map, guide = "none") +
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_br() +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA))
    )

  message(">>> Creating Bayesian regression plots (no ribbons)...")
  p_br_norib <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = br_pred_df,
      ggplot2::aes(x = YearNum, y = mean, color = Variable,
                   group = Variable, linetype = Variable, linewidth = Variable)
    ) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE, labels = lab_map) +
    ggplot2::scale_linetype_manual(values = lt_vals, guide = "none") +
    ggplot2::scale_linewidth_manual(values = lw_map, guide = "none") +
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_br() +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA))
    )

  p_br_pts_norib <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = long,
      ggplot2::aes(x = YearNum, y = percent, color = Variable, shape = Variable),
      size = 2.6, stroke = 0.7, fill = "white", show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = br_pred_df,
      ggplot2::aes(x = YearNum, y = mean, color = Variable,
                   group = Variable, linetype = Variable, linewidth = Variable)
    ) +
    ggplot2::scale_x_continuous(breaks = year_nums) +
    ggplot2::scale_color_manual(values = pal, drop = FALSE, labels = lab_map) +
    ggplot2::scale_linetype_manual(values = lt_vals, guide = "none") +
    ggplot2::scale_linewidth_manual(values = lw_map, guide = "none") +
    ggplot2::labs(title = plot_title, x = "Year", y = "% Contribution") +
    trend_theme_br() +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(linewidth = 2.4, shape = NA))
    )

  # ---- [12/13] Frequentist LR stats -----------------------------------------
  message(">>> [12/13] Computing frequentist LR stats...")
  stats_lr_df <- data.frame(
    Variable   = character(),
    n_points   = integer(),
    intercept  = double(),
    slope      = double(),
    r_squared  = double(),
    equation   = character(),
    stringsAsFactors = FALSE
  )
  if (length(split_list)) {
    stats_lr_list <- lapply(names(split_list), function(v) {
      d <- split_list[[v]]
      d <- d[!is.na(d$percent) & !is.na(d$YearNum), ]
      if (nrow(d) >= 2 && length(unique(d$YearNum)) >= 2) {
        ok <- try(stats::lm(percent ~ YearNum, data = d), silent = TRUE)
        if (!inherits(ok, "try-error")) {
          co  <- stats::coef(ok); a <- unname(co[1]); b <- unname(co[2])
          r2  <- unname(summary(ok)$r.squared)
          eqn <- paste0("y = ", sprintf("%.3f", a),
                        ifelse(b >= 0, " + ", " - "),
                        sprintf("%.3f", abs(b)), "x; R^2 = ", sprintf("%.3f", r2))
          return(data.frame(Variable=v, n_points=nrow(d),
                            intercept=a, slope=b, r_squared=r2, equation=eqn,
                            stringsAsFactors=FALSE))
        }
      }
      data.frame(Variable=v, n_points=nrow(d), intercept=NA_real_, slope=NA_real_,
                 r_squared=NA_real_, equation=NA_character_, stringsAsFactors=FALSE)
    })
    stats_lr_df <- do.call(rbind, stats_lr_list)
  }

  # ---- [13/13] save plots + stats + log -------------------------------------
  message(">>> [13/13] Saving plots and stats...")
  if (!dir.exists(vars_dir)) dir.create(vars_dir, recursive = TRUE, showWarnings = FALSE)

  dpi_out   <- 72
  height_px <- 750
  height_in <- height_px / dpi_out
  width_in  <- height_in * (9/5)

  ggplot2::ggsave(png_heatmap,      p_heatmap,      width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_lines,        p_lines,        width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_lr,           p_lr,           width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_lr_pts,       p_lr_pts,       width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_br,           p_br,           width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_br_pts,       p_br_pts,       width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_br_norib,     p_br_norib,     width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")
  ggplot2::ggsave(png_br_pts_norib, p_br_pts_norib, width = width_in, height = height_in, dpi = dpi_out, units = "in", bg = "white")

  if (nrow(stats_lr_df) > 0) {
    readr$write_csv(stats_lr_df, stats_lr_path)
  } else {
    readr$write_csv(stats_lr_df[0, ], stats_lr_path)
  }
  if (nrow(stats_br_df) > 0) {
    readr$write_csv(stats_br_df, stats_br_path)
  } else {
    readr$write_csv(stats_br_df[0, ], stats_br_path)
  }

  message("    Saved: ", png_heatmap,      " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_lines,        " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_lr,           " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_lr_pts,       " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_br,           " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_br_pts,       " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_br_norib,     " (1350x750 px @ 72 dpi)")
  message("    Saved: ", png_br_pts_norib, " (1350x750 px @ 72 dpi)")
  message("    Saved: ", stats_lr_path)
  message("    Saved: ", stats_br_path)

  # ---- logging ---------------------------------------------------------------
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  fmt_hms <- function(sec) { h <- floor(sec/3600); m <- floor((sec%%3600)/60); s <- sec%%60; sprintf("%02d:%02d:%05.2f", h, m, s) }

  per_year_lines <- vapply(year_nums, function(y) {
    vals <- long$percent[long$Year == as.character(y)]
    n_non_na <- sum(!is.na(vals))
    sv <- if (all(is.na(vals))) NA_real_ else sum(vals, na.rm = TRUE)
    rng <- if (all(is.na(vals))) c(NA_real_, NA_real_) else range(vals, na.rm = TRUE)
    sprintf("  %d  status=plotted vars=%d sum=%s   range=%s-%s   notes=-",
            y, n_non_na,
            ifelse(is.na(sv), "NA", sprintf("%.2f", sv)),
            ifelse(is.na(rng[1]), "NA", sprintf("%.2f", rng[1])),
            ifelse(is.na(rng[2]), "NA", sprintf("%.2f", rng[2])))
  }, FUN.VALUE = character(1))

  log_block <- paste0(
    "\n",
    "------------------------------------------------------------------------\n",
    "Processing summary (summarize_variable_contributions)\n",
    "Timestamp:          ", format(end_time, "%Y-%m-%d %H:%M:%S %Z"), "\n",
    "Alpha code:         ", alpha_code, "\n",
    "Years processed:    ", length(year_nums), "\n",
    "Eligibility rule:   variables present in >= 3 years (1980-2020)\n",
    "Eligible variables: ", sprintf("%d", n_available), "\n",
    "Top_n requested:    ", if (is.null(top_n)) "all" else as.character(n_requested), "\n",
    "Variables plotted:  ", n_use, "\n",
    "Completed:          1\n",
    "Failed:             0\n",
    "Total elapsed:      ", sprintf("%s (%.2f s)", fmt_hms(elapsed), elapsed), "\n",
    "Outputs saved:      heatmap+lines+lr+lr_points+br+br_noribbon+br_points_noribbon+lr_stats+br_stats\n",
    "Per-year status (sum of % contributions):\n",
    paste(per_year_lines, collapse = "\n"),
    "\n"
  )
  try({
    if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
    cat(log_block, file = log_path, append = TRUE)
  }, silent = TRUE)
  message("    Log appended: ", log_path)

  invisible(list(
    heatmap              = p_heatmap,
    lineplot             = p_lines,
    lrplot               = p_lr,
    lrplot_pts           = p_lr_pts,
    brplot               = p_br,
    brplot_pts           = p_br_pts,
    brplot_noribbon      = p_br_norib,
    brplot_pts_noribbon  = p_br_pts_norib,
    data                 = long,
    top_vars             = top_vars,
    stats_lr             = stats_lr_df,
    stats_br             = stats_br_df,
    paths                = list(
      heatmap      = png_heatmap,
      lines        = png_lines,
      lr           = png_lr,
      lr_pts       = png_lr_pts,
      br           = png_br,
      br_pts       = png_br_pts,
      br_norib     = png_br_norib,
      br_pts_norib = png_br_pts_norib,
      stats_lr     = stats_lr_path,
      stats_br     = stats_br_path
    )
  ))
}

## Global variable declarations for R CMD check (not part of roxygen docs)
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "Variable", "Year", "YearNum",
    "mean_pct", "percent",
    "q2.5", "q97.5",
    "years_present"
  ))
}
