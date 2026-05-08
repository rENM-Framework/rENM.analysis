utils::globalVariables(c(
  ## common ggplot / terra column names
  "x", "y", "value", "trend",
  ## dplyr / tidy eval
  ".data",
  ## summarize_variable_contributions
  "Variable", "Year", "YearNum",
  "mean_pct", "percent",
  "q2.5", "q97.5",
  "years_present"
))

#' @importFrom rENM.core rENM_project_dir get_species_info show_species show_variables
NULL

## Infix null-coalescing helper used by plot_suitability_change_trend and
## plot_trend_with_centroids.
`%||%` <- function(a, b) if (is.null(a)) b else a
