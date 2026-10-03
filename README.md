# rENM.analysis

![rENM](https://img.shields.io/badge/rENM-framework-blue) ![module](https://img.shields.io/badge/module-analysis-informational)[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.20798280.svg)](https://doi.org/10.5281/zenodo.20798280)

**Trend analysis and derived metrics for the rENM Framework**

## Overview

`rENM.analysis` computes the core analytical products of the rENM Framework. It transforms modeled suitability outputs into interpretable trends, spatial metrics, and ecological signals.

This package depends on `rENM.core` for project-directory resolution and species metadata access. Functions find the project directory through `rENM.core::rENM_project_dir()`; see `?rENM_project_dir` for configuration options.

## Key functions

| Function | Description |
|------------------------------------|------------------------------------|
| `find_suitability_trend()` | Compute per-cell Theil-Sen trends and Mann-Kendall statistics |
| `find_suitability_change_trend()` | Compute the suitability change trend (trend in successive changes) |
| `find_trend_percentages()` | Summarize positive, negative, and zero trend proportions, overall and by region |
| `find_range_change_percentages()` | Quantify trend sign proportions within the GAP range |
| `find_weighted_centroid()` | Compute suitability-weighted spatial centroids |
| `analyze_weighted_centroids()` | Fit Bayesian trends to centroid latitude and longitude |
| `find_bioclimatic_velocity()` | Estimate the distance, bearing, and velocity of the suitability centroid's shift |
| `find_hot_spots()` | Identify cells with accelerating suitability declines |
| `create_hot_spot_map()` | Map hot spots within states intersecting the GAP range |
| `find_boundary_trend_statistics()` | Compare trends inside the GAP range with the 250 km ring around it |
| `create_state_trend_analysis()` | Per-state suitability trend map and statistics |
| `create_suitability_change_map()` | Full workflow for suitability change trend visualization |
| `gather_variable_contributions()` | Consolidate per-year variable importance files |
| `summarize_variable_contributions()` | Frequentist and Bayesian trend analysis of variable contributions |
| `plot_trend()` | Plot a climatic suitability trend raster |
| `plot_trend_with_centroids()` | Plot suitability trend with centroid shift overlay |
| `plot_suitability_change_trend()` | Plot suitability change trend raster |
| `save_trend_plot()` | Save a suitability trend map to PNG |
| `save_trend_plot_with_centroids()` | Save a suitability trend map with centroid arrow to PNG |

## Installation

``` r
# From GitHub
remotes::install_github("rENM-Framework/rENM.analysis")

# From a local source directory
remotes::install_local("rENM.analysis")
```

## Getting started

Set up a project directory and generate modeled suitability surfaces first (see `rENM.model`), then run the analysis pipeline in order:

``` r
library(rENM.analysis)

# set once per session, or set RENM_PROJECT_DIR in ~/.Renviron
options(rENM.project_dir = "/path/to/your/rENM/project")

# 1. Suitability trend and summary statistics
find_suitability_trend("CASP")
find_trend_percentages("CASP")
find_range_change_percentages("CASP")
create_state_trend_analysis("CASP")

# 2. Centroid shift and bioclimatic velocity
analyze_weighted_centroids("CASP")
find_bioclimatic_velocity("CASP")
save_trend_plot_with_centroids("CASP")

# 3. Variable contributions
gather_variable_contributions("CASP")
summarize_variable_contributions("CASP")

# 4. Change trend, hot spots, and boundary statistics
create_suitability_change_map("CASP")
find_trend_percentages("CASP", layer = "Suitability-Change-Trend")
create_hot_spot_map("CASP")
find_boundary_trend_statistics("CASP")
```

## Analysis pipeline

```         
find_suitability_trend()
        ↓
find_trend_percentages()
find_range_change_percentages()
create_state_trend_analysis()
        ↓
analyze_weighted_centroids()      ← calls find_weighted_centroid()
find_bioclimatic_velocity()
save_trend_plot_with_centroids()
        ↓
gather_variable_contributions()
summarize_variable_contributions()
        ↓
create_suitability_change_map()   ← calls find_suitability_change_trend()
find_trend_percentages(layer = "Suitability-Change-Trend")
        ↓
create_hot_spot_map()             ← calls find_hot_spots()
find_boundary_trend_statistics()
```

Trend rasters are written to `<run_dir>/Trends/suitability/`. Centroid and velocity outputs go to `<run_dir>/Trends/centroids/`. Variable contribution summaries go to `<run_dir>/Trends/variables/`. Most functions append a structured summary block to `<run_dir>/_log.txt`.

## Role in the rENM Framework

`rENM.analysis` is the fourth stage in the pipeline:

```         
rENM.core → rENM.data → rENM.model → rENM.analysis → rENM.ai → rENM.reports
```

It consumes the modeled suitability surfaces produced by `rENM.model` and generates the quantitative trends, spatial metrics, and derived signals consumed by `rENM.ai` and `rENM.reports`.

## License

See `LICENSE` for details.

------------------------------------------------------------------------

**rENM Framework** — A modular system for reconstructing and analyzing long-term ecological niche dynamics.
