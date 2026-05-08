# rENM.analysis

![rENM](https://img.shields.io/badge/rENM-framework-blue)
![module](https://img.shields.io/badge/module-analysis-informational)

**Trend analysis and derived metrics for the rENM Framework**

## Overview

`rENM.analysis` computes the core analytical products of the rENM Framework.
It transforms modeled suitability outputs into interpretable trends, spatial
metrics, and ecological signals.

This package depends on `rENM.core` for project-directory resolution and
species metadata access. All functions accept an optional `project_dir`
argument; see `?rENM_project_dir` for configuration options.

## Key functions

| Function | Description |
|---|---|
| `find_suitability_trend()` | Compute per-cell Theil-Sen trends and Mann-Kendall statistics |
| `find_suitability_change_trend()` | Compute trends in suitability rate-of-change (acceleration/deceleration) |
| `find_trend_percentages()` | Summarize positive, negative, and zero trend proportions |
| `find_range_change_percentages()` | Quantify trend sign proportions within the GAP range |
| `find_weighted_centroid()` | Compute suitability-weighted spatial centroids |
| `analyze_weighted_centroids()` | Fit Bayesian trends to centroid latitude and longitude |
| `find_bioclimatic_velocity()` | Estimate climate-space displacement between 1980 and 2020 |
| `find_hot_spots()` | Identify cells with accelerating suitability declines |
| `create_hot_spot_map()` | Map hot spots within states intersecting the GAP range |
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

```r
# From GitHub
devtools::install_github("rENM-Framework/rENM.analysis")

# From a local source directory
devtools::install_local("rENM.analysis")
```

## Getting started

Set up a project directory and generate modeled suitability surfaces first
(see `rENM.model`), then run the analysis pipeline in order:

```r
library(rENM.analysis)

proj <- "/path/to/your/rENM/project"

# 1. Compute temporal suitability trends
find_suitability_trend("CASP")
find_suitability_change_trend("CASP")

# 2. Summarize trend statistics
find_trend_percentages("CASP")
find_range_change_percentages("CASP")

# 3. Spatial centroid analysis
find_weighted_centroid("CASP")
analyze_weighted_centroids("CASP")
find_bioclimatic_velocity("CASP")

# 4. Hot spots and state-level summaries
find_hot_spots("CASP")
create_hot_spot_map("CASP")
create_state_trend_analysis("CASP")

# 5. Variable contribution analysis
gather_variable_contributions("CASP")
summarize_variable_contributions("CASP")
```

For interactive work, configure the project directory once per session to
avoid passing it to every function:

```r
options(rENM.project_dir = "/path/to/your/rENM/project")

find_suitability_trend("CASP")
find_weighted_centroid("CASP")
# ...
```

## Analysis pipeline

```
find_suitability_trend()
        ↓
find_suitability_change_trend()
        ↓
find_trend_percentages()
find_range_change_percentages()
        ↓
find_weighted_centroid()
        ↓
analyze_weighted_centroids()
find_bioclimatic_velocity()
        ↓
find_hot_spots()
create_hot_spot_map()
create_state_trend_analysis()
        ↓
gather_variable_contributions()
        ↓
summarize_variable_contributions()
```

Trend rasters are written to `<run_dir>/Trends/suitability/`. Centroid and
velocity outputs go to `<run_dir>/Trends/centroids/`. Variable contribution
summaries go to `<run_dir>/Trends/variables/`. All functions append a
structured summary block to `<run_dir>/_log.txt`.

## Role in the rENM framework

`rENM.analysis` is the fourth stage in the pipeline:

```
rENM.core → rENM.data → rENM.model → rENM.analysis → rENM.ai → rENM.reports
```

It consumes the modeled suitability surfaces produced by `rENM.model` and
generates the quantitative trends, spatial metrics, and derived signals
consumed by `rENM.ai` and `rENM.reports`.

## License

See `LICENSE` for details.

---

**rENM Framework** — A modular system for reconstructing and analyzing
long-term ecological niche dynamics.
