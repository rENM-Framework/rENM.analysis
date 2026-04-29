# rENM.analysis

![rENM](https://img.shields.io/badge/rENM-framework-blue) ![module](https://img.shields.io/badge/module-analysis-informational)

**Trend analysis and derived metrics for the rENM Framework**

## Overview

`rENM.analysis` computes the core analytical products of the rENM Framework. It transforms modeled suitability outputs into interpretable trends, spatial metrics, and ecological signals.

This package focuses on **trend detection, spatial analysis, and derived ecological metrics**.

## Role in the rENM Framework

Within the modular rENM ecosystem, `rENM.analysis`: - Quantifies **temporal trends in climatic suitability** - Computes **range change and hotspot dynamics** - Derives **centroid movement and bioclimatic velocity** - Summarizes **variable contributions and meta-patterns** - Produces analysis-ready datasets for reporting and AI synthesis

It converts modeled outputs into biologically meaningful signals.

## Key Functions

-   `find_suitability_trend()` — Estimate temporal suitability trends
-   `find_suitability_change_trend()` — Detect acceleration/deceleration
-   `find_range_change_percentages()` — Quantify range dynamics
-   `find_trend_percentages()` — Summarize trend categories
-   `find_hot_spots()` / `create_hot_spot_map()` — Identify key regions
-   `find_weighted_centroid()` — Compute spatial centroids
-   `analyze_weighted_centroids()` — Track centroid movement
-   `find_bioclimatic_velocity()` — Estimate climate-space movement
-   `gather_variable_contributions()` — Compile model contributions
-   `summarize_variable_contributions()` — Aggregate importance metrics
-   `plot_trend()` / `plot_trend_with_centroids()` — Visualization tools

## Installation

``` r
devtools::install_local("rENM.analysis")
```

## Example

``` r
library(rENM.analysis)

# compute trends
find_suitability_trend("CASP")

# analyze centroids
analyze_weighted_centroids("CASP")

# derive velocity
find_bioclimatic_velocity("CASP")
```

## Relationship to Other Packages

`rENM.analysis` provides the quantitative backbone for interpretation and reporting.

## License

See `LICENSE` for details.

------------------------------------------------------------------------

**rENM Framework**\
A modular system for reconstructing and analyzing long-term ecological niche dynamics.
