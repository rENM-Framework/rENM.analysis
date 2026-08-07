# rENM.analysis 0.2.0.9000

* Fixed `create_state_trend_analysis()` erroring with `[crop] extents do not
  overlap` for species whose trend raster's real data footprint (the
  modeling extent used to crop predictor variables) is smaller than the
  GAP.RANGE polygon. State selection now requires intersecting both
  GAP.RANGE and `_occs/extent.txt`, and the per-state raster crop falls back
  to `NA` trend percentages (rather than aborting) if a selected state still
  has no raster overlap.
* Fixed `create_hot_spot_map()` including states in its map labels and
  per-state stats CSV that intersect the GAP.RANGE polygon but fall outside
  the trend raster's actual modeling extent (same root cause as the
  `create_state_trend_analysis()` fix above). State selection now also
  requires intersecting `_occs/extent.txt`.

# rENM.analysis 0.1.0

* Initial release.
* Added `find_suitability_trend()` to compute per-cell Theil-Sen trends and
  Mann-Kendall statistics across the rENM time series.
* Added `find_suitability_change_trend()` to compute trends in suitability
  rate-of-change (acceleration/deceleration).
* Added `find_trend_percentages()` to summarize positive, negative, and zero
  trend proportions across the study extent.
* Added `find_range_change_percentages()` to quantify trend sign proportions
  within the USGS GAP range polygon.
* Added `find_weighted_centroid()` to compute suitability-weighted spatial
  centroids for each time bin.
* Added `analyze_weighted_centroids()` to fit Bayesian trends to centroid
  latitude and longitude over time.
* Added `find_bioclimatic_velocity()` to estimate climate-space displacement
  between 1980 and 2020.
* Added `find_hot_spots()` to identify cells with accelerating suitability
  declines.
* Added `create_hot_spot_map()` to map hot spots within states intersecting
  the GAP range.
* Added `create_state_trend_analysis()` to produce per-state suitability trend
  maps and statistics.
* Added `create_suitability_change_map()` as a full workflow function for
  suitability change trend visualization.
* Added `gather_variable_contributions()` to consolidate per-year variable
  importance files across the time series.
* Added `summarize_variable_contributions()` for frequentist and Bayesian trend
  analysis of variable contributions.
* Added `plot_trend()` to plot a climatic suitability trend raster.
* Added `plot_trend_with_centroids()` to plot a suitability trend raster with
  centroid shift overlay.
* Added `plot_suitability_change_trend()` to plot a suitability change trend
  raster.
* Added `save_trend_plot()` to save a suitability trend map to PNG.
* Added `save_trend_plot_with_centroids()` to save a suitability trend map with
  centroid arrow to PNG.
