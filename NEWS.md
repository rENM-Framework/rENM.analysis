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
* Fixed `create_hot_spot_map()` summing hot-spot area over an entire state
  instead of the state's GAP.RANGE portion, which could make reported
  hot-spot area exceed range area for states where GAP.RANGE covers only a
  small part of the state (observed for Idaho and Oregon in the Pinyon Jay
  pilot run). Hot-spot area is now masked to GAP.RANGE within the state
  before summation.
* Fixed `create_hot_spot_map()` counting each cell's full area when summing
  hot-spot area, so cells straddling the GAP.RANGE boundary contributed
  area lying outside the range. Cells are now weighted by the fraction of
  their area inside the polygon. The former behavior overshot in proportion
  to boundary length against cell size, which is negligible for a large
  range portion but substantial where a state holds only a handful of cells'
  worth of range; it was enough to report Oregon above 100% of its range
  area for Pinyon Jay. Reported hot-spot areas shift slightly for all
  states, and appreciably for small-range ones.

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
