# rENM.analysis 0.2.0.9000
* `summarize_variable_contributions()` records convergence diagnostics for
  each per-variable Bayesian fit in the BR-Stats file: `rhat_max`,
  `ess_bulk_min`, `ess_tail_min`, `n_divergent` and `fit_ok`, the last
  `TRUE` when R-hat is at most 1.01, both effective sample sizes are at
  least 400 and no transition diverged (Vehtari et al. 2021). These fits use
  three to nine points; in the 35 validation runs Stan warned of divergent
  transitions or low effective sample size with no record of the variable
  affected. Adds `posterior` and `rstan` to Imports.
* The `analyze_weighted_centroids()` help said its fixed seed (1234) served
  reproducibility of the plotted ribbons. It governs the MCMC draws, and so
  every slope, interval, PD and ROPE figure; the text now says so and that
  the seed is independent of the one passed to `rENM()`.

* `summarize_variable_contributions()` now ranks the top variables by their
  average contribution over all intervals, an interval in which a variable
  was not selected counting as zero. It ranked by the mean over selected
  intervals only (`mean_pct`), so a variable chosen in three intervals
  outranked one chosen in all nine at a similar level, and the report
  captions, which describe "average contribution to the overall time
  series", did not match. In the ten CASP runs swgnt was selected in 7 to 9
  of 9 intervals at about 10 percent yet reached the top 10 in only 4; on
  seed 42 it now ranks seventh, and qv2m, selected in three intervals,
  drops out. The three-interval eligibility rule is unchanged, since the
  per-variable regressions need at least three points.
* `gather_variable_contributions()` trims variable names. The PI-Ranked
  files pad each name with a space before the tab, so every name in
  `<CODE>-Variables-AllYears.csv` ended in one and exact matching failed.

* `summarize_variable_contributions()` no longer marks variables whose
  contribution slope has a probability of direction of at least 85 percent.
  Legend labels carried stars and a (+)/(-) sign and those lines were drawn
  thicker. Across six species and 35 seeded runs, only one
  flagged trend recurred in at least 80 percent of a species' runs (CASP
  bio8, declining in 9 of 10); the rest appeared at one seed and not the
  next, so marking them presented a single draw as a finding. The marking remains available as
  `mark_trends = TRUE`, for a later multi-run version that can mark only
  trends recurring across seeds. Every line is now drawn at one width.

* `find_hot_spots()` mapped the wrong cells. It masked a negative
  suitability trend with a positive change trend, `(A < 0) & (B > 0)`. The
  change trend is the Theil-Sen slope of successive differences, each a
  later suitability minus an earlier one, so a positive value on a
  declining cell means losses are easing. The mask selected declines that
  were slowing, the opposite of the documented definition. It now tests
  `(A < 0) & (B < 0)`: suitability falling, and each change more negative
  than the last. The error dates from the first commit and is in v0.1.0.
  On the CASP, EAME and GRRO pilot runs the old and new masks share no
  cells. Cells under the old mask lost suitability fastest before 2000 and
  then flattened (CASP -0.0026 per year in 1980-2000, -0.0002 in
  2000-2020); cells under the new mask were flat or rising before 2000 and
  lost suitability after it (CASP +0.00001, then -0.0022). State figures
  move substantially: Eastern Meadowlark in Florida from 7.23 to 92.29
  percent of range, Greater Roadrunner in Arizona from 61.74 to 26.10.
  Every consumer reads the mask file, so `create_hot_spot_map()`, its state
  statistics, and the hot-spot columns of
  `find_boundary_trend_statistics()` change with it. Range areas and all
  other boundary columns are unchanged. The map subtitle, the report
  caption and the GenAI prompts already described accelerating declines
  and now match what is mapped.
* The change trend is no longer described as "acceleration" and
  "deceleration". Those words were used in the signed sense, positive for
  acceleration, while hot spots use them in the magnitude sense, so a hot
  spot sat in an area the same report called decelerating. Help text and
  the change-trend map subtitle now say what the sign means: positive when
  successive changes grew more positive (gains strengthening or losses
  easing), negative when they grew more negative (gains fading or losses
  steepening). The `find_hot_spots()` help also named a nonexistent input
  file and the wrong output file; both are corrected.

* `find_boundary_trend_statistics()` gained `med_slope` and `med_abs_slope`
  per zone. `pos_pct` counts the area whose slope is above zero and says
  nothing about how far above, so a zone whose slopes hover around zero
  reports a precise-looking percentage near 50 that carries no signal.
  Across the twelve pilot species the ring's median magnitude ran from 4 to
  70 percent of the interior's: Loggerhead Shrike's ring is 20 percent
  positive at 70 percent of the interior's magnitude, a declining periphery,
  while Eastern Meadowlark's is 51 percent positive at a median slope of
  -3e-07, which over the study window is no change at all. The narrative for
  that species had already shipped saying "Both zones are majority positive".
  The percentages are unchanged; what was missing was the means to tell the
  two cases apart. Positivity was also found to fall steadily from the range
  core to its boundary rather than being uniform inside, so the interior
  figure is an average over a gradient; reporting that profile is left to a
  later release.
* `create_state_trend_analysis()` — `GAP.RANGE.POS.PCT` and
  `GAP.RANGE.NEG.PCT` are now taken over the state's range area that
  carries trend data, not over its whole range. The denominator summed
  every range cell while the numerators dropped cells where the trend
  raster is NA, so the two percentages fell short of 100 with nothing
  saying why. Across twelve species that affected 220 of 650 rows, worst
  case Western Meadowlark in Michigan at 96.153, and the shortfall tracks
  water: Michigan, New York and Rhode Island are the extremes. Not merely
  cosmetic, because `create_suitability_trend_summary_table()` appends the
  interior and ring rows from `find_boundary_trend_statistics()` into the
  same "Positive %" and "Negative %" columns, and those were already taken
  over the data-carrying area. One printed table was mixing two
  conventions, showing state rows summing to 96 above boundary rows
  summing to 100, in the same columns.
  Gained `GAP.RANGE.DATA.PCT`, the share of the state's range carrying
  trend data, so the coverage is visible rather than merely absent; it is
  the denominator the other two are taken over. A state whose range has no
  trend data at all now reports NA for both percentages rather than 0,
  matching the existing no-overlap branch and the boundary function: the
  shares are undefined there, not zero.
  Because no cell in any of the twelve trend rasters has a slope of
  exactly zero, positive and negative areas sum to the data area
  identically, so the new percentages sum to exactly 100 by construction
  and the change is a pure renormalization of the old pair. Existing
  values shift: 78 of 275 state rows move, 3 by more than a point, the
  largest 1.764 points. Reports built from these CSVs should be
  regenerated.
* `find_boundary_trend_statistics()` now skips with a warning when the
  buffered range polygon is absent, rather than stopping. That polygon is
  written only by `find_range_extent()`; an extent set by
  `find_occurrence_extent()` or `set_extent()` produces no range polygon, so
  there is no interior to difference a ring against and the comparison is
  undefined rather than merely unavailable. Stopping aborted the whole run
  from inside `rENM()`, discarding a finished model fit over an optional
  statistic, and it did so after the expensive modeling had completed.
  `create_suitability_trend_summary_table()` already guarded the same file and
  omitted the boundary block, so the report degrades consistently now.
* `find_boundary_trend_statistics()` — the run log entry never appeared. The
  log block referenced `gap_range`, which is defined nowhere; the variable
  holding the shapefile path is `gap_path`. Because the block is wrapped in
  `try(..., silent = TRUE)`, evaluating that argument threw while the vector
  was being assembled, which aborted the whole `cat()` before anything was
  written and discarded the error. The function produced its CSV correctly
  throughout, so nothing computed was affected; what was missing was the
  audit trail, on every run for every species since the function was added.
  Found by `R CMD check`, which reported it as a variable with no visible
  binding.
* `find_trend_percentages()` — `percent_positive`, `percent_negative` and
  `percent_zero` are now fractions of area rather than of cell count. The
  area columns were added earlier in this cycle without rebasing the
  percentages, so each report printed an area beside a percentage that was
  not that area over the stated denominator. Cells range from about 62 to
  72 km² across these extents, which moved the two bases apart by up to 1.4
  points: Cassin's Sparrow read 66.86% on counts against 65.43% on area.
  Area is the basis everywhere else the framework reports a percentage, so a
  single report was mixing both conventions, `find_boundary_trend_statistics()`
  being already area-based. The cell counts are unchanged and remain in the
  file as counts. Percentages in existing CSVs shift slightly; reports built
  from them should be regenerated.
* `find_suitability_change_trend()` — a failed ASCII Grid write no longer
  aborts the run. The `.asc` is a convenience copy; every consumer reads the
  GeoTIFF and falls back to `.asc` only when the `.tif` is absent, so losing
  it costs nothing. Greater Roadrunner died here fifteen minutes into its run
  with a complete, valid GeoTIFF and a complete, valid `.asc` both sitting on
  disk. Three faults in one block: the retry passed `gdal = "AAIGrid"`, but
  `gdal` takes creation options as `KEY=VALUE` and the driver belongs in
  `filetype`, so the retry could never have succeeded; the handler discarded
  the real condition and reported a missing AAIGrid driver, which was a guess
  and a wrong one, since the driver is present; and a secondary artifact was
  allowed to end the run. Now writes with `filetype`, reports the actual
  error, and warns.
* `find_suitability_trend()` — removed the `asc_out` list, which was built,
  never written, and then logged under "Outputs (AAIGrid)". Every run log for
  every species named three files that had never existed. The returned
  `paths` list loses its `asc` element for the same reason; nothing in the
  framework read it.
* `find_trend_percentages()` — gained a `layer` argument, defaulting to
  `"Suitability-Trend"` and also accepting `"Suitability-Change-Trend"`.
  The function was hardcoded to the suitability trend raster, so the
  accelerating and decelerating figures the narrative reports had no
  computed source anywhere in the pipeline. The argument names both the
  input raster and the output CSV, so the two runs do not collide.
* `find_trend_percentages()` — added `valid_area_km2`,
  `positive_area_km2`, `negative_area_km2` and `zero_area_km2`. The
  function reported cell counts and one total extent area, but the report
  asks for the area of each trend class, and on a lon/lat grid that cannot
  be recovered from a count and a mean cell size. `valid_area_km2` is also
  the denominator the existing percentages are taken against, which was
  not previously stated anywhere: `extent_area_km2` covers the full
  rectangle including cells with no data, and describing a percentage as a
  share of the modeled extent was therefore wrong.
* `create_hot_spot_map()` and `create_state_trend_analysis()` — the GAP range
  shapefile is now resolved through the `GAP.RANGE` column of
  `data/_species.csv` rather than built from the alpha code as
  `b<CODE>x_CONUS_Range_2001v1`. That pattern is a convention and not a
  rule, and both functions stopped with "GAP range shapefile not found" for
  any species whose range file departs from it. Mexican Spotted Owl, alpha
  code `MSOW`, ships as `bSPOWl_CONUS_Range_2001v1`. The more serious
  problem was that the three functions reading this file did not agree on
  how to find it: `find_boundary_trend_statistics()` already used the
  species table, so a single run could resolve the same input two different
  ways. All three now call one internal helper, `.gap_range_path()`.
* Fixed `create_hot_spot_map()` aborting with `replacement has 2 rows, data
  has 1` for species whose range boundary coincides with a state boundary.
  Cassin's Sparrow failed at the western tip of Texas. The range is clipped
  to the union of the states before use, so its edges already run along
  state lines; intersecting it again with an individual state lays two
  edges on top of each other, and GEOS reports that contact as points and
  slivers alongside the shared area. sf cannot fit the resulting geometry
  collection onto a one-row data frame. Intersections in this function now
  work on geometry rather than on `sf` objects, and keep only the polygonal
  part of any result, which is the only part that carries area. Whether a
  species trips this depends on where its range edge falls, so it was
  absent from the three species the clipping fix was developed against.
  Dissolving those parts is done with `terra::aggregate()` rather than
  `sf::st_union()`, because a sliver can carry a duplicate vertex and s2
  rejects that as a degenerate edge on the sphere. terra is planar and
  rasterizes in the raster's CRS in any case. Per-state areas are also
  summed across parts now, rather than taking the first, which was wrong
  for any multi-part clip.
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
* Added `find_boundary_trend_statistics()`, which compares trend behavior
  inside the GAP range against the buffer ring surrounding it — the zone
  between the raw range polygon and the buffered polygon written by
  `find_range_extent()`. A static range polygon cannot reveal a range edge
  that is moving; strong positive trends concentrated just outside the
  historic boundary are the signature of an advancing leading edge, and a
  range-based statistic never looks there. Reports area, positive and
  negative percentages, and hot spot area and percentage for each zone,
  writing `<CODE>-Suitability-Trend-Boundary-Statistics.csv`.

  Reported range-wide rather than per state, for two reasons: the ring
  extends beyond the range into states holding none of it, which have no
  counterpart row in the range-based table; and a statistic over few raster
  cells is unstable between model realizations, while the range-wide ring
  and interior figures reproduce closely. Areas use the same coverage
  weighting as `create_hot_spot_map()`, and percentages are taken over the
  area carrying trend data, which can be smaller than the zone where the
  model produced no prediction.
* `create_hot_spot_map()` now reports `range_area_km2` per state, measured as
  a coverage-weighted cell sum over the GAP range within that state — the
  same basis as `hotspot_area_km2`. The percentage column is computed against
  it and renamed from `hotspot_pct_of_state` to `hotspot_pct_of_range`.
  Because numerator and denominator are now sums over the same cells with the
  same weights, hot-spot area cannot exceed range area and the percentage
  cannot pass 100; previously the two were measured on different bases (a
  raster sum against a vector polygon area in EPSG:5070), which left a
  saturated state able to report above its own range. The vector area
  remains available as `GAP.RANGE.AREA` from `create_state_trend_analysis()`.
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
