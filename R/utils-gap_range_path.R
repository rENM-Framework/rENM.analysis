#' Resolve a species' GAP range shapefile
#'
#' Look up a species' GAP range directory in \code{data/_species.csv} and
#' return the path to the shapefile inside it.
#'
#' @details
#' The shapefile name is read from the \code{GAP.RANGE} column rather than
#' built from the alpha code. Most entries follow the pattern
#' \code{b<CODE>x_CONUS_Range_2001v1}, but the pattern is a convention and
#' not a rule: Mexican Spotted Owl, alpha code \code{MSOW}, is distributed
#' as \code{bSPOWl_CONUS_Range_2001v1}. The alpha code identifies the
#' species in this study and the range file identifies a GAP product, and
#' the two are set independently. The species table is the mapping between
#' them, so it is the only place either name should be read from.
#'
#' @param project_dir Character. rENM project directory.
#' @param alpha_code Character. Four-letter species code.
#'
#' @return Character. Path to the \code{.shp} file. Not checked for
#'   existence; callers report that in their own terms.
#'
#' @keywords internal
#' @noRd
.gap_range_path <- function(project_dir, alpha_code) {
  code <- toupper(trimws(alpha_code))

  species_csv <- file.path(project_dir, "data", "_species.csv")
  if (!file.exists(species_csv)) {
    stop("Species table not found at: ", species_csv, call. = FALSE)
  }

  sp       <- utils::read.csv(species_csv, stringsAsFactors = FALSE,
                              check.names = FALSE)
  norm     <- function(x) gsub("[^A-Z0-9]", "", toupper(x))
  col_norm <- norm(names(sp))
  a_idx    <- match("ALPHACODE", col_norm, nomatch = 0L)
  g_idx    <- match("GAPRANGE",  col_norm, nomatch = 0L)
  if (a_idx == 0L || g_idx == 0L) {
    stop("Species table must contain ALPHA.CODE and GAP.RANGE columns: ",
         species_csv, call. = FALSE)
  }

  row_idx <- which(toupper(trimws(sp[[names(sp)[a_idx]]])) == code)
  if (!length(row_idx)) {
    stop("No row found with alpha code '", code, "' in: ", species_csv,
         call. = FALSE)
  }

  gap_range <- trimws(sp[[names(sp)[g_idx]]][row_idx[1L]])
  if (!nzchar(gap_range) || is.na(gap_range)) {
    stop("GAP.RANGE is empty for alpha code '", code, "' in: ", species_csv,
         call. = FALSE)
  }

  file.path(project_dir, "data", "shapefiles", gap_range,
            paste0(gap_range, ".shp"))
}
