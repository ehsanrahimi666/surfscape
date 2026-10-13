#' Extract decomposed change around point locations at several radii
#'
#' Summarises the output of [ss_decompose()] in circular buffers around field
#' sites, at several radii, and returns a tidy long table ready for
#' multi-scale modelling.
#'
#' @param x A `SpatRaster`, typically the output of [ss_decompose()].
#' @param pts An `sf` object, a `SpatVector` of points, a two-column matrix, or
#'   a `data.frame` with coordinate columns.
#' @param radii Numeric vector of buffer radii in map units.
#' @param fun Summary function applied within each buffer.
#' @param coords Column names holding x and y when `pts` is a `data.frame`.
#' @param id Optional column name in `pts` to carry through as a site
#'   identifier.
#'
#' @return A long `data.frame` with columns `site`, `radius`, `layer`, `value`.
#'
#' @section Buffer geometry:
#' Buffers are circular and defined by radius, so buffer area scales with
#' `radius^2`. When comparing extraction shapes rather than sizes, match on
#' area rather than on a linear dimension, otherwise any apparent effect of
#' shape is confounded with a threefold difference in the area sampled.
#'
#' @examples
#' set.seed(1)
#' z1 <- terra::rast(nrows = 50, ncols = 50, xmin = 0, xmax = 500,
#'                   ymin = 0, ymax = 500)
#' terra::values(z1) <- stats::rnorm(2500)
#' z2 <- terra::rast(z1); terra::values(z2) <- stats::rnorm(2500)
#' d <- ss_decompose(z1, z2, metrics = "sr1", w = 5)
#' pts <- data.frame(x = c(100, 300), y = c(100, 400))
#' ss_extract(d, pts, radii = c(50, 100))
#' @export
ss_extract <- function(x, pts, radii, fun = mean,
                       coords = c("x", "y"), id = NULL) {
  x <- .ss_rast(x)
  if (inherits(pts, "sf")) {
    v <- terra::vect(pts)
  } else if (inherits(pts, "SpatVector")) {
    v <- pts
  } else if (is.matrix(pts) || is.data.frame(pts)) {
    df <- as.data.frame(pts)
    if (!all(coords %in% names(df))) {
      if (ncol(df) >= 2L) names(df)[1:2] <- coords else
        stop("`pts` needs two coordinate columns.", call. = FALSE)
    }
    v <- terra::vect(df, geom = coords, crs = terra::crs(x))
  } else {
    stop("`pts` must be sf, SpatVector, matrix or data.frame.", call. = FALSE)
  }
  sid <- if (!is.null(id) && id %in% names(v)) as.character(v[[id]][, 1]) else
    as.character(seq_len(nrow(v)))

  out <- list()
  for (rad in radii) {
    b <- terra::buffer(v, width = rad)
    e <- terra::extract(x, b, fun = fun, na.rm = TRUE, ID = FALSE)
    for (ly in names(e)) {
      out[[length(out) + 1L]] <- data.frame(
        site = sid, radius = rad, layer = ly, value = e[[ly]],
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
