#' @keywords internal
#' @noRd
.ss_rast <- function(x, arg = "x") {
  if (inherits(x, "SpatRaster")) return(x)
  if (is.matrix(x)) return(terra::rast(x))
  if (is.character(x) && length(x) == 1L) return(terra::rast(x))
  stop("`", arg, "` must be a SpatRaster, a matrix or a raster file path.",
       call. = FALSE)
}

#' @keywords internal
#' @noRd
.ss_layer <- function(x, arg = "x") {
  x <- .ss_rast(x, arg)
  if (terra::nlyr(x) != 1L) {
    stop("`", arg, "` must have exactly one layer; got ", terra::nlyr(x), ".",
         call. = FALSE)
  }
  x
}

#' @keywords internal
#' @noRd
.ss_window <- function(w) {
  if (is.matrix(w)) w <- dim(w)
  if (length(w) == 1L) w <- c(w, w)
  if (length(w) != 2L) stop("`w` must be a number or c(nrow, ncol).", call. = FALSE)
  w <- as.integer(round(w))
  if (any(w < 3L)) stop("`w` must be at least 3.", call. = FALSE)
  if (any(w %% 2L == 0L)) stop("`w` must be odd.", call. = FALSE)
  w
}

#' @keywords internal
#' @noRd
.ss_pars <- function(x, nbins = 32L, range = NULL, nlevels = 16L,
                     lacbox = 3L) {
  if (is.null(range)) {
    mm <- terra::minmax(x)
    if (!all(is.finite(mm))) {
      x <- terra::setMinMax(x)
      mm <- terra::minmax(x)
    }
    range <- c(min(mm[1, ]), max(mm[2, ]))
  }
  if (length(range) != 2L || !all(is.finite(range))) {
    stop("`range` must be a finite numeric vector of length 2.", call. = FALSE)
  }
  if (range[2] <= range[1]) range[2] <- range[1] + 1
  res <- terra::res(x)
  list(nbins = as.integer(nbins), gmin = range[1], gmax = range[2],
       nlevels = as.integer(nlevels), lacbox = as.integer(lacbox),
       resx = res[1], resy = res[2])
}

#' Global rank transfer
#'
#' Builds a counterfactual surface that keeps the spatial rank ordering of
#' `arrangement` and takes the value distribution of `values`. Only cells that
#' are non-missing in both inputs are used.
#'
#' @param arrangement,values Single-layer `SpatRaster` objects with identical
#'   geometry.
#' @return A `SpatRaster`.
#' @examples
#' r1 <- terra::rast(nrows = 20, ncols = 20)
#' r2 <- terra::rast(nrows = 20, ncols = 20)
#' set.seed(1)
#' terra::values(r1) <- stats::rnorm(400)
#' terra::values(r2) <- stats::rnorm(400, 2)
#' ss_rank_transfer(r2, r1)
#' @export
ss_rank_transfer <- function(arrangement, values) {
  arrangement <- .ss_layer(arrangement, "arrangement")
  values <- .ss_layer(values, "values")
  a <- terra::values(arrangement)[, 1]
  v <- terra::values(values)[, 1]
  ok <- !is.na(a) & !is.na(v)
  if (!any(ok)) stop("No cells are valid in both layers.", call. = FALSE)
  s <- sort(v[ok])
  out <- rep(NA_real_, length(a))
  out[ok] <- s[rank(a[ok], ties.method = "first")]
  r <- terra::rast(arrangement)
  terra::values(r) <- out
  names(r) <- "rank_transfer"
  r
}

#' @keywords internal
#' @noRd
.ss_blocks <- function(x, chunk_rows = NULL) {
  nr <- terra::nrow(x)
  if (is.null(chunk_rows)) {
    b <- terra::blocks(x, n = 4)
    return(list(row = b$row, nrows = b$nrows, n = b$n))
  }
  chunk_rows <- max(1L, as.integer(chunk_rows))
  starts <- seq.int(1L, nr, by = chunk_rows)
  list(row = starts, nrows = pmin(chunk_rows, nr - starts + 1L),
       n = length(starts))
}

#' Block-wise application of a window function to one or more aligned layers
#'
#' `fun(vlist, nr, nc)` receives a list of row-major numeric vectors, one per
#' input layer, and must return a matrix with `nr * nc` rows.
#' @keywords internal
#' @noRd
.ss_block_apply <- function(xlist, halo, fun, out_names,
                            filename = "", overwrite = FALSE,
                            chunk_rows = NULL, wopt = list()) {
  x1 <- xlist[[1]]
  nc <- terra::ncol(x1)
  nr_tot <- terra::nrow(x1)
  out <- terra::rast(x1, nlyrs = length(out_names))
  names(out) <- out_names

  for (x in xlist) terra::readStart(x)
  on.exit(for (x in xlist) try(terra::readStop(x), silent = TRUE), add = TRUE)

  bt <- terra::writeStart(out, filename = filename, overwrite = overwrite,
                          wopt = wopt)
  b <- if (is.null(chunk_rows)) bt else .ss_blocks(x1, chunk_rows)

  for (i in seq_len(b$n)) {
    r1 <- b$row[i]
    nrb <- b$nrows[i]
    rs <- max(1L, r1 - halo)
    re <- min(nr_tot, r1 + nrb - 1L + halo)
    nread <- re - rs + 1L
    vlist <- lapply(xlist, function(x)
      as.numeric(terra::readValues(x, rs, nread, 1L, nc)))
    res <- fun(vlist, nread, nc)
    top <- r1 - rs
    keep <- seq.int(top * nc + 1L, (top + nrb) * nc)
    terra::writeValues(out, res[keep, , drop = FALSE], r1, nrb)
  }
  out <- terra::writeStop(out)
  names(out) <- out_names
  out
}

#' @keywords internal
#' @noRd
.ss_check_pair <- function(z1, z2) {
  z1 <- .ss_layer(z1, "z1")
  z2 <- .ss_layer(z2, "z2")
  if (!terra::compareGeom(z1, z2, stopOnError = FALSE)) {
    stop("`z1` and `z2` must have identical geometry (extent, resolution, CRS).",
         call. = FALSE)
  }
  m <- !is.na(z1) & !is.na(z2)
  list(z1 = terra::mask(z1, m, maskvalues = c(0, NA)),
       z2 = terra::mask(z2, m, maskvalues = c(0, NA)))
}
