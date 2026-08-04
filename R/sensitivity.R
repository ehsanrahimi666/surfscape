#' Audit how much a metric responds to spatial arrangement
#'
#' Quantifies the share of a metric's variation that comes from the spatial
#' arrangement of values rather than from the values themselves. Windows are
#' sampled from a real surface, the metric is evaluated on the observed
#' arrangement and on random permutations of the same values, and the two are
#' compared.
#'
#' @param x A single-layer `SpatRaster`.
#' @param metrics Character vector of metric names; defaults to every metric in
#'   [ss_metric_list()].
#' @param w Window size in cells.
#' @param n Number of windows to sample.
#' @param nperm Number of within-window permutations per sampled window.
#' @param seed Optional integer seed.
#' @param range,nbins,nlevels,lacbox Metric settings, see [ss_metrics()].
#'
#' @return A `data.frame` with one row per metric:
#' \describe{
#'   \item{csi}{configurational sensitivity index, \eqn{1 - R^2} of the
#'     observed metric regressed on its permutation mean. Zero means the metric
#'     is entirely determined by the value multiset and cannot see arrangement;
#'     one means arrangement explains all of its variation.}
#'   \item{effect}{mean absolute difference between the observed value and the
#'     permutation mean, in units of the between-window standard deviation.}
#'   \item{exceed}{proportion of windows whose observed value falls outside the
#'     central 95 percent of its own permutation distribution.}
#'   \item{arrangement}{whether the metric can respond to arrangement at all,
#'     as a matter of construction.}
#' }
#'
#' @details
#' `csi` separates two questions that the literature has treated as one. The
#' `arrangement` column states whether a metric *can* respond to configuration;
#' `csi` states how much it *does* in a particular landscape. A metric with
#' `arrangement = FALSE` always returns `csi = 0` and `effect = 0` exactly,
#' which is a useful check that the audit is behaving.
#'
#' @examples
#' set.seed(1)
#' r <- terra::rast(nrows = 60, ncols = 60)
#' terra::values(r) <- stats::rnorm(3600)
#' r <- terra::focal(r, 5, "mean")
#' ss_sensitivity(r, metrics = c("sq", "sr1", "gcon"), w = 9,
#'                n = 200, nperm = 19, seed = 1)
#' @export
ss_sensitivity <- function(x, metrics = ss_metric_list()$metric, w = 9,
                           n = 1000, nperm = 99, seed = NULL,
                           range = NULL, nbins = 32L, nlevels = 16L,
                           lacbox = 3L) {
  x <- .ss_layer(x)
  w <- .ss_window(w)
  metrics <- unique(metrics)
  ids <- .ss_metric_id(metrics)
  pars <- .ss_pars(x, nbins, range, nlevels, lacbox)
  if (!is.null(seed)) set.seed(seed)

  nr <- terra::nrow(x); nc <- terra::ncol(x)
  h <- (w - 1L) %/% 2L
  if (nr <= 2 * h[1] || nc <= 2 * h[2]) {
    stop("The raster is too small for a ", w[1], " x ", w[2], " window.",
         call. = FALSE)
  }
  rows <- sample(seq.int(h[1], nr - h[1] - 1L), n, replace = TRUE)
  cols <- sample(seq.int(h[2], nc - h[2] - 1L), n, replace = TRUE)
  cells <- as.integer(rows * nc + cols)

  v <- terra::values(x)[, 1]
  M <- ss_sensitivity_cpp(v, nr, nc, w[1], w[2], ids, cells,
                          as.integer(nperm), pars)
  nm <- length(metrics)
  obs <- M[, seq_len(nm), drop = FALSE]
  mu <- M[, nm + seq_len(nm), drop = FALSE]
  sdv <- M[, 2 * nm + seq_len(nm), drop = FALSE]

  info <- ss_metric_list()
  res <- lapply(seq_len(nm), function(j) {
    o <- obs[, j]; m <- mu[, j]; s <- sdv[, j]
    ok <- is.finite(o) & is.finite(m)
    if (sum(ok) < 10L) {
      return(data.frame(metric = metrics[j], n = sum(ok), csi = NA_real_,
                        effect = NA_real_, exceed = NA_real_,
                        stringsAsFactors = FALSE))
    }
    o <- o[ok]; m <- m[ok]; s <- s[ok]
    vo <- stats::var(o)
    csi <- if (vo <= .Machine$double.eps) 0 else {
      max(0, 1 - suppressWarnings(stats::cor(o, m))^2)
    }
    if (!is.finite(csi)) csi <- 0
    sdo <- sqrt(vo)
    eff <- if (sdo <= 0) 0 else mean(abs(o - m)) / sdo
    tol <- 1e-9 * max(abs(o), na.rm = TRUE)
    exc <- if (max(s, na.rm = TRUE) <= tol) 0 else
      mean(abs(o - m) > pmax(1.96 * s, tol), na.rm = TRUE)
    data.frame(metric = metrics[j], n = length(o), csi = csi,
               effect = eff, exceed = exc, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  res$family <- info$family[match(res$metric, info$metric)]
  res$arrangement <- info$arrangement[match(res$metric, info$metric)]
  res[order(-res$csi), c("metric", "family", "arrangement", "n", "csi",
                         "effect", "exceed")]
}

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
