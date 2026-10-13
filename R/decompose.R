#' Decompose change in a continuous landscape surface
#'
#' Splits the change in a moving-window landscape metric between two dates
#' into three additive components: a change in the landscape-wide value
#' distribution, a redistribution of values among landscape windows, and a
#' rearrangement of values within windows. The three components sum to the
#' observed change at every cell, up to floating-point rounding.
#'
#' @param z1,z2 Single-layer `SpatRaster` objects for the two dates, with
#'   identical geometry. Cells missing in either date are dropped from both.
#' @param metrics Character vector of metric names, see [ss_metric_list()].
#' @param w Window size in cells: an odd integer or `c(nrow, ncol)`.
#' @param components Which components to return. Any subset of `"total"`,
#'   `"global"`, `"local"`, `"config"`.
#' @param range Value range shared by both dates, used for binning. Defaults to
#'   the range of the two dates combined, which is the correct choice for a
#'   temporal comparison.
#' @param nbins,nlevels,lacbox Metric settings, see [ss_metrics()].
#' @param filename,overwrite,chunk_rows,wopt Output and blocking control.
#'
#' @return A `SpatRaster` with one layer per metric and component, named
#'   `<metric>_<component>`.
#'
#' @section The decomposition:
#' For every window, the change in a metric between the two dates is split
#' into three parts that add up to the total change:
#' \describe{
#'   \item{global}{the part due to the change in the landscape-wide value
#'     distribution;}
#'   \item{local}{the part due to values moving between windows, so that a
#'     window contains a different set of values even when the landscape-wide
#'     distribution is unchanged;}
#'   \item{config}{the part due to the rearrangement of values inside the
#'     window.}
#' }
#' The components are obtained from counterfactual surfaces built by rank
#' transfer (see [ss_rank_transfer()]), which keeps the spatial ordering of
#' values of one date and gives it the values of the other date. Inside each
#' window, rank transfer between the two dated windows separates the change
#' due to the window's values (marginal part) from the change due to their
#' arrangement (`config`). Over the whole landscape, rank transfer between the
#' two dated surfaces gives the change due to the landscape-wide distribution
#' (`global`), and `local` is the marginal part minus `global`. Each
#' counterfactual comparison is averaged over the two possible orderings of
#' the dates (the two-factor Shapley value), so the result does not depend on
#' which date is taken first; where composition and configuration change
#' together, their interaction is shared equally between the components.
#' The formulas are given in the package vignette.
#'
#' Two limiting cases follow exactly from this construction. If the second
#' date is a strictly increasing transformation of the first, the spatial
#' ordering of values is unchanged everywhere and the `config` and `local`
#' components are zero. If the two dates contain the same values in a
#' different arrangement, the landscape-wide distribution is unchanged and
#' the `global` component is zero, while the `local` and `config`
#' components are generally not zero because values have moved between and
#' within windows.
#'
#' For a metric that depends only on the values present in the window (a
#' compositional metric, see [ss_benchmark()]) the `config` component is zero
#' at every cell and at every window size.
#'
#' @seealso [ss_benchmark()] to classify a metric before using it,
#'   [ss_scale_profile()] to run the decomposition across scales.
#'
#' @examples
#' set.seed(1)
#' z1 <- terra::rast(nrows = 40, ncols = 40)
#' terra::values(z1) <- stats::rnorm(1600)
#' z1 <- terra::focal(z1, 5, "mean")
#' z2 <- terra::rast(z1)
#' terra::values(z2) <- stats::rnorm(1600)
#' z2 <- terra::focal(z2, 9, "mean")
#' d <- ss_decompose(z1, z2, metrics = c("sq", "sr1"), w = 7)
#' names(d)
#' @export
ss_decompose <- function(z1, z2, metrics = c("sq", "sr1"), w = 7,
                         components = c("total", "global", "local", "config"),
                         range = NULL, nbins = 32L, nlevels = 16L,
                         lacbox = 3L, filename = "", overwrite = FALSE,
                         chunk_rows = NULL, wopt = list()) {

  components <- match.arg(components, several.ok = TRUE)
  zz <- .ss_check_pair(z1, z2)
  z1 <- zz$z1; z2 <- zz$z2
  w <- .ss_window(w)
  metrics <- unique(metrics)
  ids <- .ss_metric_id(metrics)
  nm <- length(metrics)

  if (is.null(range)) {
    mm <- cbind(terra::minmax(z1), terra::minmax(z2))
    range <- c(min(mm[1, ]), max(mm[2, ]))
  }
  pars <- .ss_pars(z1, nbins, range, nlevels, lacbox)
  halo <- (w[1] - 1L) %/% 2L

  ## within-window decomposition: total, local marginal, arrangement
  fun_d <- function(vlist, nr, nc)
    ss_decompose_cpp(vlist[[1]], vlist[[2]], nr, nc, w[1], w[2], ids, pars)
  dnames <- c(paste0(metrics, "_total"), paste0(metrics, "_marginal"),
              paste0(metrics, "_config"))
  D <- .ss_block_apply(list(z1, z2), halo, fun_d, dnames,
                       chunk_rows = chunk_rows, wopt = wopt)

  out <- list()
  if ("total" %in% components) out$total <- D[[paste0(metrics, "_total")]]
  if ("config" %in% components) out$config <- D[[paste0(metrics, "_config")]]

  if (any(c("global", "local") %in% components)) {
    c2m1 <- ss_rank_transfer(z2, z1)
    c1m2 <- ss_rank_transfer(z1, z2)
    fm <- function(z) {
      f <- function(vlist, nr, nc)
        ss_focal_metrics_cpp(vlist[[1]], nr, nc, w[1], w[2], ids, pars)
      .ss_block_apply(list(z), halo, f, metrics, chunk_rows = chunk_rows)
    }
    M1 <- fm(z1); M2 <- fm(z2); MA <- fm(c2m1); MB <- fm(c1m2)
    G <- 0.5 * ((M2 - MA) + (MB - M1))
    names(G) <- paste0(metrics, "_global")
    if ("global" %in% components) out$global <- G
    if ("local" %in% components) {
      L <- D[[paste0(metrics, "_marginal")]] - G
      names(L) <- paste0(metrics, "_local")
      out$local <- L
    }
  }

  ord <- unlist(lapply(metrics, function(m)
    paste0(m, "_", components)))
  res <- terra::rast(unname(out))
  res <- res[[match(ord, names(res))]]
  names(res) <- ord
  if (nzchar(filename)) {
    res <- terra::writeRaster(res, filename, overwrite = overwrite, wopt = wopt)
  }
  res
}

#' Scale profile of decomposed landscape change
#'
#' Runs [ss_decompose()] over a series of window sizes and summarises how the
#' three components of change vary with the size of the landscape window.
#'
#' @inheritParams ss_decompose
#' @param w_seq Vector of odd window sizes in cells.
#' @param rasters Logical. If `TRUE` the decomposition rasters for each scale
#'   are returned alongside the summary.
#' @param progress Logical, print progress.
#'
#' @return A `data.frame` with one row per window size and metric: `n`
#'   (cells), the mean absolute value of the total change and of each
#'   component (`mad_*`), `share_config`, the mean absolute configurational
#'   component divided by the sum of the mean absolute values of the three
#'   components, the correlations across cells between the compositional
#'   (global plus local) and configurational components (`r_comp_config`) and
#'   between the total change and the configurational component
#'   (`r_total_config`), and `sign_disagree`, the proportion of cells with a
#'   non-zero configurational component whose sign differs from that of the
#'   total change. Correlations are `NA` when a component is constant, and
#'   `sign_disagree` is `NA` when the configurational component is zero
#'   everywhere. When `rasters = TRUE`, a list with elements `summary` and
#'   `rasters`.
#'
#' @examples
#' set.seed(1)
#' z1 <- terra::rast(nrows = 60, ncols = 60)
#' terra::values(z1) <- stats::rnorm(3600)
#' z1 <- terra::focal(z1, 5, "mean")
#' z2 <- terra::rast(z1)
#' terra::values(z2) <- stats::rnorm(3600)
#' z2 <- terra::focal(z2, 9, "mean")
#' ss_scale_profile(z1, z2, metrics = "sr1", w_seq = c(5, 9))
#' @export
ss_scale_profile <- function(z1, z2, metrics = c("sq", "sr1"),
                             w_seq = c(3, 7, 15, 31), range = NULL,
                             nbins = 32L, nlevels = 16L, lacbox = 3L,
                             rasters = FALSE, progress = TRUE) {

  zz <- .ss_check_pair(z1, z2)
  if (is.null(range)) {
    mm <- cbind(terra::minmax(zz$z1), terra::minmax(zz$z2))
    range <- c(min(mm[1, ]), max(mm[2, ]))
  }
  out <- list(); rs <- list()
  for (w in w_seq) {
    d <- ss_decompose(zz$z1, zz$z2, metrics = metrics, w = w, range = range,
                      nbins = nbins, nlevels = nlevels, lacbox = lacbox)
    if (rasters) rs[[as.character(w)]] <- d
    for (m in metrics) {
      v <- terra::values(d[[paste0(m, c("_total", "_global", "_local",
                                        "_config"))]])
      v <- v[stats::complete.cases(v) & is.finite(rowSums(v)), , drop = FALSE]
      if (!nrow(v)) next
      tot <- v[, 1]; gl <- v[, 2]; lo <- v[, 3]; cf <- v[, 4]
      comp <- gl + lo
      tolc <- 1e-12 * max(abs(tot))
      keep <- abs(cf) > tolc
      out[[length(out) + 1L]] <- data.frame(
        w = w, metric = m, n = nrow(v),
        mad_total = mean(abs(tot)),
        mad_global = mean(abs(gl)),
        mad_local = mean(abs(lo)),
        mad_config = mean(abs(cf)),
        share_config = mean(abs(cf)) /
          (mean(abs(gl)) + mean(abs(lo)) + mean(abs(cf))),
        r_comp_config = .ss_cor(comp, cf),
        r_total_config = .ss_cor(tot, cf),
        sign_disagree = if (any(keep))
          mean(sign(tot[keep]) != sign(cf[keep])) else NA_real_,
        stringsAsFactors = FALSE)
    }
    if (progress) message("scale ", w, " done")
  }
  s <- do.call(rbind, out)
  if (rasters) list(summary = s, rasters = rs) else s
}

#' Pearson correlation, NA when either variable is constant
#' @keywords internal
#' @noRd
.ss_cor <- function(a, b) {
  if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
  stats::cor(a, b)
}

#' Decompose change in landscape-level metrics
#'
#' Splits the change in metrics computed over the whole surface between two
#' dates into a compositional component, due to the change in the value
#' distribution, and a configurational component, due to the change in the
#' spatial arrangement of values. At this level the landscape window is the
#' whole surface, so the local component of [ss_decompose()] does not arise.
#'
#' @param z1,z2 Single-layer `SpatRaster` objects for the two dates, with
#'   identical geometry and no missing values.
#' @param metrics Built-in metric names, see [ss_metric_list()].
#' @param fun Alternatively, a function that takes a `SpatRaster` and returns
#'   a named numeric vector, for metrics implemented elsewhere.
#' @param range,nbins,nlevels,lacbox Settings of the built-in metrics, see
#'   [ss_metrics()]. `range` defaults to the joint range of both dates.
#'
#' @return A `data.frame` with one row per metric: the values at both dates
#'   (`value1`, `value2`), `total` change, `composition` and `configuration`
#'   components (which sum to `total`), and `share_config`, the absolute
#'   configurational component divided by the sum of the absolute values of
#'   both components.
#'
#' @details
#' The components are averages over the two orderings of the dates of the
#' differences between the observed surfaces and the rank-transfer
#' counterfactuals (one date's arrangement with the other date's values), as
#' in [ss_decompose()]. If the second date is a strictly increasing
#' transformation of the first, the configurational component is zero; if
#' both dates contain the same values, the compositional component is zero.
#'
#' @examples
#' z1 <- ss_simulate(60, 60, range = 6, seed = 1)
#' z2 <- ss_simulate(60, 60, range = 12, seed = 2)
#' ss_decompose_landscape(z1, z2, metrics = c("sq", "sr1", "gcon"))
#' ss_decompose_landscape(z1, z2, fun = function(r)
#'   c(moran = terra::autocor(r, method = "moran", global = TRUE)))
#' @export
ss_decompose_landscape <- function(z1, z2, metrics = NULL, fun = NULL,
                                   range = NULL, nbins = 32L, nlevels = 16L,
                                   lacbox = 3L) {
  zz <- .ss_check_pair(z1, z2)
  z1 <- zz$z1
  z2 <- zz$z2
  v1 <- terra::values(z1)[, 1]
  v2 <- terra::values(z2)[, 1]
  if (anyNA(v1)) {
    stop("Landscape-level decomposition requires surfaces without missing ",
         "values.", call. = FALSE)
  }
  if (!is.null(metrics) && !is.null(fun)) {
    stop("Supply either `metrics` or `fun`, not both.", call. = FALSE)
  }
  if (is.null(fun) && is.null(metrics)) metrics <- ss_metric_list()$metric
  c21 <- ss_rank_transfer(z2, z1)
  c12 <- ss_rank_transfer(z1, z2)
  if (is.null(fun)) {
    ids <- .ss_metric_id(metrics)
    if (is.null(range)) range <- range(c(v1, v2))
    pars <- .ss_pars(z1, nbins, range, nlevels, lacbox)
    ev <- function(r) {
      out <- ss_whole_cpp(terra::values(r)[, 1], terra::nrow(r),
                          terra::ncol(r), ids, 0L, pars)[1, ]
      stats::setNames(out, metrics)
    }
  } else {
    ev <- function(r) {
      out <- fun(r)
      if (!is.numeric(out)) stop("`fun` must return a numeric vector.",
                                 call. = FALSE)
      if (is.null(names(out))) names(out) <- paste0("V", seq_along(out))
      out
    }
  }
  M1 <- ev(z1)
  M2 <- ev(z2)
  MA <- ev(c21)
  MB <- ev(c12)
  comp <- 0.5 * ((M2 - MA) + (MB - M1))
  conf <- 0.5 * ((MA - M1) + (M2 - MB))
  data.frame(metric = names(M1), value1 = unname(M1), value2 = unname(M2),
             total = unname(M2 - M1), composition = unname(comp),
             configuration = unname(conf),
             share_config = unname(abs(conf) / (abs(comp) + abs(conf))),
             stringsAsFactors = FALSE)
}
