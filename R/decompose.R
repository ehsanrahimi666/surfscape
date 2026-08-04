#' Decompose change in a continuous landscape surface
#'
#' Splits the change in any continuous landscape metric between two dates into
#' three additive components: a shift in the landscape-wide value distribution,
#' a spatial redistribution of values among neighbourhoods, and a rearrangement
#' of values within neighbourhoods. The three components sum exactly to the
#' observed change at every cell.
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
#' Let \eqn{M} be a metric evaluated on a neighbourhood, \eqn{W_1} and
#' \eqn{W_2} the neighbourhood at the two dates, and let \eqn{A} denote the
#' rank transfer that keeps one arrangement and substitutes the other value
#' multiset. Symmetrising over the two possible orderings (the two-factor
#' Shapley solution) gives
#'
#' \deqn{\Delta M = \underbrace{\Delta_{\mathrm{global}}}_{\text{distribution}}
#'  + \underbrace{\Delta_{\mathrm{local}}}_{\text{redistribution}}
#'  + \underbrace{\Delta_{\mathrm{config}}}_{\text{arrangement}}}
#'
#' where the arrangement component is obtained by transferring ranks *within*
#' each window, the global component by transferring ranks over the *whole*
#' landscape, and the local component as the remainder. The identity holds
#' exactly, not approximately: if the two dates differ only by a monotone
#' transformation of values the arrangement and local components are exactly
#' zero, and if they share an identical value multiset the global component is
#' exactly zero.
#'
#' For every metric with `arrangement = FALSE` in [ss_metric_list()] the
#' configurational component is exactly zero by construction, at every scale.
#' That is the point: such metrics cannot detect spatial rearrangement, however
#' large it is.
#'
#' @seealso [ss_sensitivity()] to audit a metric before using it,
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
#' Runs [ss_decompose()] over a ladder of window sizes and summarises how the
#' three components of change scale. The scale at which the configurational
#' component peaks answers a question the total change cannot: at what grain
#' did the landscape actually rearrange.
#'
#' @inheritParams ss_decompose
#' @param w_seq Vector of odd window sizes in cells.
#' @param rasters Logical. If `TRUE` the decomposition rasters for each scale
#'   are returned alongside the summary.
#' @param progress Logical, print progress.
#'
#' @return A `data.frame` with one row per window size and metric, containing
#'   the mean absolute magnitude of each component, the configurational share,
#'   the correlation between components across cells, and the proportion of
#'   cells where the sign of the total change differs from the sign of the
#'   configurational change. When `rasters = TRUE`, a list with elements
#'   `summary` and `rasters`.
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
      out[[length(out) + 1L]] <- data.frame(
        w = w, metric = m, n = nrow(v),
        mad_total = mean(abs(tot)),
        mad_global = mean(abs(gl)),
        mad_local = mean(abs(lo)),
        mad_config = mean(abs(cf)),
        share_config = mean(abs(cf)) / (mean(abs(comp)) + mean(abs(cf))),
        r_comp_config = suppressWarnings(stats::cor(comp, cf)),
        r_total_config = suppressWarnings(stats::cor(tot, cf)),
        sign_disagree = {
          tolc <- 1e-12 * max(abs(tot), na.rm = TRUE)
          keep <- abs(cf) > tolc
          if (!any(keep)) NA_real_ else
            mean(sign(tot[keep]) != sign(cf[keep]))
        },
        stringsAsFactors = FALSE)
    }
    if (progress) message("scale ", w, " done")
  }
  s <- do.call(rbind, out)
  if (rasters) list(summary = s, rasters = rs) else s
}
