#' Audit any metric function for configurational sensitivity
#'
#' The open counterpart of [ss_sensitivity()]. Instead of the metrics built into
#' \pkg{surfscape}, this audits an arbitrary R function that maps a
#' neighbourhood to one or more numbers, which makes it possible to test metrics
#' implemented in other packages (\pkg{geodiv}, \pkg{GLCMTextures},
#' \pkg{rasterdiv}) or defined by the user on exactly the same footing.
#'
#' @param x A single-layer `SpatRaster`.
#' @param fun A function taking a `w x w` numeric matrix and returning a numeric
#'   vector. If the vector is named, the names are used as metric names;
#'   otherwise names are generated. The function must return a vector of the
#'   same length for every window.
#' @param w Window size in cells (single odd integer, or `c(nrow, ncol)`).
#' @param n Number of windows to sample.
#' @param nperm Number of within-window permutations per sampled window.
#' @param seed Optional integer seed.
#' @param complete_only Logical. Sample only windows with no missing values
#'   (default `TRUE`), which is what most external metric implementations
#'   require.
#' @param max_tries Maximum number of candidate windows examined while searching
#'   for `n` complete ones.
#'
#' @return A `data.frame` with the same columns as [ss_sensitivity()]:
#'   `metric`, `n`, `csi`, `effect`, `exceed`, plus `invariant`, which is `TRUE`
#'   when the observed value equalled the permuted value in every window to
#'   within numerical tolerance. `invariant = TRUE` is the operational
#'   definition of an aspatial metric: it cannot detect spatial rearrangement,
#'   in this landscape or any other.
#'
#' @details
#' The permutation reshuffles the values inside the neighbourhood and leaves the
#' multiset of values untouched. Any metric that is a function of that multiset
#' alone is therefore exactly invariant, and both `csi` and `effect` are zero to
#' machine precision. The audit makes no assumption about how the metric is
#' implemented, so a published implementation can be tested directly rather than
#' reasoned about from its formula.
#'
#' @seealso [ss_sensitivity()] for the built-in metrics, which is much faster
#'   because it avoids the R callback.
#'
#' @examples
#' set.seed(1)
#' r <- terra::rast(nrows = 80, ncols = 80)
#' terra::values(r) <- stats::rnorm(6400)
#' r <- terra::focal(r, 5, "mean")
#'
#' f <- function(m) {
#'   v <- as.numeric(m)
#'   d <- v - mean(v)
#'   c(sd = sqrt(mean(d^2)),                              # aspatial
#'     lag1 = mean((m[, -ncol(m)] - m[, -1])^2))          # spatial
#' }
#' ss_sensitivity_fun(r, f, w = 9, n = 100, nperm = 19, seed = 1)
#' @export
ss_sensitivity_fun <- function(x, fun, w = 9, n = 300, nperm = 99,
                               seed = NULL, complete_only = TRUE,
                               max_tries = 50L) {

  x <- .ss_layer(x)
  w <- .ss_window(w)
  if (!is.function(fun)) stop("`fun` must be a function.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)

  nr <- terra::nrow(x); nc <- terra::ncol(x)
  h <- (w - 1L) %/% 2L
  if (nr <= 2 * h[1] || nc <= 2 * h[2]) {
    stop("The raster is too small for a ", w[1], " x ", w[2], " window.",
         call. = FALSE)
  }
  M <- terra::as.matrix(x, wide = TRUE)

  draw <- function(k) {
    i <- sample(seq.int(h[1] + 1L, nr - h[1]), k, replace = TRUE)
    j <- sample(seq.int(h[2] + 1L, nc - h[2]), k, replace = TRUE)
    cbind(i, j)
  }
  wins <- vector("list", n)
  got <- 0L; tries <- 0L
  while (got < n && tries < max_tries) {
    ij <- draw(n - got)
    for (r in seq_len(nrow(ij))) {
      win <- M[(ij[r, 1] - h[1]):(ij[r, 1] + h[1]),
               (ij[r, 2] - h[2]):(ij[r, 2] + h[2]), drop = FALSE]
      if (complete_only && anyNA(win)) next
      got <- got + 1L
      wins[[got]] <- win
      if (got == n) break
    }
    tries <- tries + 1L
  }
  wins <- wins[seq_len(got)]
  if (got < 10L) {
    stop("Fewer than ten usable windows were found; check the missing-data ",
         "mask or reduce `w`.", call. = FALSE)
  }

  first <- fun(wins[[1]])
  if (!is.numeric(first)) stop("`fun` must return a numeric vector.",
                               call. = FALSE)
  nm <- length(first)
  mnames <- names(first)
  if (is.null(mnames)) mnames <- paste0("V", seq_len(nm))

  obs <- matrix(NA_real_, got, nm)
  mu <- matrix(NA_real_, got, nm)
  sdv <- matrix(NA_real_, got, nm)
  maxdev <- rep(0, nm)

  for (i in seq_along(wins)) {
    win <- wins[[i]]
    o <- as.numeric(fun(win))
    P <- matrix(NA_real_, nperm, nm)
    for (p in seq_len(nperm)) {
      pw <- win
      pw[] <- sample(as.numeric(win))
      P[p, ] <- as.numeric(fun(pw))
    }
    obs[i, ] <- o
    mu[i, ] <- colMeans(P, na.rm = TRUE)
    sdv[i, ] <- apply(P, 2, stats::sd, na.rm = TRUE)
    dv <- abs(o - apply(P, 2, max, na.rm = TRUE))
    dv2 <- abs(o - apply(P, 2, min, na.rm = TRUE))
    maxdev <- pmax(maxdev, pmax(dv, dv2), na.rm = TRUE)
  }

  scale0 <- pmax(apply(abs(obs), 2, mean, na.rm = TRUE), .Machine$double.eps)
  res <- lapply(seq_len(nm), function(j) {
    o <- obs[, j]; m <- mu[, j]; s <- sdv[, j]
    ok <- is.finite(o) & is.finite(m)
    if (sum(ok) < 10L) {
      return(data.frame(metric = mnames[j], n = sum(ok), csi = NA_real_,
                        effect = NA_real_, exceed = NA_real_,
                        invariant = NA, stringsAsFactors = FALSE))
    }
    o <- o[ok]; m <- m[ok]; s <- s[ok]
    vo <- stats::var(o)
    csi <- if (vo <= .Machine$double.eps) 0 else
      max(0, 1 - suppressWarnings(stats::cor(o, m))^2)
    if (!is.finite(csi)) csi <- 0
    sdo <- sqrt(vo)
    eff <- if (sdo <= 0) 0 else mean(abs(o - m)) / sdo
    tol <- 1e-9 * max(abs(o), na.rm = TRUE)
    exc <- if (max(s, na.rm = TRUE) <= tol) 0 else
      mean(abs(o - m) > pmax(1.96 * s, tol), na.rm = TRUE)
    data.frame(metric = mnames[j], n = length(o), csi = csi, effect = eff,
               exceed = exc,
               invariant = maxdev[j] <= 1e-9 * scale0[j],
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  res[order(-res$csi), ]
}
