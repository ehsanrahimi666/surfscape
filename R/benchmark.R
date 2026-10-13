#' Benchmark landscape metrics against controlled changes of composition and
#' configuration
#'
#' Tests how a continuous landscape metric responds to changes in the values
#' present (composition) and to changes in their spatial arrangement
#' (configuration). Landscape windows are sampled from a surface, each test is
#' applied, and the metric is evaluated on the same windows before and after
#' the change. The same design is applied to the built-in metrics and to any
#' metric function supplied by the user, including implementations from other
#' packages, so that they can be compared on identical windows.
#'
#' @param x A single-layer `SpatRaster`.
#' @param metrics Built-in metric names, see [ss_metric_list()]. Defaults to
#'   all of them when `fun` is not supplied.
#' @param fun Alternatively, a function that takes one landscape window and
#'   returns a named numeric vector of metric values, the same length for
#'   every window. If `fun` has an argument called `range`, it receives the
#'   value range used for binning under the chosen `protocol`; if it has an
#'   argument called `surface`, it receives the values that define that range
#'   (all values of both surfaces compared under `"fixed"`, of the surface the
#'   window is taken from under `"own"`), for binning schemes based on the
#'   value distribution, such as equal-probability quantization.
#' @param w Window size in cells, a single odd integer or `c(nrow, ncol)`, or
#'   `"landscape"` to treat the whole surface as a single window.
#' @param n Number of windows sampled at random among the complete windows
#'   (no missing values) of `x`. Ignored when `cells` is supplied or
#'   `w = "landscape"`.
#' @param nperm Number of random spatial permutations per window for the
#'   `"permute"` test.
#' @param tests Any of `"permute"`, `"shift"`, `"scale"`, `"reshape"`,
#'   `"rearrange"` and `"target"`; see Details.
#' @param shift,scale,skew Sizes of the `"shift"` (in standard deviations of
#'   `x`), `"scale"` (stretch factor) and `"reshape"` (target skewness) tests.
#' @param target A `SpatRaster` with the geometry of `x` whose value
#'   distribution is imposed on `x` by the `"target"` test, for example the
#'   same landscape at a second date.
#' @param arrangement A `SpatRaster` with the geometry of `x` whose spatial
#'   ranks receive the values of `x` in the `"rearrange"` test. If `NULL`, a
#'   surface is simulated with [ss_simulate()] using practical range
#'   `sim_range` (default: the window size).
#' @param sim_range See `arrangement`.
#' @param protocol How binned metrics (entropy, Rao's Q, lacunarity and grey
#'   level co-occurrence textures) set their value range when the surface
#'   changes. `"fixed"`: both surfaces of a comparison share their joint value
#'   range, as required when values are compared between dates. `"own"`: each
#'   surface uses its own value range.
#' @param cells Optional cell numbers of the window centres, for example the
#'   `cells` element of a previous result, so that different metric sets can
#'   be benchmarked on the same windows.
#' @param seed Optional integer seed. The random-number state of the session
#'   is restored on exit.
#' @param tol Relative numerical tolerance. A test leaves a metric invariant
#'   when no window changes by more than `tol` times the largest absolute value
#'   of the metric over all evaluated windows (observed, permuted and changed
#'   by any test), which absorbs floating-point rounding only.
#' @param alpha Significance level of the per-window permutation test.
#' @param fun_input Whether `fun` receives each window as a numeric
#'   `"matrix"` or as a `"raster"` (`SpatRaster` with the window's geometry and
#'   coordinate reference system).
#' @param nbins,nlevels,lacbox Settings of the built-in metrics, see
#'   [ss_metrics()].
#' @param details Logical. If `TRUE`, the per-window values are returned.
#'
#' @return An object of class `ss_benchmark`, a list with
#' \describe{
#'   \item{summary}{one row per metric: `category`, `configuration_response`
#'     (response to `"permute"`), `composition_response` (largest response
#'     among the composition tests), `signif` (proportion of windows in which
#'     the observed arrangement differs significantly from random
#'     arrangements) and, for each test, whether the metric was invariant.
#'     Responses to tests to which the metric is invariant are reported as 0;
#'     the unrounded values are in `responses`.}
#'   \item{responses}{one row per metric and test with `n` (windows),
#'     `invariant`, `max_abs` and `mean_abs` (largest and mean absolute change),
#'     `mean_change` (mean signed change), `response` (`mean_abs` divided by
#'     the standard deviation of the metric among windows; `NA` when the
#'     metric does not vary among windows or only one window is evaluated),
#'     and, for
#'     `"permute"`, `signif` and `ses` (mean standardised effect size,
#'     observed minus permutation mean over permutation standard deviation).}
#'   \item{cells}{the window centres.}
#'   \item{settings}{the settings used.}
#'   \item{values}{per-window values, when `details = TRUE`.}
#' }
#'
#' @details
#' The tests separate the two properties of a landscape:
#' \describe{
#'   \item{permute}{the values inside each window are randomly rearranged.
#'     Composition is unchanged and configuration is destroyed.}
#'   \item{shift}{a constant is added to every cell.}
#'   \item{scale}{values are stretched about their mean, which changes their
#'     spread but not their mean.}
#'   \item{reshape}{the shape of the value distribution changes while the
#'     mean, the spread and the spatial ordering of the values are kept.}
#'   \item{target}{the value distribution of `target` is imposed while the
#'     spatial ordering of `x` is kept.}
#'   \item{rearrange}{the values of the whole surface are moved to new
#'     locations while the landscape-wide value distribution is kept. Inside
#'     a window this changes both which values are present and how they are
#'     arranged, so it is reported but not used for classification.}
#' }
#' `shift`, `scale`, `reshape` and `target` change composition and keep the
#' spatial ordering of values. A metric is classified as
#' **compositional** when it is invariant to `permute` and responds to at
#' least one composition test or to `rearrange` (a metric that is invariant to
#' `permute` can only respond to `rearrange` through the change in the values
#' present), **configurational** when it responds to `permute` and is
#' invariant to every composition test, **mixed** when it responds to
#' `permute` and to at least one composition test, and **unresponsive** when
#' it is invariant to all tests.
#'
#' The response to a test is the mean absolute change of the metric divided by
#' its standard deviation among the sampled windows, which expresses the
#' change on the scale of the metric's natural variation across the landscape.
#' In the permutation test each window is tested separately, two-sided, by
#' comparing its observed value with its `nperm` permuted values; with
#' `nperm = 39` and `alpha = 0.05` a window is significant only when its
#' observed value lies beyond all permuted values in one direction (Hope
#' 1968). Values equal within the numerical tolerance count as ties. Fewer
#' than `2 / alpha - 1` permutations cannot give a significant result, and a
#' warning is issued.
#'
#' The tests differ in size: the permutation removes all spatial structure
#' inside the window, whereas the composition tests are moderate changes.
#' Responses are therefore comparable among metrics for the same test, not
#' between tests; the classification uses only whether a metric is invariant.
#'
#' @references
#' Hope ACA (1968) A simplified Monte Carlo significance test procedure.
#' *Journal of the Royal Statistical Society B* 30:582-598.
#'
#' @seealso [ss_perturb()] for the tests themselves, [ss_simulate()] for
#'   landscapes with known structure.
#'
#' @examples
#' z <- ss_simulate(80, 80, range = 10, seed = 1)
#' b <- ss_benchmark(z, metrics = c("sq", "sent", "sr1", "gcon"), w = 9,
#'                   n = 30, nperm = 39, seed = 1)
#' b
#'
#' # any metric function, evaluated on the same windows
#' f <- function(m) c(sd = stats::sd(as.vector(m)),
#'                    lag1 = mean((m[, -ncol(m)] - m[, -1])^2))
#' ss_benchmark(z, fun = f, w = 9, nperm = 39, cells = b$cells, seed = 1)
#' @export
ss_benchmark <- function(x, metrics = NULL, fun = NULL, w = 15, n = 100,
                         nperm = 39,
                         tests = c("permute", "shift", "scale", "reshape"),
                         shift = 0.5, scale = 1.5, skew = 2,
                         target = NULL, arrangement = NULL, sim_range = NULL,
                         protocol = c("fixed", "own"), cells = NULL,
                         seed = NULL, tol = 1e-9, alpha = 0.05,
                         fun_input = c("matrix", "raster"),
                         nbins = 32L, nlevels = 16L, lacbox = 3L,
                         details = FALSE) {

  x <- .ss_layer(x)
  protocol <- match.arg(protocol)
  fun_input <- match.arg(fun_input)
  all_tests <- c("permute", "shift", "scale", "reshape", "rearrange", "target")
  tests <- unique(match.arg(tests, all_tests, several.ok = TRUE))
  if (!is.null(metrics) && !is.null(fun)) {
    stop("Supply either `metrics` or `fun`, not both.", call. = FALSE)
  }
  if (!is.null(fun) && !is.function(fun)) {
    stop("`fun` must be a function.", call. = FALSE)
  }
  if (is.null(fun) && is.null(metrics)) metrics <- ss_metric_list()$metric
  if ("target" %in% tests && is.null(target)) {
    stop("`target` is required for the \"target\" test.", call. = FALSE)
  }
  landscape <- identical(w, "landscape")
  if (!landscape) w <- .ss_window(w)
  nperm <- as.integer(nperm)
  if ("permute" %in% tests && (length(nperm) != 1L || is.na(nperm) ||
                               nperm < 1L)) {
    stop("`nperm` must be a positive integer for the \"permute\" test.",
         call. = FALSE)
  }
  if (!"permute" %in% tests) nperm <- 0L
  if (nperm > 0L && 2 / (nperm + 1) > alpha) {
    warning("With nperm = ", nperm, " the smallest attainable two-sided ",
            "p-value is ", signif(2 / (nperm + 1), 3), ", larger than alpha; ",
            "no window can be significant. Use nperm >= ",
            ceiling(2 / alpha - 1), ".", call. = FALSE)
  }

  nr <- terra::nrow(x)
  nc <- terra::ncol(x)
  v <- terra::values(x)[, 1]
  if (landscape && anyNA(v)) {
    stop("`w = \"landscape\"` requires a surface without missing values.",
         call. = FALSE)
  }

  .ss_with_seed(seed, {
    if (landscape) {
      cells <- NA_integer_
    } else if (is.null(cells)) {
      cells <- .ss_sample_windows(x, w, n)
    } else {
      cells <- .ss_check_cells(x, w, cells)
    }

    pert <- list()
    for (tt in setdiff(tests, "permute")) {
      pert[[tt]] <- switch(tt,
        shift = .ss_perturb_values(v, "shift", shift),
        scale = .ss_perturb_values(v, "scale", scale),
        reshape = .ss_perturb_values(v, "reshape", skew),
        rearrange = {
          sr <- if (!is.null(sim_range)) sim_range else
            if (landscape) 10 else w[1]
          terra::values(ss_perturb(x, "rearrange", arrangement = arrangement,
                                   sim_range = sr))[, 1]
        },
        target = terra::values(ss_perturb(x, "target", target = target))[, 1])
    }

    ev <- if (is.null(fun)) {
      .ss_bench_builtin(x, v, pert, metrics, w, cells, nperm, protocol,
                        landscape, nbins, nlevels, lacbox)
    } else {
      .ss_bench_fun(x, v, pert, fun, w, cells, nperm, protocol, landscape,
                    fun_input)
    }

    resp <- .ss_bench_stats(ev, tests, tol, alpha)
    out <- list(summary = .ss_bench_summary(resp, tests),
                responses = resp,
                cells = cells,
                settings = list(w = if (landscape) "landscape" else w,
                                n = if (landscape) 1L else length(cells),
                                nperm = nperm, tests = tests,
                                shift = shift, scale = scale, skew = skew,
                                protocol = protocol, tol = tol, alpha = alpha,
                                nbins = nbins, nlevels = nlevels,
                                lacbox = lacbox))
    if (details) out$values <- ev
    class(out) <- "ss_benchmark"
    out
  })
}

#' @keywords internal
#' @noRd
.ss_sample_windows <- function(x, w, n) {
  bad <- terra::focal(is.na(x), w = c(w[1], w[2]), fun = "sum",
                      fillvalue = 1)
  cand <- which(terra::values(bad)[, 1] == 0)
  if (length(cand) < n) {
    stop("Only ", length(cand), " complete ", w[1], " x ", w[2],
         " windows are available; reduce `n` or `w`.", call. = FALSE)
  }
  cand[sample.int(length(cand), n)]
}

#' @keywords internal
#' @noRd
.ss_check_cells <- function(x, w, cells) {
  cells <- as.integer(cells)
  nr <- terra::nrow(x)
  nc <- terra::ncol(x)
  h <- (w - 1L) %/% 2L
  r <- (cells - 1L) %/% nc
  cc <- (cells - 1L) %% nc
  if (anyNA(cells) || any(r < h[1] | r >= nr - h[1] | cc < h[2] |
                          cc >= nc - h[2])) {
    stop("Some `cells` are too close to the edge for a ", w[1], " x ", w[2],
         " window.", call. = FALSE)
  }
  v <- terra::values(x)[, 1]
  for (k in seq_along(cells)) {
    if (anyNA(v[.ss_window_index(cells[k], nc, h)])) {
      stop("The window centred on cell ", cells[k],
           " contains missing values.", call. = FALSE)
    }
  }
  cells
}

#' Row-major indices of the cells of a window centred on `cell`
#' @keywords internal
#' @noRd
.ss_window_index <- function(cell, nc, h) {
  r <- (cell - 1L) %/% nc
  cc <- (cell - 1L) %% nc
  as.vector(t(outer((r - h[1]):(r + h[1]), (cc - h[2]):(cc + h[2]),
                    function(a, b) a * nc + b + 1L)))
}

#' @keywords internal
#' @noRd
.ss_vrange <- function(v) range(v, na.rm = TRUE)

#' Evaluate the built-in metrics for every test
#' @keywords internal
#' @noRd
.ss_bench_builtin <- function(x, v, pert, metrics, w, cells, nperm, protocol,
                              landscape, nbins, nlevels, lacbox) {
  ids <- .ss_metric_id(metrics)
  nr <- terra::nrow(x)
  nc <- terra::ncol(x)
  c0 <- if (landscape) NULL else as.integer(cells - 1L)
  ns <- if (landscape) 1L else length(cells)
  run <- function(vals, rng, np) {
    p <- .ss_pars(x, nbins, rng, nlevels, lacbox)
    M <- if (landscape) ss_whole_cpp(vals, nr, nc, ids, np, p) else
      ss_permute_cpp(vals, nr, nc, w[1], w[2], ids, c0, np, p)
    colnames(M) <- metrics
    M
  }
  rng0 <- .ss_vrange(v)
  M <- run(v, rng0, nperm)
  obs <- M[seq_len(ns), , drop = FALSE]
  perm <- NULL
  if (nperm > 0L) {
    perm <- array(M[-seq_len(ns), , drop = FALSE], c(ns, nperm, length(ids)),
                  dimnames = list(NULL, NULL, metrics))
  }
  tl <- lapply(names(pert), function(tt) {
    vt <- pert[[tt]]
    if (protocol == "fixed") {
      rj <- range(c(rng0, .ss_vrange(vt)))
      ob <- if (identical(rj, rng0)) obs else run(v, rj, 0L)
      list(obs = ob, pert = run(vt, rj, 0L))
    } else {
      list(obs = obs, pert = run(vt, .ss_vrange(vt), 0L))
    }
  })
  names(tl) <- names(pert)
  list(obs = obs, perm = perm, tests = tl, metrics = metrics)
}

#' Evaluate a user metric function for every test
#' @keywords internal
#' @noRd
.ss_bench_fun <- function(x, v, pert, fun, w, cells, nperm, protocol,
                          landscape, fun_input) {
  fmls <- names(formals(fun))
  has_range <- "range" %in% fmls
  has_surface <- "surface" %in% fmls
  nr <- terra::nrow(x)
  nc <- terra::ncol(x)
  if (landscape) {
    wd <- c(nr, nc)
    idx <- list(seq_len(nr * nc))
    tmpl <- list(terra::rast(x))
  } else {
    wd <- w
    h <- (w - 1L) %/% 2L
    idx <- lapply(cells, .ss_window_index, nc = nc, h = h)
    tmpl <- if (fun_input == "raster") {
      e <- terra::ext(x)
      rs <- terra::res(x)
      lapply(cells, function(cl) {
        r <- (cl - 1L) %/% nc
        cc <- (cl - 1L) %% nc
        terra::rast(nrows = w[1], ncols = w[2],
                    xmin = e$xmin + (cc - h[2]) * rs[1],
                    xmax = e$xmin + (cc + h[2] + 1L) * rs[1],
                    ymin = e$ymax - (r + h[1] + 1L) * rs[2],
                    ymax = e$ymax - (r - h[1]) * rs[2],
                    crs = terra::crs(x))
      })
    } else {
      vector("list", length(cells))
    }
  }
  ns <- length(idx)
  refinfo <- function(ref) {
    list(range = if (has_range) .ss_vrange(ref) else NULL,
         surface = if (has_surface) ref[!is.na(ref)] else NULL)
  }
  call1 <- function(vals, s, ri) {
    arg <- if (fun_input == "matrix") {
      matrix(vals, wd[1], wd[2], byrow = TRUE)
    } else {
      terra::setValues(tmpl[[s]], vals)
    }
    a <- list(arg)
    if (has_range) a$range <- ri$range
    if (has_surface) a$surface <- ri$surface
    out <- do.call(fun, a)
    if (!is.numeric(out)) stop("`fun` must return a numeric vector.",
                               call. = FALSE)
    out
  }

  first <- call1(v[idx[[1]]], 1L, refinfo(v))
  mnames <- names(first)
  if (is.null(mnames)) mnames <- paste0("V", seq_along(first))
  nm <- length(first)
  fill <- function(vals, ref, np = 0L) {
    ri <- refinfo(ref)
    obs <- matrix(NA_real_, ns, nm, dimnames = list(NULL, mnames))
    perm <- if (np > 0L) array(NA_real_, c(ns, np, nm),
                               dimnames = list(NULL, NULL, mnames)) else NULL
    for (s in seq_len(ns)) {
      vw <- vals[idx[[s]]]
      o <- call1(vw, s, ri)
      if (length(o) != nm) stop("`fun` returned ", length(o), " values for ",
                                "one window and ", nm, " for another.",
                                call. = FALSE)
      obs[s, ] <- o
      for (b in seq_len(np)) perm[s, b, ] <- call1(vw[sample.int(length(vw))],
                                                   s, ri)
    }
    list(obs = obs, perm = perm)
  }

  base <- fill(v, v, nperm)
  binned <- has_range || has_surface
  tl <- lapply(names(pert), function(tt) {
    vt <- pert[[tt]]
    same_values <- tt == "rearrange"
    if (protocol == "fixed" && !same_values) {
      ref <- c(v, vt)
      ob <- if (binned) fill(v, ref)$obs else base$obs
      list(obs = ob, pert = fill(vt, ref)$obs)
    } else {
      list(obs = base$obs, pert = fill(vt, if (same_values) v else vt)$obs)
    }
  })
  names(tl) <- names(pert)
  list(obs = base$obs, perm = base$perm, tests = tl, metrics = mnames)
}

#' Per-metric, per-test statistics
#' @keywords internal
#' @noRd
.ss_bench_stats <- function(ev, tests, tol, alpha) {
  rows <- list()
  sdw <- function(o) {
    o <- o[is.finite(o)]
    if (length(o) < 2L) return(NA_real_)
    s <- stats::sd(o)
    if (s > 0) s else NA_real_
  }
  atol <- function(o) {
    o <- o[is.finite(o)]
    m <- if (length(o)) max(abs(o)) else 0
    if (m == 0) tol else tol * m
  }
  for (m in ev$metrics) {
    o <- ev$obs[, m]
    P <- if ("permute" %in% tests) {
      matrix(ev$perm[, , m, drop = FALSE], nrow = dim(ev$perm)[1])
    } else NULL
    a <- atol(c(o, P, unlist(lapply(ev$tests, function(x)
      c(x$obs[, m], x$pert[, m])))))
    if ("permute" %in% tests) {
      D <- P - o
      fin <- is.finite(D)
      und <- sum(xor(is.finite(P), is.finite(o)))
      ok <- is.finite(o)
      lo <- rowSums(P <= o + a, na.rm = TRUE)
      hi <- rowSums(P >= o - a, na.rm = TRUE)
      np <- rowSums(is.finite(P))
      pv <- pmin(1, 2 * pmin(lo + 1, hi + 1) / (np + 1))
      mu <- rowMeans(P, na.rm = TRUE)
      sp <- apply(P, 1, function(z) if (sum(is.finite(z)) > 1L)
        stats::sd(z, na.rm = TRUE) else NA_real_)
      ses <- ifelse(is.finite(sp) & sp > a, (o - mu) / sp, NA_real_)
      mx <- if (any(fin)) max(abs(D[fin])) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        metric = m, test = "permute", n = sum(ok & np > 0L),
        invariant = und == 0L && (is.na(mx) || mx <= a),
        max_abs = mx,
        mean_abs = if (any(fin)) mean(abs(D[fin])) else NA_real_,
        mean_change = if (any(fin)) mean(D[fin]) else NA_real_,
        response = if (any(fin)) mean(abs(D[fin])) / sdw(o) else NA_real_,
        signif = mean(pv[ok & np > 0L] <= alpha * (1 + 1e-9)),
        ses = if (any(is.finite(ses))) mean(ses, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE)
    }
    for (tt in names(ev$tests)) {
      ot <- ev$tests[[tt]]$obs[, m]
      pt <- ev$tests[[tt]]$pert[, m]
      D <- pt - ot
      fin <- is.finite(D)
      und <- sum(xor(is.finite(pt), is.finite(ot)))
      mx <- if (any(fin)) max(abs(D[fin])) else NA_real_
      rows[[length(rows) + 1L]] <- data.frame(
        metric = m, test = tt, n = sum(fin),
        invariant = und == 0L && (is.na(mx) || mx <= a),
        max_abs = mx,
        mean_abs = if (any(fin)) mean(abs(D[fin])) else NA_real_,
        mean_change = if (any(fin)) mean(D[fin]) else NA_real_,
        response = if (any(fin)) mean(abs(D[fin])) / sdw(ot) else NA_real_,
        signif = NA_real_, ses = NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  res <- do.call(rbind, rows)
  rownames(res) <- NULL
  res
}

#' One row per metric with the classification
#' @keywords internal
#' @noRd
.ss_bench_summary <- function(resp, tests) {
  comp <- intersect(tests, c("shift", "scale", "reshape", "target"))
  ms <- unique(resp$metric)
  out <- lapply(ms, function(m) {
    r <- resp[resp$metric == m, , drop = FALSE]
    inv <- stats::setNames(r$invariant, r$test)
    rp <- stats::setNames(ifelse(r$invariant, 0, r$response), r$test)
    conf_inv <- if ("permute" %in% tests) inv[["permute"]] else NA
    comp_inv <- if (length(comp)) all(inv[comp]) else NA
    rear_inv <- if ("rearrange" %in% tests) inv[["rearrange"]] else TRUE
    category <- if (is.na(conf_inv) || is.na(comp_inv)) NA_character_ else
      if (conf_inv && comp_inv && rear_inv) "unresponsive" else
        if (conf_inv) "compositional" else
          if (comp_inv) "configurational" else "mixed"
    d <- data.frame(
      metric = m, category = category,
      configuration_response = if ("permute" %in% tests)
        unname(rp[["permute"]]) else NA_real_,
      composition_response = if (!length(comp)) NA_real_ else
        if (all(inv[comp])) 0 else .ss_max_finite(rp[comp][!inv[comp]]),
      signif = if ("permute" %in% tests)
        r$signif[r$test == "permute"] else NA_real_,
      stringsAsFactors = FALSE)
    for (tt in tests) d[[paste0(tt, "_invariant")]] <- unname(inv[[tt]])
    d
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' @keywords internal
#' @noRd
.ss_max_finite <- function(z) {
  z <- z[is.finite(z)]
  if (length(z)) max(z) else NA_real_
}

#' Combine benchmark results computed on separate sets of windows
#'
#' Combines [ss_benchmark()] results obtained with identical settings on
#' different windows, for example windows processed in parallel or windows
#' sampled from replicate landscapes, and recomputes all statistics on the
#' pooled windows.
#'
#' @param ... `ss_benchmark` objects computed with `details = TRUE`, or a
#'   single list of them.
#'
#' @return An object of class `ss_benchmark` with `details` included.
#'
#' @examples
#' z1 <- ss_simulate(60, 60, range = 8, seed = 1)
#' z2 <- ss_simulate(60, 60, range = 8, seed = 2)
#' b1 <- ss_benchmark(z1, metrics = c("sq", "sr1"), w = 9, n = 20,
#'                    seed = 1, details = TRUE)
#' b2 <- ss_benchmark(z2, metrics = c("sq", "sr1"), w = 9, n = 20,
#'                    seed = 2, details = TRUE)
#' ss_benchmark_combine(b1, b2)
#' @export
ss_benchmark_combine <- function(...) {
  bl <- list(...)
  if (length(bl) == 1L && !inherits(bl[[1]], "ss_benchmark")) bl <- bl[[1]]
  if (!length(bl) || !all(vapply(bl, inherits, logical(1), "ss_benchmark"))) {
    stop("All inputs must be `ss_benchmark` objects.", call. = FALSE)
  }
  if (!all(vapply(bl, function(b) !is.null(b$values), logical(1)))) {
    stop("All inputs must be computed with `details = TRUE`.", call. = FALSE)
  }
  keys <- c("w", "nperm", "tests", "shift", "scale", "skew", "protocol",
            "tol", "alpha", "nbins", "nlevels", "lacbox")
  s1 <- bl[[1]]$settings
  for (b in bl[-1]) {
    if (!identical(b$settings[keys], s1[keys]) ||
        !identical(b$values$metrics, bl[[1]]$values$metrics)) {
      stop("All inputs must use the same metrics and settings.", call. = FALSE)
    }
  }
  vals <- lapply(bl, `[[`, "values")
  ev <- list(obs = do.call(rbind, lapply(vals, `[[`, "obs")),
             perm = NULL, tests = list(), metrics = vals[[1]]$metrics)
  if (s1$nperm > 0L) {
    ns <- vapply(vals, function(v) dim(v$perm)[1], integer(1))
    P <- array(NA_real_, c(sum(ns), s1$nperm, length(ev$metrics)),
               dimnames = list(NULL, NULL, ev$metrics))
    at <- cumsum(c(0L, ns))
    for (k in seq_along(vals)) {
      P[(at[k] + 1L):at[k + 1L], , ] <- vals[[k]]$perm
    }
    ev$perm <- P
  }
  for (tt in names(vals[[1]]$tests)) {
    ev$tests[[tt]] <- list(
      obs = do.call(rbind, lapply(vals, function(v) v$tests[[tt]]$obs)),
      pert = do.call(rbind, lapply(vals, function(v) v$tests[[tt]]$pert)))
  }
  resp <- .ss_bench_stats(ev, s1$tests, s1$tol, s1$alpha)
  settings <- s1
  settings$n <- nrow(ev$obs)
  out <- list(summary = .ss_bench_summary(resp, s1$tests), responses = resp,
              cells = unlist(lapply(bl, `[[`, "cells")), settings = settings,
              values = ev)
  class(out) <- "ss_benchmark"
  out
}

#' @export
print.ss_benchmark <- function(x, digits = 3, ...) {
  s <- x$settings
  wtxt <- if (identical(s$w, "landscape")) "the whole landscape" else
    paste0(s$n, " windows of ", s$w[1], " x ", s$w[2], " cells")
  cat("<ss_benchmark> ", nrow(x$summary), " metrics, ", wtxt, "\n", sep = "")
  cat("tests: ", paste(s$tests, collapse = ", "),
      if ("permute" %in% s$tests) paste0(" (", s$nperm, " permutations)"),
      "; protocol: ", s$protocol, " range\n\n", sep = "")
  sm <- x$summary
  num <- vapply(sm, is.numeric, logical(1))
  sm[num] <- lapply(sm[num], signif, digits = digits)
  print(sm, row.names = FALSE)
  invisible(x)
}

#' @export
plot.ss_benchmark <- function(x, label = TRUE, ...) {
  s <- x$summary
  comp_tests <- intersect(names(s), paste0(c("shift", "scale", "reshape",
                                             "target"), "_invariant"))
  if (!"permute_invariant" %in% names(s) || !length(comp_tests)) {
    stop("The plot needs the \"permute\" test and at least one composition ",
         "test.", call. = FALSE)
  }
  xv <- s$composition_response
  yv <- s$configuration_response
  xinv <- apply(as.matrix(s[comp_tests]), 1, all)
  yinv <- s$permute_invariant
  ok <- (is.finite(xv) | xinv) & (is.finite(yv) | yinv)
  pos <- c(xv[ok & !xinv], yv[ok & !yinv])
  pos <- pos[is.finite(pos) & pos > 0]
  floor0 <- if (length(pos)) min(pos) / 10 else 1e-3
  px <- ifelse(xinv, floor0, xv)[ok]
  py <- ifelse(yinv, floor0, yv)[ok]
  s <- s[ok, , drop = FALSE]
  cols <- c(compositional = "#1b9e77", configurational = "#d95f02",
            mixed = "#7570b3", unresponsive = "grey50")
  graphics::plot(px, py, log = "xy", pch = 19, col = cols[s$category],
                 xlab = "Response to composition change (SD units)",
                 ylab = "Response to spatial permutation (SD units)", ...)
  graphics::abline(h = floor0, v = floor0, lty = 3, col = "grey60")
  graphics::mtext("invariant", side = 2, at = floor0, line = 0.2, cex = 0.7,
                  las = 1, adj = 1)
  graphics::mtext("invariant", side = 1, at = floor0, line = 0.2, cex = 0.7)
  if (label) graphics::text(px, py, s$metric, pos = 3, cex = 0.7, xpd = NA)
  graphics::legend("bottomright", legend = names(cols), col = cols, pch = 19,
                   bty = "n", cex = 0.8)
  invisible(x)
}
