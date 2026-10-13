#' Controlled changes of composition or configuration
#'
#' Applies one controlled change to a continuous surface. Each change alters
#' either the values present (composition) while keeping their spatial
#' ordering, or the spatial arrangement while keeping the values present.
#' These are the tests used by [ss_benchmark()], exposed so that they can be
#' inspected, mapped or combined.
#'
#' @param x A single-layer `SpatRaster`.
#' @param type The change to apply:
#'   \describe{
#'     \item{`"shift"`}{adds a constant, `amount` standard deviations of `x`,
#'       to every cell (composition change; arrangement unchanged).}
#'     \item{`"scale"`}{stretches values about the mean by the factor
#'       `amount`; the mean is unchanged and the spread changes (composition
#'       change; arrangement unchanged).}
#'     \item{`"reshape"`}{a rank-preserving change of distribution shape: the
#'       values are replaced by those of a gamma distribution with skewness
#'       `amount` (a normal distribution for `amount = 0`) rescaled to the mean
#'       and standard deviation of `x`. Mean, spread and the spatial ordering
#'       of values are unchanged (composition change; arrangement unchanged).}
#'     \item{`"permute"`}{randomly permutes all values across the surface
#'       (configuration change; values unchanged).}
#'     \item{`"rearrange"`}{reassigns the values of `x` to the spatial ranks of
#'       `arrangement`, or of a simulated surface with practical range
#'       `sim_range` (configuration change at the extent of `x`; values
#'       unchanged).}
#'     \item{`"target"`}{keeps the spatial ranks of `x` and gives it the value
#'       distribution of `target` (composition change; arrangement unchanged).}
#'   }
#' @param amount Size of the change for `"shift"` (default 0.5), `"scale"`
#'   (default 1.5) and `"reshape"` (default 2).
#' @param target A `SpatRaster` with the geometry of `x`, for `"target"`.
#' @param arrangement A `SpatRaster` with the geometry of `x`, for
#'   `"rearrange"`.
#' @param sim_range Practical range, in cells, of the surface simulated with
#'   [ss_simulate()] for `"rearrange"` when `arrangement` is not supplied.
#' @param seed Optional integer seed for `"permute"` and `"rearrange"`. The
#'   random-number state of the session is restored on exit.
#'
#' @return A `SpatRaster` with the geometry of `x`.
#'
#' @details
#' Ties in `x` are kept as ties by `"reshape"`, so the transformation is a
#' non-decreasing function of the original values and the spatial ordering of
#' values inside any window is unchanged.
#'
#' @examples
#' z <- ss_simulate(60, 60, range = 8, seed = 1)
#' zs <- ss_perturb(z, "shift")
#' zr <- ss_perturb(z, "reshape", amount = 2)
#' cor(terra::values(z)[, 1], terra::values(zr)[, 1], method = "spearman")
#' zp <- ss_perturb(z, "rearrange", sim_range = 20, seed = 2)
#' all.equal(sort(terra::values(z)[, 1]), sort(terra::values(zp)[, 1]))
#' @export
ss_perturb <- function(x, type = c("shift", "scale", "reshape", "permute",
                                   "rearrange", "target"),
                       amount = NULL, target = NULL, arrangement = NULL,
                       sim_range = NULL, seed = NULL) {
  x <- .ss_layer(x)
  type <- match.arg(type)
  if (is.null(amount)) {
    amount <- switch(type, shift = 0.5, scale = 1.5, reshape = 2, NA_real_)
  }
  r <- .ss_with_seed(seed, {
    switch(type,
      target = {
        if (is.null(target)) stop("`target` is required.", call. = FALSE)
        .ss_same_geom(x, target, "target")
        ss_rank_transfer(x, target)
      },
      rearrange = {
        if (is.null(arrangement)) {
          if (is.null(sim_range)) {
            stop("Supply `arrangement` or `sim_range`.", call. = FALSE)
          }
          arrangement <- ss_simulate(terra::nrow(x), terra::ncol(x),
                                     range = sim_range)
          arrangement <- terra::rast(x, vals = terra::values(arrangement))
        }
        .ss_same_geom(x, arrangement, "arrangement")
        ss_rank_transfer(arrangement, x)
      },
      {
        v <- terra::values(x)[, 1]
        out <- terra::rast(x)
        terra::values(out) <- .ss_perturb_values(v, type, amount)
        out
      })
  })
  names(r) <- paste0(names(x), "_", type)
  r
}

#' @keywords internal
#' @noRd
.ss_perturb_values <- function(v, type, amount) {
  ok <- !is.na(v)
  z <- v[ok]
  v[ok] <- switch(type,
    shift = z + amount * .ss_sd(z),
    scale = mean(z) + amount * (z - mean(z)),
    reshape = .ss_reshape(z, amount),
    permute = z[sample.int(length(z))],
    stop("Unknown type.", call. = FALSE))
  v
}

#' @keywords internal
#' @noRd
.ss_reshape <- function(z, skew) {
  u <- (rank(z, ties.method = "average") - 0.5) / length(z)
  q <- if (skew == 0) {
    stats::qnorm(u)
  } else if (skew > 0) {
    stats::qgamma(u, shape = 4 / skew^2)
  } else {
    -stats::qgamma(1 - u, shape = 4 / skew^2)
  }
  mean(z) + .ss_sd(z) * (q - mean(q)) / .ss_sd(q)
}

#' @keywords internal
#' @noRd
.ss_same_geom <- function(x, y, arg) {
  if (!inherits(y, "SpatRaster") || terra::nlyr(y) != 1L ||
      !terra::compareGeom(x, y, stopOnError = FALSE)) {
    stop("`", arg, "` must be a single-layer SpatRaster with the geometry ",
         "of `x`.", call. = FALSE)
  }
  invisible(TRUE)
}
