#' Simulate a continuous neutral landscape
#'
#' Generates a stationary Gaussian random field on a regular grid by circulant
#' embedding. Such fields are the continuous counterpart of neutral landscape
#' models: their spatial structure is set by a single parameter, the practical
#' range, and their value distribution can be replaced by any other through
#' rank transfer without altering the spatial arrangement. They provide
#' landscapes with known properties for controlled tests of landscape metrics.
#'
#' @param nrow,ncol Grid size in cells.
#' @param range Practical range in cells, the distance at which the
#'   correlation between two cells falls to 0.05.
#' @param model Correlation model: `"exponential"`, `exp(-3 h / range)`, or
#'   `"gaussian"`, `exp(-3 (h / range)^2)`, where `h` is the distance in cells.
#' @param values Optional numeric vector. If supplied, the simulated field is
#'   used only for its spatial ranks and receives the value distribution of
#'   `values`: the sorted values when `length(values)` equals the number of
#'   cells, otherwise the corresponding sample quantiles.
#' @param res Cell size in map units.
#' @param crs Coordinate reference system, passed to [terra::rast()].
#' @param seed Optional integer seed. The random-number state of the session
#'   is restored on exit.
#'
#' @return A single-layer `SpatRaster`. Without `values`, the field is
#'   standardised to mean 0 and standard deviation 1 over the grid.
#'
#' @details
#' The covariance is embedded in a torus of twice the grid size in each
#' direction and the field is drawn with two fast Fourier transforms (Wood and
#' Chan 1994; Dietrich and Newsam 1997). Negative eigenvalues of the embedding,
#' which can occur for the Gaussian model at long ranges, are set to zero; the
#' exponential model does not produce them at the ranges used in practice.
#'
#' @references
#' Wood ATA, Chan G (1994) Simulation of stationary Gaussian processes in
#' \eqn{[0, 1]^d}. *Journal of Computational and Graphical Statistics*
#' 3:409-432.
#'
#' Dietrich CR, Newsam GN (1997) Fast and exact simulation of stationary
#' Gaussian processes through circulant embedding of the covariance matrix.
#' *SIAM Journal on Scientific Computing* 18:1088-1107.
#'
#' @examples
#' z <- ss_simulate(100, 100, range = 10, seed = 1)
#' z
#' # the same arrangement with a skewed value distribution
#' zg <- ss_simulate(100, 100, range = 10, seed = 1,
#'                   values = stats::rgamma(10000, shape = 2))
#' @export
ss_simulate <- function(nrow, ncol, range,
                        model = c("exponential", "gaussian"),
                        values = NULL, res = 1, crs = "local", seed = NULL) {
  model <- match.arg(model)
  nrow <- as.integer(nrow)
  ncol <- as.integer(ncol)
  if (length(nrow) != 1L || length(ncol) != 1L || nrow < 2L || ncol < 2L) {
    stop("`nrow` and `ncol` must be single integers of at least 2.",
         call. = FALSE)
  }
  if (length(range) != 1L || !is.finite(range) || range <= 0) {
    stop("`range` must be a single positive number.", call. = FALSE)
  }
  f <- .ss_with_seed(seed, .ss_grf(nrow, ncol, range, model))
  v <- as.vector(t(f))
  if (!is.null(values)) {
    values <- values[is.finite(values)]
    if (!length(values)) stop("`values` has no finite elements.", call. = FALSE)
    rk <- rank(v, ties.method = "first")
    v <- if (length(values) == length(v)) {
      sort(values)[rk]
    } else {
      stats::quantile(values, (rk - 0.5) / length(v), names = FALSE)
    }
  }
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol * res,
                   ymin = 0, ymax = nrow * res, crs = crs)
  terra::values(r) <- v
  names(r) <- "z"
  r
}

#' @keywords internal
#' @noRd
.ss_grf <- function(nr, nc, range, model) {
  m1 <- 2L * nr
  m2 <- 2L * nc
  i <- 0:(m1 - 1L)
  j <- 0:(m2 - 1L)
  dy <- pmin(i, m1 - i)
  dx <- pmin(j, m2 - j)
  h <- sqrt(outer(dy^2, dx^2, "+"))
  cv <- if (model == "exponential") exp(-3 * h / range) else
    exp(-3 * (h / range)^2)
  lam <- Re(stats::fft(cv))
  lam[lam < 0] <- 0
  eps <- complex(real = stats::rnorm(m1 * m2), imaginary = stats::rnorm(m1 * m2))
  f <- Re(stats::fft(sqrt(lam / (m1 * m2)) * matrix(eps, m1, m2)))
  f <- f[seq_len(nr), seq_len(nc)]
  (f - mean(f)) / .ss_sd(f)
}

#' @keywords internal
#' @noRd
.ss_sd <- function(x) {
  x <- as.vector(x)
  sqrt(mean((x - mean(x))^2))
}

#' @keywords internal
#' @noRd
.ss_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
    old <- get(".Random.seed", envir = genv, inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = genv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = genv), add = TRUE)
  }
  set.seed(seed)
  expr
}
