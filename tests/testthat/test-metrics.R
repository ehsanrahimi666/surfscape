mk <- function(nr = 31, nc = 37, seed = 1, smooth = 0) {
  set.seed(seed)
  r <- terra::rast(nrows = nr, ncols = nc, xmin = 0, xmax = nc, ymin = 0, ymax = nr)
  terra::values(r) <- stats::rnorm(terra::ncell(r))
  if (smooth > 0) r <- terra::focal(r, smooth, "mean", na.rm = TRUE)
  r
}
ref_win <- function(M, w, f) {
  h <- (w - 1) / 2; nr <- nrow(M); nc <- ncol(M)
  o <- matrix(NA_real_, nr, nc)
  for (i in (h + 1):(nr - h)) for (j in (h + 1):(nc - h))
    o[i, j] <- f(M[(i - h):(i + h), (j - h):(j + h)])
  o
}
ref_glcm <- function(W, L, gmin, gmax) {
  bw <- (gmax - gmin) / L
  q <- pmin(pmax(floor((W - gmin) / bw), 0), L - 1)
  P <- matrix(0, L, L); wr <- nrow(q); wc <- ncol(q)
  dr <- c(0, -1, -1, -1); dc <- c(1, 1, 0, -1)
  for (d in 1:4) for (r in 1:wr) for (cc in 1:wc) {
    r2 <- r + dr[d]; c2 <- cc + dc[d]
    if (r2 < 1 || r2 > wr || c2 < 1 || c2 > wc) next
    a <- q[r, cc] + 1; b <- q[r2, c2] + 1
    P[a, b] <- P[a, b] + 1; P[b, a] <- P[b, a] + 1
  }
  P <- P / sum(P)
  i <- row(P) - 1; j <- col(P) - 1
  mu <- sum(P * i)
  list(gcon = sum(P * (i - j)^2), gdis = sum(P * abs(i - j)),
       ghom = sum(P / (1 + (i - j)^2)), gasm = sum(P^2),
       gent = -sum(P[P > 0] * log(P[P > 0])),
       gcor = sum(P * (i - mu) * (j - mu)) / sum(P * (i - mu)^2))
}

test_that("composition-family metrics match a direct reference", {
  r <- mk(); M <- terra::as.matrix(r, wide = TRUE); w <- 5
  m <- ss_metrics(r, c("smean","sa","sq","ssk","sku","sp","sv","sz"), w = w)
  expect_equal(terra::as.matrix(m[["sq"]], wide = TRUE),
               ref_win(M, w, function(v) sqrt(mean((v - mean(v))^2))), tolerance = 1e-10)
  expect_equal(terra::as.matrix(m[["sa"]], wide = TRUE),
               ref_win(M, w, function(v) mean(abs(v - mean(v)))), tolerance = 1e-10)
  expect_equal(terra::as.matrix(m[["ssk"]], wide = TRUE),
               ref_win(M, w, function(v) { d <- v - mean(v); mean(d^3)/mean(d^2)^1.5 }),
               tolerance = 1e-8)
  expect_equal(terra::as.matrix(m[["sku"]], wide = TRUE),
               ref_win(M, w, function(v) { d <- v - mean(v); mean(d^4)/mean(d^2)^2 }),
               tolerance = 1e-8)
  expect_equal(terra::as.matrix(m[["sz"]], wide = TRUE),
               ref_win(M, w, function(v) max(v) - min(v)), tolerance = 1e-10)
})

test_that("lag-1 autocorrelation matches a direct reference", {
  r <- mk(smooth = 3); M <- terra::as.matrix(r, wide = TRUE); w <- 7
  m <- ss_metrics(r, "sr1", w = w)
  expect_equal(terra::as.matrix(m[["sr1"]], wide = TRUE),
    ref_win(M, w, function(W) {
      dx <- (W[, -ncol(W)] - W[, -1])^2; dy <- (W[-nrow(W), ] - W[-1, ])^2
      g <- sum(c(dx, dy)) / (2 * length(c(dx, dy)))
      1 - g / mean((W - mean(W))^2) }), tolerance = 1e-9)
})

test_that("GLCM metrics match a direct reference", {
  r <- mk(smooth = 3); M <- terra::as.matrix(r, wide = TRUE); w <- 7
  rng <- range(terra::values(r), na.rm = TRUE)
  m <- ss_metrics(r, c("gcon","gdis","ghom","gasm","gent","gcor"),
                  w = w, nlevels = 8, range = rng)
  for (nm in c("gcon","gdis","ghom","gasm","gent","gcor")) {
    expect_equal(terra::as.matrix(m[[nm]], wide = TRUE),
                 ref_win(M, w, function(W) ref_glcm(W, 8, rng[1], rng[2])[[nm]]),
                 tolerance = 1e-9, info = nm)
  }
})

test_that("gliding-box lacunarity matches a direct reference", {
  r <- mk(smooth = 3); M <- terra::as.matrix(r, wide = TRUE); w <- 7
  rng <- range(terra::values(r), na.rm = TRUE)
  m <- ss_metrics(r, "slac", w = w, lacbox = 3, range = rng)
  expect_equal(terra::as.matrix(m[["slac"]], wide = TRUE),
    ref_win(M, w, function(W) {
      b <- 3; ms <- c()
      for (i in 1:(nrow(W)-b+1)) for (j in 1:(ncol(W)-b+1))
        ms <- c(ms, sum(W[i:(i+b-1), j:(j+b-1)] - rng[1]))
      mean(ms^2) / mean(ms)^2 }), tolerance = 1e-9)
})

test_that("block-wise processing is exact and NA handling is strict", {
  r <- mk(nr = 41, nc = 23, smooth = 3)
  a <- ss_metrics(r, c("sq","sr1","gcon","srao"), w = 7)
  b <- ss_metrics(r, c("sq","sr1","gcon","srao"), w = 7, chunk_rows = 3)
  expect_identical(terra::values(a), terra::values(b))
  set.seed(2); terra::values(r)[sample(terra::ncell(r), 20)] <- NA
  m <- ss_metrics(r, "sq", w = 5)
  expect_gt(terra::global(is.na(m), "sum")[1, 1],
            terra::global(is.na(r), "sum")[1, 1])
})

test_that("the metric catalogue is internally consistent", {
  info <- ss_metric_list()
  expect_identical(info$binned, info$metric %in% c("sent", "srao", "slac",
                                                   "gcon", "gdis", "ghom",
                                                   "gasm", "gent", "gcor"))
  expect_error(ss_metrics(mk(), "nonsense"), "Unknown metric")
  expect_error(ss_metrics(mk(), "sq", w = 4), "odd")
  expect_equal(sum(info$arrangement), 10L)
})
