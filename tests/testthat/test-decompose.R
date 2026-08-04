mk <- function(nr = 45, nc = 45, seed = 1, smooth = 5) {
  set.seed(seed)
  r <- terra::rast(nrows = nr, ncols = nc, xmin = 0, xmax = nc, ymin = 0, ymax = nr)
  terra::values(r) <- stats::rnorm(terra::ncell(r))
  terra::focal(r, smooth, "mean", na.rm = TRUE)
}
ALL <- ss_metric_list()
COMP <- ALL$metric[!ALL$arrangement]
CONF <- ALL$metric[ALL$arrangement]

test_that("the three components sum exactly to the total change", {
  z1 <- mk(seed = 1, smooth = 5); z2 <- mk(seed = 2, smooth = 9)
  d <- ss_decompose(z1, z2, metrics = c("sq", "sent", "sr1", "gcon", "slac"), w = 9)
  for (m in c("sq", "sent", "sr1", "gcon", "slac")) {
    v <- terra::values(d[[paste0(m, c("_total", "_global", "_local", "_config"))]])
    v <- v[stats::complete.cases(v), , drop = FALSE]
    expect_gt(nrow(v), 100)
    expect_equal(v[, 1], v[, 2] + v[, 3] + v[, 4], tolerance = 1e-10, info = m)
  }
})

test_that("a monotone value change is attributed entirely to the global component", {
  z1 <- mk(seed = 3)
  z2 <- z1; terra::values(z2) <- terra::values(z1)[, 1]^3 * 5 + 2   # strictly monotone
  d <- ss_decompose(z1, z2, metrics = c("sq", "sr1", "gcon"), w = 9)
  for (m in c("sq", "sr1", "gcon")) {
    v <- terra::values(d[[paste0(m, c("_total", "_global", "_local", "_config"))]])
    v <- v[stats::complete.cases(v), , drop = FALSE]
    expect_equal(max(abs(v[, 4])), 0, tolerance = 1e-9, info = paste(m, "config"))
    expect_equal(max(abs(v[, 3])), 0, tolerance = 1e-9, info = paste(m, "local"))
    expect_equal(v[, 2], v[, 1], tolerance = 1e-9, info = paste(m, "global"))
  }
})

test_that("a pure spatial rearrangement has a zero global component", {
  z1 <- mk(seed = 4, smooth = 5)
  z2 <- ss_rank_transfer(mk(seed = 5, smooth = 11), z1)  # values of z1, new arrangement
  d <- ss_decompose(z1, z2, metrics = c("sq", "sr1", "gcon"), w = 9)
  for (m in c("sq", "sr1", "gcon")) {
    v <- terra::values(d[[paste0(m, "_global")]])
    v <- v[is.finite(v)]
    expect_equal(max(abs(v)), 0, tolerance = 1e-9, info = m)
  }
})

test_that("configurational change is identically zero for multiset metrics", {
  z1 <- mk(seed = 6, smooth = 5); z2 <- mk(seed = 7, smooth = 13)
  d <- ss_decompose(z1, z2, metrics = COMP, w = 9, components = "config")
  v <- terra::values(d)
  v <- v[stats::complete.cases(v), , drop = FALSE]
  expect_gt(nrow(v), 100)
  mx <- apply(abs(v), 2, max)
  # multiset metrics are invariant to rearrangement in exact arithmetic; in
  # floating point most are bit-exact, s10z is limited by summation order and
  # scv by a window mean that can approach zero
  expect_lt(max(mx[paste0(setdiff(COMP, "scv"), "_config")]), 1e-12)
  expect_lt(unname(mx["scv_config"]), 1e-6)
})

test_that("configurational change is non-zero for arrangement metrics", {
  z1 <- mk(seed = 8, smooth = 5); z2 <- mk(seed = 9, smooth = 13)
  d <- ss_decompose(z1, z2, metrics = CONF, w = 9, components = "config")
  for (m in CONF) {
    v <- terra::values(d[[paste0(m, "_config")]]); v <- v[is.finite(v)]
    expect_gt(mean(abs(v)), 0, label = m)
  }
})

test_that("the decomposition is invariant to block size", {
  z1 <- mk(nr = 51, nc = 33, seed = 10); z2 <- mk(nr = 51, nc = 33, seed = 11, smooth = 9)
  a <- ss_decompose(z1, z2, metrics = c("sq", "sr1"), w = 7)
  b <- ss_decompose(z1, z2, metrics = c("sq", "sr1"), w = 7, chunk_rows = 4)
  expect_identical(terra::values(a), terra::values(b))
})

test_that("the sensitivity audit separates what can from what cannot see arrangement", {
  r <- mk(seed = 12, smooth = 7)
  s <- ss_sensitivity(r, w = 9, n = 150, nperm = 25, seed = 1)
  z <- s[!s$arrangement, ]
  a <- s[s$arrangement, ]
  expect_equal(max(z$csi), 0, tolerance = 1e-12)
  expect_equal(max(z$effect), 0, tolerance = 1e-12)
  expect_true(all(a$csi > 0))
  expect_true(all(a$effect > 0))
  expect_true(all(s$csi >= 0 & s$csi <= 1))
})

test_that("scale profile and point extraction return tidy output", {
  z1 <- mk(seed = 13); z2 <- mk(seed = 14, smooth = 11)
  p <- ss_scale_profile(z1, z2, metrics = c("sq", "sr1"), w_seq = c(5, 9),
                        progress = FALSE)
  expect_equal(nrow(p), 4L)
  expect_true(all(c("share_config", "sign_disagree") %in% names(p)))
  expect_equal(p$mad_config[p$metric == "sq"], c(0, 0), tolerance = 1e-12)

  d <- ss_decompose(z1, z2, metrics = "sr1", w = 5)
  pts <- data.frame(x = c(10, 25, 35), y = c(12, 22, 30))
  e <- ss_extract(d, pts, radii = c(3, 6))
  expect_true(all(c("site", "radius", "layer", "value") %in% names(e)))
  expect_equal(length(unique(e$radius)), 2L)
  expect_equal(length(unique(e$site)), 3L)
})
