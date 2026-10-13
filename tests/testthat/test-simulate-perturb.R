psd <- function(v) sqrt(mean((v - mean(v))^2))

test_that("simulated landscapes are standardised and have the requested structure", {
  z <- ss_simulate(150, 150, range = 12, seed = 3)
  expect_s4_class(z, "SpatRaster")
  expect_equal(dim(z)[1:2], c(150, 150))
  v <- terra::values(z)[, 1]
  expect_equal(mean(v), 0, tolerance = 1e-12)
  expect_equal(psd(v), 1, tolerance = 1e-12)
  m <- terra::as.matrix(z, wide = TRUE)
  r1 <- stats::cor(as.vector(m[, -1]), as.vector(m[, -ncol(m)]))
  expect_equal(r1, exp(-3 / 12), tolerance = 0.05)
  z30 <- ss_simulate(150, 150, range = 30, seed = 3)
  m30 <- terra::as.matrix(z30, wide = TRUE)
  expect_gt(stats::cor(as.vector(m30[, -1]), as.vector(m30[, -ncol(m30)])), r1)
})

test_that("simulation is reproducible and leaves the session's random numbers alone", {
  set.seed(10); a <- stats::runif(1)
  set.seed(10); z1 <- ss_simulate(30, 30, 5, seed = 1); b <- stats::runif(1)
  expect_identical(a, b)
  z2 <- ss_simulate(30, 30, 5, seed = 1)
  expect_identical(terra::values(z1), terra::values(z2))
})

test_that("supplied values change the distribution but not the arrangement", {
  z <- ss_simulate(40, 40, 6, seed = 2)
  set.seed(1); vals <- stats::rexp(1600)
  zg <- ss_simulate(40, 40, 6, seed = 2, values = vals)
  expect_equal(sort(terra::values(zg)[, 1]), sort(vals))
  expect_equal(rank(terra::values(zg)[, 1]), rank(terra::values(z)[, 1]))
})

test_that("composition changes keep the spatial ordering of values", {
  z <- ss_simulate(50, 50, 8, seed = 4)
  v <- terra::values(z)[, 1]
  sh <- terra::values(ss_perturb(z, "shift", 0.5))[, 1]
  expect_equal(sh - v, rep(0.5 * psd(v), length(v)))
  sc <- terra::values(ss_perturb(z, "scale", 1.5))[, 1]
  expect_equal(mean(sc), mean(v))
  expect_equal(psd(sc), 1.5 * psd(v))
  rs <- terra::values(ss_perturb(z, "reshape", 2))[, 1]
  expect_equal(rank(rs), rank(v))
  expect_equal(mean(rs), mean(v))
  expect_equal(psd(rs), psd(v))
  expect_gt(mean((rs - mean(rs))^3) / psd(rs)^3, 1.5)
  set.seed(2); tg <- terra::rast(z, vals = stats::rexp(2500))
  tv <- terra::values(ss_perturb(z, "target", target = tg))[, 1]
  expect_equal(sort(tv), sort(terra::values(tg)[, 1]))
  expect_equal(rank(tv), rank(v))
})

test_that("reshape keeps tied values tied", {
  z <- terra::rast(nrows = 10, ncols = 10)
  terra::values(z) <- rep(1:20, 5)
  rs <- terra::values(ss_perturb(z, "reshape", 2))[, 1]
  expect_equal(length(unique(rs)), 20L)
  expect_equal(rank(rs), rank(terra::values(z)[, 1]))
})

test_that("configuration changes keep the values", {
  z <- ss_simulate(50, 50, 8, seed = 4)
  v <- terra::values(z)[, 1]
  pm <- terra::values(ss_perturb(z, "permute", seed = 1))[, 1]
  expect_equal(sort(pm), sort(v))
  ar <- terra::rast(z, vals = terra::values(ss_simulate(50, 50, 20, seed = 9)))
  ra <- terra::values(ss_perturb(z, "rearrange", arrangement = ar))[, 1]
  expect_equal(sort(ra), sort(v))
  expect_equal(rank(ra), rank(terra::values(ar)[, 1]))
  rb <- ss_perturb(z, "rearrange", sim_range = 20, seed = 3)
  expect_equal(sort(terra::values(rb)[, 1]), sort(v))
  expect_error(ss_perturb(z, "rearrange"), "sim_range")
  expect_error(ss_perturb(z, "target"), "target")
})
