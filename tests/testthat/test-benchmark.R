z <- ss_simulate(70, 70, range = 8, seed = 1)
info <- ss_metric_list()

test_that("metrics of the window values only are invariant to spatial permutation", {
  b <- ss_benchmark(z, w = 9, n = 40, nperm = 39, seed = 1)
  s <- b$summary
  comp <- s$metric %in% info$metric[!info$arrangement]
  expect_true(all(s$permute_invariant[comp]))
  expect_true(all(!s$permute_invariant[!comp]))
  expect_true(all(s$category[comp] == "compositional"))
  expect_true(all(s$signif[comp] == 0))
  expect_true(all(s$signif[!comp] > 0))
  expect_equal(s$signif[s$metric == "sr1"], 1)
  expect_equal(length(b$cells), 40L)
})

test_that("lag-1 autocorrelation responds to permutation and distribution shape only", {
  b <- ss_benchmark(z, metrics = c("sr1", "sq"), w = 9, n = 40, nperm = 39,
                    seed = 1)
  s <- b$summary
  expect_true(s$shift_invariant[s$metric == "sr1"])
  expect_true(s$scale_invariant[s$metric == "sr1"])
  expect_false(s$reshape_invariant[s$metric == "sr1"])
  expect_identical(s$category, c("mixed", "compositional"))
  expect_true(s$shift_invariant[s$metric == "sq"])
  expect_false(s$scale_invariant[s$metric == "sq"])
})

test_that("the binning protocol decides the response of binned metrics", {
  m <- c("sent", "srao", "gcon", "slac")
  own <- ss_benchmark(z, metrics = m, w = 9, n = 40, nperm = 39,
                      tests = c("permute", "shift", "scale"),
                      protocol = "own", seed = 2)$summary
  fix <- ss_benchmark(z, metrics = m, w = 9, n = 40, nperm = 39,
                      tests = c("permute", "shift", "scale"),
                      protocol = "fixed", seed = 2)$summary
  expect_true(all(own$shift_invariant))
  expect_identical(own$scale_invariant, c(TRUE, FALSE, TRUE, TRUE))
  expect_false(any(fix$shift_invariant))
})

test_that("a user function and the built-in kernel agree", {
  f <- function(m) {
    d <- m - mean(m)
    dx <- (m[, -ncol(m)] - m[, -1])^2
    dy <- (m[-nrow(m), ] - m[-1, ])^2
    c(sq = sqrt(mean(d^2)), sr1 = 1 - mean(c(dx, dy)) / 2 / mean(d^2))
  }
  b1 <- ss_benchmark(z, metrics = c("sq", "sr1"), w = 9, n = 30, nperm = 39,
                     seed = 2)
  b2 <- ss_benchmark(z, fun = f, w = 9, nperm = 39, cells = b1$cells,
                     seed = 2)
  expect_identical(b1$summary$category, b2$summary$category)
  r1 <- b1$responses[b1$responses$test != "permute", ]
  r2 <- b2$responses[b2$responses$test != "permute", ]
  expect_equal(r1$response, r2$response, tolerance = 1e-8)
  b3 <- ss_benchmark(z, fun = function(m) c(sd = stats::sd(terra::values(m))),
                     w = 9, nperm = 39, cells = b1$cells, fun_input = "raster",
                     tests = c("permute", "shift"))
  expect_true(b3$summary$permute_invariant)
})

test_that("landscape level, explicit cells and missing values are handled", {
  bl <- ss_benchmark(z, metrics = c("sq", "sr1"), w = "landscape",
                     nperm = 39, seed = 1)
  expect_identical(bl$summary$category, c("compositional", "mixed"))
  expect_equal(bl$summary$configuration_response[1], 0)
  expect_true(is.na(bl$summary$configuration_response[2]))
  expect_error(ss_benchmark(z, metrics = "sq", w = 9, cells = 1L), "edge")
  zna <- z
  zna[35, 35] <- NA
  b <- ss_benchmark(zna, metrics = "sq", w = 9, n = 200, nperm = 39,
                    seed = 3)
  rc <- terra::rowColFromCell(zna, b$cells)
  near <- abs(rc[, 1] - 35) <= 4 & abs(rc[, 2] - 35) <= 4
  expect_false(any(near))
  expect_error(ss_benchmark(zna, metrics = "sq", w = "landscape"), "missing")
  expect_error(ss_benchmark(z, metrics = "sq", w = 9, nperm = 0), "nperm")
  expect_warning(ss_benchmark(z, metrics = "sq", w = 9, n = 5, nperm = 19,
                              tests = "permute"), "smallest attainable")
})

test_that("results are reproducible and printable", {
  a <- ss_benchmark(z, metrics = c("sq", "gcon"), w = 9, n = 20, nperm = 39,
                    tests = c("permute", "shift", "rearrange"), seed = 5)
  b <- ss_benchmark(z, metrics = c("sq", "gcon"), w = 9, n = 20, nperm = 39,
                    tests = c("permute", "shift", "rearrange"), seed = 5)
  expect_identical(a$responses, b$responses)
  expect_output(print(a), "ss_benchmark")
  expect_false(a$summary$rearrange_invariant[a$summary$metric == "sq"])
})

test_that("results on separate windows combine into the pooled result", {
  b <- ss_benchmark(z, metrics = c("sq", "sr1", "gcon"), w = 9, n = 30,
                    nperm = 39, seed = 7, details = TRUE)
  b1 <- ss_benchmark(z, metrics = c("sq", "sr1", "gcon"), w = 9, nperm = 39,
                     cells = b$cells[1:12], seed = 8, details = TRUE)
  b2 <- ss_benchmark(z, metrics = c("sq", "sr1", "gcon"), w = 9, nperm = 39,
                     cells = b$cells[13:30], seed = 9, details = TRUE)
  cb <- ss_benchmark_combine(b1, b2)
  expect_equal(cb$settings$n, 30L)
  expect_identical(cb$cells, b$cells)
  keep <- b$responses$test != "permute"
  expect_equal(cb$responses[keep, ], b$responses[keep, ])
  expect_identical(cb$summary$category, b$summary$category)
  expect_error(ss_benchmark_combine(b1, ss_benchmark(z, metrics = "sq", w = 9,
    nperm = 39, n = 5, details = TRUE)), "same metrics")
})

test_that("user functions receive the binning reference of the protocol", {
  q16 <- function(m, surface) {
    br <- stats::quantile(surface, seq(0, 1, length.out = 17), names = FALSE)
    g <- matrix(findInterval(m, br[2:16]), nrow(m))
    c(levels_mean = mean(g), levels_lag = mean(abs(g[, -1] - g[, -ncol(g)])))
  }
  own <- ss_benchmark(z, fun = q16, w = 9, n = 20, nperm = 39,
                      protocol = "own", seed = 4)$summary
  expect_true(all(own$shift_invariant & own$scale_invariant &
                    own$reshape_invariant))
  expect_identical(own$category, c("unresponsive", "configurational"))
  fix <- ss_benchmark(z, fun = q16, w = 9, n = 20, nperm = 39,
                      protocol = "fixed", seed = 4)$summary
  expect_false(any(fix$shift_invariant))
})

test_that("rounding in a near-zero landscape mean is not taken as a response", {
  z0 <- ss_simulate(40, 40, range = 5, seed = 3)
  b <- ss_benchmark(z0, metrics = c("smean", "sq"), w = "landscape", nperm = 39,
                    tests = c("permute", "shift", "scale", "reshape"), seed = 1)
  s <- b$summary
  expect_true(all(s$permute_invariant))
  expect_identical(s$category, c("compositional", "compositional"))
  expect_equal(s$signif, c(0, 0))
  expect_true(is.na(s$composition_response[s$metric == "smean"]))
})
