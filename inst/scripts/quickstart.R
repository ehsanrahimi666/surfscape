library(terra)
library(surfscape)

## ---- 1. a simulated landscape with known spatial structure ------------------
z <- ss_simulate(120, 120, range = 10, seed = 1)

## ---- 2. classify the built-in metrics ---------------------------------------
b <- ss_benchmark(z, w = 15, n = 100, nperm = 39,
                  tests = c("permute", "shift", "scale", "reshape"), seed = 1)
b
plot(b)

## ---- 3. a metric from another package, on the same windows ------------------
moran <- function(r) c(moran = unname(terra::autocor(r, method = "moran", global = TRUE)))
ss_benchmark(z, fun = moran, fun_input = "raster", w = 15, nperm = 39,
             tests = c("permute", "shift", "scale", "reshape"),
             cells = b$cells, seed = 1)

## ---- 4. decompose a change between two dates --------------------------------
z1 <- rast(system.file("extdata", "ndvi_2013.tif", package = "surfscape"))
z2 <- rast(system.file("extdata", "ndvi_2023.tif", package = "surfscape"))
d <- ss_decompose(z1, z2, metrics = c("sq", "sr1"), w = 15)
plot(d[[c("sr1_total", "sr1_global", "sr1_local", "sr1_config")]])
ss_decompose_landscape(z1, z2, metrics = c("sq", "sent", "sr1", "gcon"))
