## Extent and grain analysis. Nested landscape windows of 5, 9, 15, 31 and 61
## cells around common centres, and the whole landscape as one window, for the
## simulated landscapes and the 2013 NDVI surface; and the NDVI surface
## aggregated to 60 m and 120 m with windows of 7 and 15 cells.
## Run as: Rscript R/05_extent_grain_units.R <worker> <nworkers>

args <- as.integer(commandArgs(trailingOnly = TRUE))
worker <- args[1]
nworkers <- args[2]

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(surfscape)
})
source(here("R", "external_metrics.R"))

sim <- read.csv(here("output", "tables", "simulated_landscapes.csv"))
extents <- c(5L, 9L, 15L, 31L, 61L)
all_metrics <- ss_metric_list()$metric

units <- rbind(
  data.frame(unit = paste0("ext_sim_", sim$id), type = "sim", sim_row = seq_len(nrow(sim)),
             chunk = NA_integer_, grain = 30L,
             seed = 500000L + 100L * sim$range + sim$rep),
  data.frame(unit = sprintf("ext_ndvi_chunk%02d", 1:10), type = "ndvi", sim_row = NA_integer_,
             chunk = 1:10, grain = 30L, seed = 600000L + 1:10),
  data.frame(unit = "ext_ndvi_landscape", type = "ndvi_landscape", sim_row = NA_integer_,
             chunk = NA_integer_, grain = 30L, seed = 610000L),
  data.frame(unit = sprintf("grain%03d_chunk%02d", rep(c(60L, 120L), each = 2), rep(1:2, 2)),
             type = "grain", sim_row = NA_integer_, chunk = rep(1:2, 2),
             grain = rep(c(60L, 120L), each = 2), seed = 700000L + 1:4))
mine <- units[seq(worker, nrow(units), by = nworkers), ]

centres <- function(z, w, n, seed) {
  set.seed(seed)
  h <- (w - 1L) %/% 2L
  rows <- (h + 1L):(nrow(z) - h)
  cols <- (h + 1L):(ncol(z) - h)
  sort(sample(as.vector(outer((rows - 1L) * ncol(z), cols, "+")), n))
}

bench_sets <- function(z, w, cells, tests, skew, target, arrangement, seed) {
  run <- function(...) {
    t0 <- proc.time()[3]
    b <- ss_benchmark(z, w = w, nperm = 39L, tests = tests, shift = 0.5, scale = 1.5,
                      skew = skew, target = target, arrangement = arrangement,
                      cells = cells, seed = seed, details = TRUE, ...)
    b$elapsed <- unname(proc.time()[3] - t0)
    b
  }
  list(builtin_fixed = run(metrics = all_metrics, protocol = "fixed"),
       glcmtex_fixed = run(fun = glcm_texture_metrics, protocol = "fixed"),
       glcmtex_own = run(fun = glcm_texture_metrics, protocol = "own"),
       terra = run(fun = terra_metrics, fun_input = "raster", protocol = "fixed"))
}

for (i in seq_len(nrow(mine))) {
  u <- mine[i, ]
  f_out <- here("output", "rds", "units", paste0(u$unit, ".rds"))
  if (file.exists(f_out)) next
  ext_log_reset()
  if (u$type == "sim") {
    s <- sim[u$sim_row, ]
    z <- rast(here(s$file_landscape))
    arrangement <- rast(here(s$file_arrangement))
    target <- NULL
    tests <- c("permute", "shift", "scale", "reshape", "rearrange")
    skew <- 2
  } else {
    sfx <- if (u$grain == 30L) "aw30" else sprintf("aw%dm", u$grain)
    z <- rast(here("data", "derived", sprintf("ndvi2013_%s.tif", sfx)))
    target <- rast(here("data", "derived", sprintf("ndvi2023_%s.tif", sfx)))
    arrangement <- target
    tests <- c("permute", "shift", "scale", "reshape", "target", "rearrange")
    skew <- 0
  }
  res <- list()
  if (u$type == "sim") {
    cells <- centres(z, max(extents), 10L, u$seed)
    for (w in extents) res[[paste0("w", w)]] <- bench_sets(z, w, cells, tests, skew, target,
                                                          arrangement, u$seed)
    res$landscape <- bench_sets(z, "landscape", NULL, tests, skew, target, arrangement, u$seed)
  } else if (u$type == "ndvi") {
    cells <- centres(z, max(extents), 100L, 610001L)[(u$chunk - 1L) * 10L + 1:10]
    for (w in extents) res[[paste0("w", w)]] <- bench_sets(z, w, cells, tests, skew, target,
                                                          arrangement, u$seed)
  } else if (u$type == "ndvi_landscape") {
    cells <- NULL
    res$landscape <- bench_sets(z, "landscape", NULL, tests, skew, target, arrangement, u$seed)
  } else {
    cells <- centres(z, 15L, 100L, 710000L + u$grain)[(u$chunk - 1L) * 50L + 1:50]
    for (w in c(7L, 15L)) res[[paste0("w", w)]] <- bench_sets(z, w, cells, tests, skew, target,
                                                             arrangement, u$seed)
  }
  saveRDS(list(unit = u, cells = cells, results = res, log = ext_log()), f_out)
  message(format(Sys.time()), " ", u$unit, " done")
}
