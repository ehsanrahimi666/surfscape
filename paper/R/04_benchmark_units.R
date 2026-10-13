## Benchmark at 15 x 15 cell windows: 100 windows of the 2013 NDVI surface
## (ten chunks of ten) and 100 windows for each simulated range (ten
## landscapes of ten windows). Run as: Rscript R/04_benchmark_units.R <worker> <nworkers>

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
ndw <- read.csv(here("output", "tables", "ndvi_windows_w15.csv"))
units <- rbind(
  data.frame(unit = sprintf("ndvi_chunk%02d", 1:10), type = "ndvi", chunk = 1:10,
             sim_row = NA_integer_, seed = 400000L + 1:10),
  data.frame(unit = paste0("sim_", sim$id), type = "virtual", chunk = NA_integer_,
             sim_row = seq_len(nrow(sim)),
             seed = 300000L + 100L * sim$range + sim$rep))
mine <- units[seq(worker, nrow(units), by = nworkers), ]

w <- 15L
h <- (w - 1L) %/% 2L
all_metrics <- ss_metric_list()$metric

for (i in seq_len(nrow(mine))) {
  u <- mine[i, ]
  f_out <- here("output", "rds", "units", paste0(u$unit, ".rds"))
  if (file.exists(f_out)) next
  if (u$type == "ndvi") {
    z <- rast(here("data", "derived", "ndvi2013_aw30.tif"))
    z23 <- rast(here("data", "derived", "ndvi2023_aw30.tif"))
    cells <- ndw$cell[ndw$chunk == u$chunk]
    tests <- c("permute", "shift", "scale", "reshape", "target", "rearrange")
    skew <- 0
    target <- z23
    arrangement <- z23
  } else {
    s <- sim[u$sim_row, ]
    z <- rast(here(s$file_landscape))
    arrangement <- rast(here(s$file_arrangement))
    set.seed(u$seed)
    rows <- (h + 1L):(nrow(z) - h)
    cols <- (h + 1L):(ncol(z) - h)
    cand <- as.vector(outer((rows - 1L) * ncol(z), cols, "+"))
    cells <- sort(sample(cand, 10L))
    tests <- c("permute", "shift", "scale", "reshape", "rearrange")
    skew <- 2
    target <- NULL
  }
  run <- function(...) {
    t0 <- proc.time()[3]
    b <- ss_benchmark(z, w = w, nperm = 39L, tests = tests, shift = 0.5, scale = 1.5,
                      skew = skew, target = target, arrangement = arrangement,
                      cells = cells, seed = u$seed, details = TRUE, ...)
    b$elapsed <- unname(proc.time()[3] - t0)
    b
  }
  ext_log_reset()
  res <- list(
    builtin_fixed = run(metrics = all_metrics, protocol = "fixed"),
    builtin_own = run(metrics = all_metrics, protocol = "own"),
    texture_fixed = run(fun = texture_metrics, protocol = "fixed"),
    texture_own = run(fun = texture_metrics, protocol = "own"),
    rasterdiv_fixed = run(fun = rasterdiv_metrics, protocol = "fixed"),
    rasterdiv_own = run(fun = rasterdiv_metrics, protocol = "own"),
    unbinned = run(fun = unbinned_metrics, fun_input = "raster", protocol = "fixed"))
  saveRDS(list(unit = u, cells = cells, results = res, log = ext_log()), f_out)
  message(format(Sys.time()), " ", u$unit, " done in ",
          round(sum(vapply(res, `[[`, numeric(1), "elapsed"))), " s")
}
