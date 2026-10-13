## geodiv computes the areal autocorrelation function from the Fourier
## transform of the window values without removing their mean. This analysis
## repeats the benchmark of the autocorrelation and spectral metrics on windows
## from which the window mean has been subtracted, on the same 100 NDVI windows.

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(surfscape)
})
source(here("R", "external_metrics.R"))

z <- rast(here("data", "derived", "ndvi2013_aw30.tif"))
z23 <- rast(here("data", "derived", "ndvi2023_aw30.tif"))
ndw <- read.csv(here("output", "tables", "ndvi_windows_w15.csv"))

spectral_metrics <- function(r, centre) {
  if (centre) r <- r - mean(values(r)[, 1])
  scl <- .guard("geodiv_scl", unlist(geodiv::scl(r)), n = 4L)
  stx <- .guard("geodiv_stxr", unlist(geodiv::stxr(r)), n = 2L)
  srw <- .guard("geodiv_srw", unlist(geodiv::srw(r)), n = 3L)
  c(geodiv_scl20 = as.numeric(scl)[1], geodiv_str20 = as.numeric(stx)[1],
    geodiv_srw = as.numeric(srw)[1], geodiv_srwi = as.numeric(srw)[2],
    geodiv_shw = as.numeric(srw)[3])
}

ext_log_reset()
res <- lapply(c(raw = FALSE, centred = TRUE), function(centre) {
  ss_benchmark(z, fun = function(r) spectral_metrics(r, centre), fun_input = "raster",
               w = 15L, nperm = 39L,
               tests = c("permute", "shift", "scale", "reshape", "target", "rearrange"),
               skew = 0, target = z23, arrangement = z23, cells = ndw$cell,
               seed = 800000L, details = TRUE)
})
tab <- do.call(rbind, lapply(names(res), function(k)
  cbind(variant = k, res[[k]]$summary)))
write.csv(tab, here("output", "tables", "benchmark_geodiv_centred_summary.csv"), row.names = FALSE)
rsp <- do.call(rbind, lapply(names(res), function(k)
  cbind(variant = k, res[[k]]$responses)))
write.csv(rsp, here("output", "tables", "benchmark_geodiv_centred_responses.csv"),
          row.names = FALSE)
write.csv(ext_log(), here("logs", "benchmark_geodiv_centred_external_log.csv"), row.names = FALSE)
print(tab[, c("variant", "metric", "category", "configuration_response",
              "composition_response", "shift_invariant", "scale_invariant")])
