## Figures for the README, from the example data shipped with the package.
## Run from the package root: Rscript inst/scripts/make_figures.R

library(terra)
library(surfscape)

out <- file.path("man", "figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
div <- hcl.colors(100, "Blue-Red 3")
seq_col <- hcl.colors(100, "YlGn", rev = TRUE)

z1 <- rast(system.file("extdata", "ndvi_2013.tif", package = "surfscape"))
z2 <- rast(system.file("extdata", "ndvi_2023.tif", package = "surfscape"))

## 1. the tests on a simulated landscape
z <- ss_simulate(120, 120, range = 10, seed = 1)
png(file.path(out, "README-tests.png"), width = 1800, height = 480, res = 150)
par(mfrow = c(1, 5), mar = c(0.5, 0.5, 2, 0.5))
panels <- list(original = z,
               shift = ss_perturb(z, "shift", 0.5),
               scale = ss_perturb(z, "scale", 1.5),
               reshape = ss_perturb(z, "reshape", 2),
               rearrange = ss_perturb(z, "rearrange", sim_range = 10, seed = 2))
rng <- range(unlist(lapply(panels, function(r) minmax(r))))
for (nm in names(panels)) {
  plot(panels[[nm]], axes = FALSE, legend = FALSE, main = nm, range = rng,
       col = hcl.colors(100, "viridis"), mar = NA)
}
dev.off()

## 2. benchmark of the built-in metrics on the 2013 NDVI example
b <- ss_benchmark(z1, w = 15, n = 100, nperm = 39,
                  tests = c("permute", "shift", "scale", "reshape", "target"),
                  target = z2, skew = 0, seed = 1)
png(file.path(out, "README-benchmark.png"), width = 1400, height = 1200, res = 170)
par(mar = c(4.5, 4.5, 1, 1))
plot(b)
dev.off()

## 3. decomposition of the 2013-2023 change in lag-1 autocorrelation
d <- ss_decompose(z1, z2, metrics = c("sq", "sr1"), w = 15)
png(file.path(out, "README-decomposition.png"), width = 2000, height = 520, res = 150)
par(mfrow = c(1, 5), mar = c(0.5, 0.5, 2, 3.5))
plot(z1, axes = FALSE, main = "NDVI 2013", col = seq_col, mar = NA)
labs <- c(sr1_total = "total change", sr1_global = "global composition",
          sr1_local = "local composition", sr1_config = "configuration")
mx <- max(abs(quantile(values(d[[names(labs)]]), c(0.02, 0.98), na.rm = TRUE)))
for (k in names(labs)) {
  plot(clamp(d[[k]], -mx, mx, values = TRUE), axes = FALSE, main = labs[[k]],
       range = c(-mx, mx) * 1.001, col = div, mar = NA)
}
dev.off()

## 4. share of configuration across window sizes
sp <- ss_scale_profile(z1, z2, metrics = c("sq", "sr1", "gcon", "slac"),
                       w_seq = c(5, 9, 15, 31), progress = FALSE)
png(file.path(out, "README-scale.png"), width = 1200, height = 900, res = 170)
par(mar = c(4.5, 4.5, 1, 1))
cols <- c(sq = "grey40", sr1 = "#d95f02", gcon = "#7570b3", slac = "#1b9e77")
plot(NA, xlim = range(sp$w), ylim = c(0, 1), log = "x",
     xlab = "Window size (cells)", ylab = "Share of change due to configuration")
for (m in names(cols)) {
  s <- sp[sp$metric == m, ]
  lines(s$w, s$share_config, type = "b", pch = 19, col = cols[[m]])
}
legend("topright", legend = names(cols), col = cols, pch = 19, lty = 1, bty = "n")
dev.off()

print(b$summary[, c("metric", "category", "configuration_response",
                    "composition_response", "signif")])
print(sp[, c("w", "metric", "share_config")])
