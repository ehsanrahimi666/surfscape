## Controlled test of the decomposition on simulated landscapes: the second
## date is built from the first with a known change in composition (kappa)
## and a known change in arrangement (lambda), in a factorial design.

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(surfscape)
})

sim <- read.csv(here("output", "tables", "simulated_landscapes.csv"))
sim <- sim[sim$range == 10, ]
kappa <- c(0, 1, 2)
lambda <- c(0, 0.25, 0.5, 0.75, 1)
w <- 15L
metrics <- ss_metric_list()$metric
psd <- function(v) sqrt(mean((v - mean(v))^2))

## composition change of strength k: the mean moves by 0.5 k SD, the SD is
## multiplied by 1 + 0.5 k and the skewness becomes k; k = 0 keeps the values
new_values <- function(v, k) {
  if (k == 0) return(v)
  u <- (rank(v, ties.method = "average") - 0.5) / length(v)
  q <- stats::qgamma(u, shape = 4 / k^2)
  q <- (q - mean(q)) / psd(q)
  mean(v) + 0.5 * k * psd(v) + (1 + 0.5 * k) * psd(v) * q
}

rows <- list()
lrows <- list()
for (i in seq_len(nrow(sim))) {
  z1 <- rast(here(sim$file_landscape[i]))
  y <- rast(here(sim$file_arrangement[i]))
  v1 <- values(z1)[, 1]
  for (k in kappa) {
    v2 <- new_values(v1, k)
    for (l in lambda) {
      a <- rast(z1, vals = (1 - l) * v1 + l * values(y)[, 1])
      z2 <- ss_rank_transfer(a, rast(z1, vals = v2))
      d <- ss_decompose(z1, z2, metrics = metrics, w = w)
      dv <- values(d)
      ok <- stats::complete.cases(dv)
      for (m in metrics) {
        cm <- dv[ok, paste0(m, c("_total", "_global", "_local", "_config"))]
        rows[[length(rows) + 1L]] <- data.frame(
          landscape = sim$id[i], kappa = k, lambda = l, metric = m, n = nrow(cm),
          rank_cor = stats::cor(v1, values(z2)[, 1], method = "spearman"),
          mad_total = mean(abs(cm[, 1])), mad_global = mean(abs(cm[, 2])),
          mad_local = mean(abs(cm[, 3])), mad_config = mean(abs(cm[, 4])),
          mean_global = mean(cm[, 2]), mean_local = mean(cm[, 3]),
          mean_config = mean(cm[, 4]),
          max_abs_global = max(abs(cm[, 2])), max_abs_local = max(abs(cm[, 3])),
          max_abs_config = max(abs(cm[, 4])),
          max_residual = max(abs(cm[, 1] - rowSums(cm[, 2:4]))))
      }
      dl <- ss_decompose_landscape(z1, z2, metrics = metrics)
      lrows[[length(lrows) + 1L]] <- cbind(landscape = sim$id[i], kappa = k, lambda = l, dl)
    }
  }
  message("landscape ", sim$id[i], " done")
}
res <- do.call(rbind, rows)
lres <- do.call(rbind, lrows)
write.csv(res, here("output", "tables", "decomposition_virtual_windows.csv"), row.names = FALSE)
write.csv(lres, here("output", "tables", "decomposition_virtual_landscape.csv"), row.names = FALSE)
