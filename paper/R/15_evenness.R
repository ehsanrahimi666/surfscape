## How the room for configuration depends on composition: across the 100 NDVI
## windows, the spread of a metric among random arrangements of the same
## values, against the evenness of the values in the window.

suppressPackageStartupMessages({
  library(here)
  library(surfscape)
})

files <- list.files(here("output", "rds", "units"), pattern = "^ndvi_chunk.*\\.rds$",
                    full.names = TRUE)
bl <- lapply(files, function(f) readRDS(f)$results$builtin_fixed)
cb <- ss_benchmark_combine(bl)
obs <- cb$values$obs
perm <- cb$values$perm
evenness <- obs[, "sent"] / log(32)
rows <- lapply(c("gcon", "gcor", "sr1", "slac"), function(m) {
  P <- matrix(perm[, , m], nrow = dim(perm)[1])
  sd_perm <- apply(P, 1, stats::sd, na.rm = TRUE)
  ses <- (obs[, m] - rowMeans(P)) / sd_perm
  data.frame(metric = m, n = nrow(obs),
             n_undefined = sum(!is.finite(obs[, m])),
             n_complete = sum(is.finite(sd_perm)),
             rho_evenness_sdperm = stats::cor(evenness, sd_perm, method = "spearman",
                                              use = "complete.obs"),
             rho_evenness_abs_ses = stats::cor(evenness, abs(ses), method = "spearman",
                                               use = "complete.obs"),
             rho_sq_sdperm = stats::cor(obs[, "sq"], sd_perm, method = "spearman",
                                        use = "complete.obs"),
             rho_evenness_sq = stats::cor(evenness, obs[, "sq"], method = "spearman"))
})
res <- do.call(rbind, rows)
write.csv(res, here("output", "tables", "evenness_configuration_room.csv"), row.names = FALSE)
write.csv(data.frame(evenness = evenness, sq = obs[, "sq"],
                     sd_perm_gcon = apply(matrix(perm[, , "gcon"], nrow = dim(perm)[1]), 1,
                                          stats::sd, na.rm = TRUE)),
          here("output", "tables", "evenness_windows_ndvi.csv"), row.names = FALSE)
print(res)
