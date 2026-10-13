## Redundancy among metrics and selection of a key subset.

suppressPackageStartupMessages({
  library(here)
})

obs <- readRDS(here("output", "rds", "benchmark_main_observed.rds"))
summ_all <- read.csv(here("output", "tables", "benchmark_main_summary.csv"))
summ <- summ_all[summ_all$protocol == "fixed", ]

usable <- function(m) {
  ok <- apply(m, 2, function(v) sum(is.finite(v)) >= 0.9 * length(v) &&
                stats::sd(v, na.rm = TRUE) > 0)
  m[, ok, drop = FALSE]
}
sp_cor <- function(m) stats::cor(m, method = "spearman", use = "pairwise.complete.obs")

ndvi <- usable(obs$ndvi)
sim <- usable(rbind(obs$r03, obs$r10, obs$r30))
common <- intersect(colnames(ndvi), colnames(sim))
excluded <- setdiff(union(colnames(obs$ndvi), colnames(obs$r10)), common)
write.csv(data.frame(metric = excluded), here("output", "tables", "redundancy_excluded.csv"),
          row.names = FALSE)
c_ndvi <- sp_cor(ndvi[, common])
c_sim <- sp_cor(sim[, common])
write.csv(c_ndvi, here("output", "tables", "redundancy_spearman_ndvi.csv"))
write.csv(c_sim, here("output", "tables", "redundancy_spearman_sim.csv"))

hc <- stats::hclust(stats::as.dist(1 - abs(c_ndvi)), method = "average")
cl <- stats::cutree(hc, h = 0.3)
saveRDS(hc, here("output", "rds", "redundancy_hclust_ndvi.rds"))

cat_of <- function(m, L, protocol = "fixed") {
  x <- summ_all$category[summ_all$metric == m & summ_all$landscape == L &
                           summ_all$protocol == protocol]
  if (length(x)) x[1] else NA_character_
}
inv_of <- function(m, test) {
  x <- summ[summ$metric == m & summ$landscape == "ndvi", paste0(test, "_invariant")]
  if (length(x)) x[1] else NA
}
groups <- data.frame(metric = common, cluster = cl[common],
                     category_ndvi = vapply(common, cat_of, "", L = "ndvi"),
                     category_r10 = vapply(common, cat_of, "", L = "r10"),
                     shift_invariant = vapply(common, inv_of, NA, test = "shift"),
                     scale_invariant = vapply(common, inv_of, NA, test = "scale"),
                     reshape_invariant = vapply(common, inv_of, NA, test = "reshape"))
write.csv(groups[order(groups$cluster, groups$metric), ],
          here("output", "tables", "redundancy_clusters.csv"), row.names = FALSE)

## preference lists by role (rules ii and iv of Section 2.7: implementations
## computable both in moving windows and for a whole landscape, rasterdiv and
## glcm therefore excluded). Within each role the first candidate that is not
## redundant (|rho| >= 0.7) with an already selected metric in both the NDVI
## and the simulated windows is selected; redundancy in only one of them is
## reported. The outcome of the strict rule, |rho| < 0.7 in both, is recorded.
roles <- list(
  location = c("surfscape_smean"),
  spread = c("surfscape_sq", "surfscape_sa", "surfscape_sz"),
  shape = c("surfscape_ssk", "surfscape_sku"),
  diversity = c("surfscape_sent", "surfscape_srao"),
  configuration_autocorrelation = c("terra_moran", "surfscape_sr1", "terra_geary",
                                    "surfscape_gcor", "glcmtex_range_correlation"),
  configuration_texture = c("glcmtex_prob_contrast", "glcmtex_prob_dissimilarity",
                            "glcmtex_prob_homogeneity", "surfscape_gcon",
                            "glcmtex_range_contrast", "surfscape_slac", "geodiv_sds"))
redundant_in <- function(cand, chosen, cm) {
  if (!length(chosen)) return(character())
  chosen[abs(cm[cand, chosen]) >= 0.7]
}
sel <- character()
sel_strict <- character()
log_rows <- list()
for (r in names(roles)) {
  picked <- NA_character_
  for (cand in roles[[r]]) {
    if (!cand %in% common) next
    rn <- redundant_in(cand, sel, c_ndvi)
    rs <- redundant_in(cand, sel, c_sim)
    both <- intersect(rn, rs)
    log_rows[[length(log_rows) + 1L]] <- data.frame(
      role = r, candidate = cand,
      max_abs_rho_ndvi = if (length(sel)) max(abs(c_ndvi[cand, sel])) else 0,
      max_abs_rho_sim = if (length(sel)) max(abs(c_sim[cand, sel])) else 0,
      redundant_ndvi = paste(rn, collapse = ";"), redundant_sim = paste(rs, collapse = ";"),
      selected = !length(both))
    if (!length(both)) {
      picked <- cand
      break
    }
  }
  if (!is.na(picked)) sel <- c(sel, stats::setNames(picked, r))
  for (cand in roles[[r]]) {
    if (!cand %in% common) next
    ok <- !length(redundant_in(cand, sel_strict, c_ndvi)) &&
      !length(redundant_in(cand, sel_strict, c_sim))
    if (ok) {
      sel_strict <- c(sel_strict, stats::setNames(cand, r))
      break
    }
  }
}
write.csv(data.frame(role = names(sel_strict), metric = unname(sel_strict)),
          here("output", "tables", "subset_strict_rule.csv"), row.names = FALSE)
write.csv(do.call(rbind, log_rows), here("output", "tables", "subset_selection_log.csv"),
          row.names = FALSE)
subset_tab <- data.frame(role = names(sel), metric = unname(sel),
                         category_ndvi = vapply(sel, cat_of, "", L = "ndvi"),
                         category_r10 = vapply(sel, cat_of, "", L = "r10"),
                         category_ndvi_own = vapply(sel, cat_of, "", L = "ndvi",
                                                    protocol = "own"),
                         category_r10_own = vapply(sel, cat_of, "", L = "r10",
                                                   protocol = "own"),
                         max_abs_rho_ndvi = vapply(sel, function(m)
                           max(abs(c_ndvi[m, setdiff(sel, m)])), 0),
                         max_abs_rho_sim = vapply(sel, function(m)
                           max(abs(c_sim[m, setdiff(sel, m)])), 0),
                         redundant_with_ndvi = vapply(sel, function(m)
                           paste(redundant_in(m, setdiff(sel, m), c_ndvi), collapse = ";"), ""),
                         redundant_with_sim = vapply(sel, function(m)
                           paste(redundant_in(m, setdiff(sel, m), c_sim), collapse = ";"), ""))
write.csv(subset_tab, here("output", "tables", "subset_selected.csv"), row.names = FALSE)
print(subset_tab)
