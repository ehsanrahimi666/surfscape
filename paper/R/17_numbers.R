## Numbers quoted in the text of the manuscript and the response letter, each
## computed from a saved output and written with its source to
## output/tables/numbers_for_text.csv.

suppressPackageStartupMessages({
  library(terra)
  library(here)
})

tab <- function(f) read.csv(here("output", "tables", f))
items <- list()
add <- function(name, value, text, source) {
  text <- gsub("(^|[^0-9A-Za-z_.])-(?=[0-9.])", "\\1−", text, perl = TRUE, useBytes = TRUE)
  items[[length(items) + 1L]] <<- data.frame(name = name, value = as.character(value),
                                             text = text, source = source)
}
s2 <- function(x) formatC(signif(as.numeric(x), 2), format = "fg", digits = 2, flag = "#") |>
  sub(pattern = "\\.$", replacement = "")
s3 <- function(x) formatC(signif(as.numeric(x), 3), format = "fg", digits = 3, flag = "#") |>
  sub(pattern = "\\.$", replacement = "")
f1 <- function(x) formatC(as.numeric(x), format = "f", digits = 1)
f2 <- function(x) formatC(as.numeric(x), format = "f", digits = 2)
f3 <- function(x) formatC(as.numeric(x), format = "f", digits = 3)
big <- function(x) formatC(x, format = "d", big.mark = ",")
sci <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("NA")
  if (x == 0) return("0")
  e <- floor(log10(abs(x)))
  m <- signif(x / 10^e, 2)
  if (abs(m) >= 10) {
    m <- m / 10
    e <- e + 1
  }
  paste0(trimws(formatC(m, format = "fg", digits = 2)), " × 10^", sub("-", "−", as.character(e)), "^")
}
rng <- function(x, f = s2) {
  if (!length(x)) return("NA")
  a <- f(min(x))
  b <- f(max(x))
  if (identical(a, b)) a else paste0(a, "–", b)
}

## ---- data and analysis window ------------------------------------------------
r13 <- rast(here("data", "raw", "NDVI-2013.tif"))
r23 <- rast(here("data", "raw", "NDVI-2023.tif"))
valid <- !is.na(values(r13)[, 1]) & !is.na(values(r23)[, 1])
add("scene_nrow", nrow(r13), as.character(nrow(r13)), "data/raw/NDVI-2013.tif")
add("scene_ncol", ncol(r13), as.character(ncol(r13)), "data/raw/NDVI-2013.tif")
add("scene_valid_both", sum(valid), big(sum(valid)), "data/raw/NDVI-2013.tif, NDVI-2023.tif")
aw <- tab("analysis_window.csv")
add("aw_feasible", aw$n_feasible_positions, big(aw$n_feasible_positions), "analysis_window.csv")
add("aw_top", aw$top_row, as.character(aw$top_row), "analysis_window.csv")
add("aw_left", aw$left_col, as.character(aw$left_col), "analysis_window.csv")
for (k in c("xmin", "xmax", "ymin", "ymax")) {
  add(paste0("aw_", k), aw[[k]], sub("^-", "−", big(aw[[k]])), "analysis_window.csv")
}
add("aw_lon", aw$centre_lon, f2(aw$centre_lon), "analysis_window.csv")
add("aw_lat", aw$centre_lat, f2(aw$centre_lat), "analysis_window.csv")
add("aw_zone", aw$utm_zone_of_centre, as.character(aw$utm_zone_of_centre), "analysis_window.csv")
add("aw_scale", aw$scale_factor_zone41, formatC(aw$scale_factor_zone41, format = "f", digits = 4),
    "analysis_window.csv")

## ---- simulated landscapes ----------------------------------------------------------
sim <- tab("simulated_landscapes.csv")
lag <- aggregate(cbind(lag1_landscape, lag1_expected) ~ range, sim, mean)
add("sim_lag1_realised", paste(lag$lag1_landscape, collapse = ";"),
    paste(f3(lag$lag1_landscape), collapse = ", "), "simulated_landscapes.csv")
add("sim_lag1_expected", paste(lag$lag1_expected, collapse = ";"),
    paste(f3(lag$lag1_expected), collapse = ", "), "simulated_landscapes.csv")
## circulant embedding as in surfscape:::.ss_grf (grid of twice the size)
min_eig <- vapply(c(3, 10, 30), function(rg) {
  m1 <- 2L * 201L
  i <- 0:(m1 - 1L)
  d <- pmin(i, m1 - i)
  min(Re(stats::fft(exp(-3 * sqrt(outer(d^2, d^2, "+")) / rg))))
}, numeric(1))
add("sim_min_eigen", paste(min_eig, collapse = ";"), paste(s3(min_eig), collapse = ", "),
    "computed from the covariance of the simulation (R/17_numbers.R)")

## ---- main benchmark --------------------------------------------------------------------
summ <- tab("benchmark_main_summary.csv")
resp <- tab("benchmark_main_responses.csv")
fx <- summ[summ$protocol == "fixed", ]
own <- summ[summ$protocol == "own", ]
cats <- c("compositional", "configurational", "mixed", "unresponsive")
cat_short <- c(compositional = "comp", configurational = "conf", mixed = "mixed",
               unresponsive = "unresp")
for (L in c("r03", "r10", "r30", "ndvi")) {
  k <- table(factor(fx$category[fx$landscape == L], cats))
  for (cc in cats) add(paste0("n_", cat_short[[cc]], "_", L), k[[cc]], as.character(k[[cc]]),
                       paste0("benchmark_main_summary.csv, fixed, ", L))
  add(paste0("n_impl_", L), sum(k), as.character(sum(k)), "benchmark_main_summary.csv")
}
same <- all(vapply(unique(fx$metric), function(m)
  length(unique(fx$category[fx$metric == m & fx$landscape %in% c("r03", "r10", "r30")])) == 1L,
  logical(1)))
add("sims_same_class", same, if (same) "yes" else "no", "benchmark_main_summary.csv")
same_ndvi <- all(vapply(unique(fx$metric), function(m)
  length(unique(fx$category[fx$metric == m])) == 1L, logical(1)))
add("all_same_class", same_ndvi, if (same_ndvi) "yes" else "no", "benchmark_main_summary.csv")

rv <- function(L, m, test, field = "response", protocol = "fixed") {
  x <- resp[resp$landscape == L & resp$metric == m & resp$test == test &
              resp$protocol == protocol, ]
  if (!nrow(x)) return(NA_real_)
  if (field == "response" && isTRUE(x$invariant[1])) 0 else x[[field]][1]
}
sv <- function(L, m, field, protocol = "fixed") {
  x <- summ[summ$landscape == L & summ$metric == m & summ$protocol == protocol, ]
  if (!nrow(x)) NA else x[[field]][1]
}

## largest change of a compositional metric under permutation, relative to its values
obs <- readRDS(here("output", "rds", "benchmark_main_observed.rds"))
rel <- unlist(lapply(c("r03", "r10", "r30", "ndvi"), function(L) {
  comp <- fx$metric[fx$landscape == L & fx$category == "compositional"]
  vapply(comp, function(m) rv(L, m, "permute", "max_abs") /
           max(abs(obs[[L]][, m]), na.rm = TRUE), numeric(1))
}))
add("comp_perm_maxrel", max(rel, na.rm = TRUE), sci(max(rel, na.rm = TRUE)),
    "benchmark_main_responses.csv (max_abs, permute) / benchmark_main_observed.rds")

## responses of named metrics
named <- c("surfscape_sr1", "terra_moran", "terra_geary", "surfscape_gcor",
           "glcmtex_range_correlation", "surfscape_gcon", "glcmtex_range_contrast",
           "surfscape_slac", "surfscape_sdq", "surfscape_sdr", "geodiv_sdr",
           "glcmtex_range_mean", "glcmtex_range_variance", "glcmtex_prob_mean",
           "glcmtex_prob_variance", "glcmtex_range_entropy", "glcmtex_range_ASM",
           "geodiv_scl20", "geodiv_str20", "geodiv_s10z", "surfscape_s10z", "geodiv_sdq6",
           "glcmtex_prob_correlation")
tests_all <- c("permute", "shift", "scale", "reshape", "target", "rearrange")
for (L in c("r03", "r10", "r30", "ndvi")) for (m in named) for (tt in tests_all) {
  v <- rv(L, m, tt)
  if (!is.na(v)) add(paste(m, tt, L, sep = "__"), v, s2(v),
                     paste0("benchmark_main_responses.csv, fixed, ", L))
}
for (L in c("r03", "r10", "r30", "ndvi")) for (m in named) {
  add(paste(m, "signif", L, sep = "__"), sv(L, m, "signif"), f2(sv(L, m, "signif")),
      paste0("benchmark_main_summary.csv, fixed, ", L))
  add(paste(m, "category", L, sep = "__"), sv(L, m, "category"), sv(L, m, "category"),
      paste0("benchmark_main_summary.csv, fixed, ", L))
}

## every arrangement-sensitive implementation also responds to the reshape test
mx <- fx[fx$category == "mixed", ]
all_reshape <- all(!mx$reshape_invariant)
add("mixed_all_reshape_sensitive", all_reshape, if (all_reshape) "yes" else "no",
    "benchmark_main_summary.csv, fixed")
own_conf <- unique(own$metric[own$category == "configurational"])
add("own_conf_all_landscapes", paste(own_conf, collapse = ";"), paste(own_conf, collapse = ", "),
    "benchmark_main_summary.csv, own")
own_conf_ls <- tapply(own$category == "configurational", own$landscape, sum)
add("own_n_conf_by_landscape", paste(names(own_conf_ls), own_conf_ls, collapse = ";"),
    paste(names(own_conf_ls), own_conf_ls, collapse = ", "), "benchmark_main_summary.csv, own")

## ratio of the composition response to the permutation response among mixed metrics
for (L in c("r03", "r10", "r30", "ndvi")) {
  x <- fx[fx$landscape == L & fx$category == "mixed", ]
  ratio <- x$composition_response / x$configuration_response
  x <- x[is.finite(ratio), ]
  ratio <- ratio[is.finite(ratio)]
  o <- order(ratio)
  add(paste0("ratio_lowest5_", L), paste(x$metric[o][1:5], collapse = ";"),
      paste0(x$metric[o][1:5], " (", s2(ratio[o][1:5]), ")", collapse = ", "),
      paste0("benchmark_main_summary.csv, fixed, ", L))
  add(paste0("ratio_highest5_", L), paste(rev(x$metric[o])[1:5], collapse = ";"),
      paste0(rev(x$metric[o])[1:5], " (", s2(rev(ratio[o])[1:5]), ")", collapse = ", "),
      paste0("benchmark_main_summary.csv, fixed, ", L))
}

for (L in c("r03", "r10", "r30", "ndvi")) {
  x <- fx[fx$landscape == L & fx$category == "mixed" & fx$shift_invariant & fx$scale_invariant, ]
  add(paste0("affine_invariant_mixed_", L), paste(x$metric, collapse = ";"),
      paste(x$metric, collapse = ", "), paste0("benchmark_main_summary.csv, fixed, ", L))
}

## ranges over groups of metrics
grp_range <- function(L, mets, test, nm) {
  v <- vapply(mets, function(m) rv(L, m, test), numeric(1))
  v <- v[is.finite(v)]
  add(nm, paste(range(v), collapse = ";"), rng(v), paste0("benchmark_main_responses.csv, ", L))
  invisible(v)
}
comp_nd <- fx$metric[fx$landscape == "ndvi" & fx$category == "compositional"]
v <- grp_range("ndvi", comp_nd, "target", "target_comp_ndvi_range")
add("target_comp_ndvi_median", median(v), s2(median(v)), "benchmark_main_responses.csv, ndvi")
v <- grp_range("ndvi", comp_nd, "rearrange", "rearr_comp_ndvi_range")
add("rearr_comp_ndvi_median", median(v), s2(median(v)), "benchmark_main_responses.csv, ndvi")
for (L in c("r03", "r10", "r30")) {
  comp_L <- fx$metric[fx$landscape == L & fx$category == "compositional"]
  v <- grp_range(L, comp_L, "rearrange", paste0("rearr_comp_", L, "_range"))
}
autoc <- c("surfscape_sr1", "terra_moran", "terra_geary")
for (L in c("r03", "r10", "r30", "ndvi")) {
  grp_range(L, autoc, "permute", paste0("autoc_perm_", L, "_range"))
  grp_range(L, autoc, "reshape", paste0("autoc_reshape_", L, "_range"))
}
grp_range("ndvi", autoc, "target", "autoc_target_ndvi_range")
mixed_nd <- fx$metric[fx$landscape == "ndvi" & fx$category == "mixed"]
grp_range("ndvi", mixed_nd, "rearrange", "rearr_mixed_ndvi_range")
grm <- c("glcmtex_range_mean", "glcmtex_range_variance", "glcmpkg_mean", "glcmpkg_variance")
for (L in c("r03", "r10", "r30", "ndvi")) grp_range(L, grm, "permute", paste0("glcm_meanvar_perm_", L))
grc <- c("glcmtex_range_correlation", "glcmtex_range_contrast", "glcmpkg_correlation",
         "glcmpkg_contrast", "surfscape_gcon", "surfscape_gcor")
for (L in c("r03", "r10", "r30", "ndvi")) grp_range(L, grc, "permute", paste0("glcm_concor_perm_", L))

## own-range protocol
bin <- read.csv(here("output", "tables", "metric_catalogue.csv"))
binned <- bin$implementation[bin$binned]
own_nd <- own[own$landscape == "ndvi" & own$metric %in% binned, ]
k_own <- table(factor(own_nd$category, cats))
for (cc in cats) add(paste0("own_n_", cat_short[[cc]]), k_own[[cc]], as.character(k_own[[cc]]),
                     "benchmark_main_summary.csv, own, ndvi")
add("own_n_binned", nrow(own_nd), as.character(nrow(own_nd)), "benchmark_main_summary.csv, own, ndvi")
aff <- own_nd$metric[own_nd$shift_invariant & own_nd$scale_invariant]
add("own_n_affine_invariant", length(aff), as.character(length(aff)),
    "benchmark_main_summary.csv, own, ndvi (shift and stretch invariant)")
add("own_affine_not", paste(setdiff(own_nd$metric, aff), collapse = ";"),
    paste(setdiff(own_nd$metric, aff), collapse = ", "), "benchmark_main_summary.csv, own, ndvi")
add("own_configurational", paste(own_nd$metric[own_nd$category == "configurational"], collapse = ";"),
    paste(own_nd$metric[own_nd$category == "configurational"], collapse = ", "),
    "benchmark_main_summary.csv, own, ndvi")
for (m in c("glcmtex_prob_mean", "glcmtex_prob_variance", "glcmtex_prob_contrast",
            "glcmtex_prob_correlation")) {
  v <- rv("ndvi", m, "permute", protocol = "own")
  add(paste0(m, "__permute__ndvi_own"), v, s2(v), "benchmark_main_responses.csv, own, ndvi")
}

## rasterdiv class limit
rc <- tab("rasterdiv_class_check.csv")
add("rasterdiv_n_over_100", rc$n_over_100, as.character(rc$n_over_100), "rasterdiv_class_check.csv")
add("rasterdiv_max_distinct", rc$max_distinct, as.character(rc$max_distinct), "rasterdiv_class_check.csv")

## geodiv metrics on centred windows
gc <- tab("benchmark_geodiv_centred_summary.csv")
gr <- tab("benchmark_geodiv_centred_responses.csv")
for (vv in c("raw", "centred")) for (m in c("geodiv_scl20", "geodiv_str20")) {
  x <- gc[gc$variant == vv & gc$metric == m, ]
  add(paste(m, vv, "category", sep = "__"), x$category, x$category,
      "benchmark_geodiv_centred_summary.csv")
  inv <- c(shift = x$shift_invariant, stretch = x$scale_invariant,
           reshape = x$reshape_invariant, observed = x$target_invariant)
  add(paste(m, vv, "invariant_to", sep = "__"), paste(names(inv)[inv], collapse = ";"),
      if (any(inv)) paste(names(inv)[inv], collapse = ", ") else "none",
      "benchmark_geodiv_centred_summary.csv")
  for (tt in c("permute", "shift", "scale", "reshape", "target")) {
    y <- gr[gr$variant == vv & gr$metric == m & gr$test == tt, ]
    val <- if (isTRUE(y$invariant)) 0 else y$response
    add(paste(m, vv, tt, sep = "__"), val, s2(val), "benchmark_geodiv_centred_responses.csv")
  }
}

## ---- extent and grain ------------------------------------------------------------------
eg <- tab("benchmark_extent_grain_summary.csv")
er <- tab("benchmark_extent_grain_responses.csv")
eg <- eg[eg$protocol == "fixed", ]
er <- er[er$protocol == "fixed", ]
comp_main <- unique(fx$metric[fx$category == "compositional"])
comp_eg <- unique(eg$metric[eg$metric %in% comp_main])
inv_all <- all(eg$permute_invariant[eg$metric %in% comp_eg])
add("extent_comp_always_invariant", inv_all, if (inv_all) "yes" else "no",
    "benchmark_extent_grain_summary.csv")
cls_same <- all(vapply(unique(eg$metric), function(m) length(unique(eg$category[eg$metric == m])) == 1L,
                       logical(1)))
add("extent_class_same", cls_same, if (cls_same) "yes" else "no", "benchmark_extent_grain_summary.csv")
changed <- unique(eg$metric[vapply(seq_len(nrow(eg)), function(i) {
  m <- eg$metric[i]
  length(unique(eg$category[eg$metric == m])) > 1L
}, logical(1))])
add("extent_class_changed", paste(changed, collapse = ";"),
    if (length(changed)) paste(changed, collapse = ", ") else "none",
    "benchmark_extent_grain_summary.csv")
ev <- function(g, lv, m, field) {
  x <- eg[eg$group == g & eg$level == lv & eg$metric == m, ]
  if (!nrow(x)) NA else x[[field]][1]
}
for (g in c("r03", "r10", "r30", "ndvi30")) for (m in c("surfscape_sr1", "terra_moran",
  "surfscape_gcon", "glcmtex_range_correlation", "surfscape_slac", "surfscape_sdr")) {
  for (lv in c("w5", "w9", "w15", "w31", "w61", "landscape")) {
    add(paste(m, "perm", g, lv, sep = "__"), ev(g, lv, m, "configuration_response"),
        s2(ev(g, lv, m, "configuration_response")), "benchmark_extent_grain_summary.csv")
    add(paste(m, "signif", g, lv, sep = "__"), ev(g, lv, m, "signif"),
        f2(ev(g, lv, m, "signif")), "benchmark_extent_grain_summary.csv")
  }
}
## change of the permutation response and of the proportion of significant windows with
## window size, mixed implementations (5-61 cells)
lvs <- c("w5", "w9", "w15", "w31", "w61")
mono <- do.call(rbind, lapply(setdiff(unique(eg$metric[eg$category == "mixed"]), changed), function(m) {
  do.call(rbind, lapply(c("r03", "r10", "r30", "ndvi30"), function(g) {
    x <- eg[eg$group == g & eg$metric == m, ]
    x <- x[match(lvs, x$level), ]
    v <- x$configuration_response
    if (any(!is.finite(v))) return(NULL)
    data.frame(metric = m, group = g, up = all(diff(v) >= 0), down = all(diff(v) <= 0),
               sig_up = all(diff(x$signif) >= 0))
  }))
}))
add("ext_perm_combos_n", nrow(mono), as.character(nrow(mono)), "benchmark_extent_grain_summary.csv")
add("ext_perm_up_n", sum(mono$up), as.character(sum(mono$up)), "benchmark_extent_grain_summary.csv")
add("ext_sig_up_n", sum(mono$sig_up), as.character(sum(mono$sig_up)), "benchmark_extent_grain_summary.csv")
mv <- grepl("_(mean|variance)$", mono$metric)
add("ext_meanvar_down_all", all(mono$down[mv]), if (all(mono$down[mv])) "yes" else "no",
    "benchmark_extent_grain_summary.csv")
oth <- mono[!mono$up & !mv, ]
add("ext_perm_other_exceptions", paste(oth$metric, oth$group, sep = ":", collapse = ";"),
    paste(oth$metric, oth$group, collapse = ", "), "benchmark_extent_grain_summary.csv")
oth_s <- mono[!mono$sig_up, ]
add("ext_sig_exceptions", paste(oth_s$metric, oth_s$group, sep = ":", collapse = ";"),
    paste(oth_s$metric, oth_s$group, collapse = ", "), "benchmark_extent_grain_summary.csv")
## increase of the standardised response from 5- to 61-cell windows split, on the log
## scale, into the increase of the mean absolute change and the decrease of the
## standard deviation among windows
pr <- er[er$test == "permute" & !er$invariant, ]
parts <- do.call(rbind, lapply(c("r03", "r10", "r30", "ndvi30"), function(g) {
  do.call(rbind, lapply(setdiff(unique(pr$metric), changed), function(m) {
    a <- pr[pr$group == g & pr$metric == m & pr$level == "w5", ]
    b <- pr[pr$group == g & pr$metric == m & pr$level == "w61", ]
    if (!nrow(a) || !nrow(b) || !is.finite(a$response) || !is.finite(b$response)) return(NULL)
    data.frame(group = g, metric = m, response_ratio = b$response / a$response,
               log_abs_part = log(b$mean_abs / a$mean_abs),
               log_sd_part = -log((b$mean_abs / b$response) / (a$mean_abs / a$response)))
  }))
}))
write.csv(parts, here("output", "tables", "extent_response_parts.csv"), row.names = FALSE)
upp <- parts[parts$response_ratio > 1, ]
add("ext_up61_n", nrow(upp), as.character(nrow(upp)), "extent_response_parts.csv")
add("ext_up61_sdmain_n", sum(upp$log_sd_part > upp$log_abs_part),
    as.character(sum(upp$log_sd_part > upp$log_abs_part)), "extent_response_parts.csv")
## absolute change caused by permutation and standard deviation among windows
for (g in c("r03", "r10", "r30", "ndvi30")) for (m in c("surfscape_sr1", "surfscape_gcon")) {
  for (lv in c("w5", "w9", "w15", "w31", "w61")) {
    x <- er[er$group == g & er$level == lv & er$metric == m & er$test == "permute", ]
    if (!nrow(x)) next
    add(paste(m, "absperm", g, lv, sep = "__"), x$mean_abs, s2(x$mean_abs),
        "benchmark_extent_grain_responses.csv (mean_abs)")
    add(paste(m, "sdwin", g, lv, sep = "__"), x$mean_abs / x$response, s2(x$mean_abs / x$response),
        "benchmark_extent_grain_responses.csv (mean_abs / response)")
  }
}
## landscape-level values of the simulated landscapes, whose mean is zero by construction
lvv <- tab("benchmark_landscape_level_values.csv")
lvv <- lvv[lvv$set == "builtin_fixed" & lvv$group %in% c("r03", "r10", "r30"), ]
x <- abs(lvv$value[lvv$metric == "surfscape_smean"])
add("sim_land_kernel_absmean_max", max(x), sci(max(x)),
    "benchmark_landscape_level_values.csv (surfscape_smean, simulated landscapes)")
x <- abs(lvv$value[lvv$metric == "surfscape_scv"])
add("sim_land_scv_abs_range", paste(range(x), collapse = ";"), paste0(sci(min(x)), " to ", sci(max(x))),
    "benchmark_landscape_level_values.csv (surfscape_scv, simulated landscapes)")
## landscape-level significance of mixed metrics
for (g in c("r03", "r10", "r30", "ndvi30")) {
  x <- eg[eg$group == g & eg$level == "landscape" & eg$category == "mixed", ]
  add(paste0("land_mixed_signif_", g), paste(range(x$signif), collapse = ";"),
      rng(x$signif, f2), "benchmark_extent_grain_summary.csv")
  add(paste0("land_mixed_n_", g), nrow(x), as.character(nrow(x)), "benchmark_extent_grain_summary.csv")
  x <- x[!x$metric %in% changed, ]
  add(paste0("land_mixed_n_excl_", g), nrow(x), as.character(nrow(x)),
      "benchmark_extent_grain_summary.csv (without metrics whose category changed)")
  add(paste0("land_sig_all_n_", g), sum(x$signif == 1), as.character(sum(x$signif == 1)),
      "benchmark_extent_grain_summary.csv")
  ex <- x[x$signif < 1, ]
  add(paste0("land_sig_exceptions_", g), paste(ex$metric, ex$signif, sep = ":", collapse = ";"),
      if (nrow(ex)) paste0(ex$metric, " (", f1(ex$signif), ")", collapse = ", ") else "none",
      "benchmark_extent_grain_summary.csv")
  add(paste0("land_sig_exc_range_", g),
      if (nrow(ex)) paste(range(ex$signif), collapse = ";") else NA,
      rng(ex$signif, f1), "benchmark_extent_grain_summary.csv")
}
x <- eg[eg$group %in% c("r03", "r10", "r30") & eg$level == "landscape" & eg$category == "mixed" &
          !eg$metric %in% changed & eg$signif < 1, ]
add("land_sig_exc_range_sim", paste(range(x$signif), collapse = ";"), rng(x$signif, f1),
    "benchmark_extent_grain_summary.csv")
add("land_sig_exc_pct_sim", paste(100 * range(x$signif), collapse = ";"),
    rng(100 * x$signif, function(v) formatC(v, format = "f", digits = 0)),
    "benchmark_extent_grain_summary.csv (percentage of landscapes)")
add("land_sig_exc_metrics_sim", paste(sort(unique(x$metric)), collapse = ";"),
    paste(sort(unique(x$metric)), collapse = ", "), "benchmark_extent_grain_summary.csv")
## own-range protocol: GLCMTextures textures with equal-probability grey levels
eo <- tab("benchmark_extent_grain_summary.csv")
ro <- tab("benchmark_extent_grain_responses.csv")
eo <- eo[eo$protocol == "own" & grepl("^glcmtex_prob_", eo$metric), ]
ro <- ro[ro$protocol == "own" & grepl("^glcmtex_prob_", ro$metric), ]
ok_lv <- aggregate(conf ~ group + level, transform(eo, conf = category == "configurational"), all)
add("own_prob_conf_all_levels", all(ok_lv$conf), if (all(ok_lv$conf)) "yes" else "no",
    "benchmark_extent_grain_summary.csv (protocol own)")
add("own_prob_not_conf_at", paste(ok_lv$group[!ok_lv$conf], ok_lv$level[!ok_lv$conf], sep = ":",
                                  collapse = ";"),
    paste(ok_lv$group[!ok_lv$conf], ok_lv$level[!ok_lv$conf], collapse = ", "),
    "benchmark_extent_grain_summary.csv (protocol own)")
nc_tests <- unique(unlist(lapply(c("shift", "scale", "reshape", "target"), function(tt) {
  k <- paste0(tt, "_invariant")
  if (k %in% names(eo) && any(eo[[k]] %in% FALSE)) tt
})))
add("own_prob_noninvariant_tests", paste(nc_tests, collapse = ";"), paste(nc_tests, collapse = ", "),
    "benchmark_extent_grain_summary.csv (protocol own)")
tg <- ro[ro$test == "target" & !ro$invariant, c("group", "level", "metric", "max_abs")]
pm <- ro[ro$test == "permute", c("group", "level", "metric", "mean_abs")]
tp <- merge(tg, pm, by = c("group", "level", "metric"))
add("own_prob_target_rel_perm_max", max(tp$max_abs / tp$mean_abs),
    sci(max(tp$max_abs / tp$mean_abs)),
    "benchmark_extent_grain_responses.csv (protocol own; largest change under the target test / mean change under permutation)")
tw <- tab("benchmark_own_target_windows.csv")
tw <- tw[grepl("^glcmtex_prob_", tw$metric), ]
for (gl in list(c("ndvi30", "w61"), c("ndvi30", "landscape"), c("ndvi120", "w15"), c("ndvi120", "w7"))) {
  x <- tw$n_changed[tw$group == gl[1] & tw$level == gl[2]]
  add(paste("own_target_nchanged", gl[1], gl[2], sep = "__"), max(x), as.character(max(x)),
      "benchmark_own_target_windows.csv (largest over the eight textures)")
}
tc <- tab("target_ties_check.csv")
for (i in seq_len(nrow(tc))) {
  gr_i <- tc$grain[i]
  add(paste0("ties_cells_", gr_i), tc$cells[i], big(tc$cells[i]), "target_ties_check.csv")
  add(paste0("ties_tied2013_", gr_i), tc$cells_tied_2013[i], big(tc$cells_tied_2013[i]),
      "target_ties_check.csv")
  add(paste0("ties_split_", gr_i), tc$ties_split[i], big(tc$ties_split[i]), "target_ties_check.csv")
  add(paste0("ties_created_", gr_i), tc$ties_created[i], big(tc$ties_created[i]),
      "target_ties_check.csv")
  add(paste0("ties_glchanged_", gr_i), tc$cells_grey_level_changed[i],
      big(tc$cells_grey_level_changed[i]), "target_ties_check.csv")
}
## landscape-wide rearrangement: compositional metrics inside windows and at the landscape level
for (g in c("r03", "r10", "r30", "ndvi30")) for (m in c("surfscape_sq", "surfscape_ssk",
                                                        "surfscape_sent")) {
  for (lv in c("w5", "w15", "w61", "landscape")) {
    x <- er[er$group == g & er$level == lv & er$metric == m & er$test == "rearrange", ]
    val <- if (!nrow(x)) NA else if (isTRUE(x$invariant)) 0 else x$response
    add(paste(m, "rearr", g, lv, sep = "__"), val, s2(val), "benchmark_extent_grain_responses.csv")
    if (lv == "landscape" && nrow(x)) {
      add(paste(m, "rearr_inv", g, lv, sep = "__"), x$invariant, as.character(x$invariant),
          "benchmark_extent_grain_responses.csv")
    }
  }
}
land_comp_inv <- all(er$invariant[er$level == "landscape" & er$test == "rearrange" &
                                    er$metric %in% setdiff(comp_eg, changed)])
add("land_comp_rearr_invariant", land_comp_inv, if (land_comp_inv) "yes" else "no",
    "benchmark_extent_grain_responses.csv (without metrics whose category changed)")
## grain
for (g in c("ndvi30", "ndvi60", "ndvi120")) {
  x <- eg[eg$group == g & eg$level == "w15", ]
  k <- table(factor(x$category, cats))
  add(paste0("grain_", g, "_w15_counts"), paste(k, collapse = ";"),
      paste0(k[["compositional"]], " compositional, ", k[["mixed"]], " mixed, ",
             k[["configurational"]], " configurational"), "benchmark_extent_grain_summary.csv")
}
gm <- unique(eg$metric[eg$group == "ndvi60"])
g_same <- all(vapply(gm, function(m) length(unique(eg$category[eg$metric == m &
  eg$group %in% c("ndvi30", "ndvi60", "ndvi120")])) == 1L, logical(1)))
add("grain_class_same", g_same, if (g_same) "yes" else "no", "benchmark_extent_grain_summary.csv")
grain_m <- c("surfscape_sr1", "terra_moran", "terra_geary", "surfscape_gcon", "surfscape_slac",
             "surfscape_sdq", "glcmtex_prob_correlation", "glcmtex_range_contrast")
for (m in grain_m) {
  for (gl in list(c("ndvi30", "w15"), c("ndvi60", "w15"), c("ndvi120", "w15"),
                  c("ndvi30", "w31"), c("ndvi120", "w7"))) {
    add(paste(m, "grain", gl[1], gl[2], sep = "__"), ev(gl[1], gl[2], m, "configuration_response"),
        s2(ev(gl[1], gl[2], m, "configuration_response")), "benchmark_extent_grain_summary.csv")
    add(paste(m, "graincomp", gl[1], gl[2], sep = "__"), ev(gl[1], gl[2], m, "composition_response"),
        s2(ev(gl[1], gl[2], m, "composition_response")), "benchmark_extent_grain_summary.csv")
    add(paste(m, "grainsig", gl[1], gl[2], sep = "__"), ev(gl[1], gl[2], m, "signif"),
        f2(ev(gl[1], gl[2], m, "signif")), "benchmark_extent_grain_summary.csv")
  }
}
## ratio of the composition response to the configuration response, 15-cell windows
gr_ratio <- do.call(rbind, lapply(c("ndvi30", "ndvi60", "ndvi120"), function(g) {
  x <- eg[eg$group == g & eg$level == "w15" & eg$category == "mixed", ]
  data.frame(group = g, metric = x$metric, ratio = x$composition_response / x$configuration_response)
}))
gw <- reshape(gr_ratio, idvar = "metric", timevar = "group", direction = "wide")
gw <- gw[stats::complete.cases(gw), ]
up <- gw$ratio.ndvi120 > gw$ratio.ndvi30
add("grain_ratio_n_mixed", nrow(gw), as.character(nrow(gw)), "benchmark_extent_grain_summary.csv")
add("grain_ratio_n_up_30_120", sum(up), as.character(sum(up)),
    "benchmark_extent_grain_summary.csv (ratio composition/configuration response, w15)")
mono_up <- gw$ratio.ndvi30 < gw$ratio.ndvi60 & gw$ratio.ndvi60 < gw$ratio.ndvi120
mono_dn <- gw$ratio.ndvi30 > gw$ratio.ndvi60 & gw$ratio.ndvi60 > gw$ratio.ndvi120
add("grain_ratio_n_mono_up", sum(mono_up), as.character(sum(mono_up)),
    "benchmark_extent_grain_summary.csv")
add("grain_ratio_n_mono_down", sum(mono_dn), as.character(sum(mono_dn)),
    "benchmark_extent_grain_summary.csv")
add("grain_ratio_n_nonmono", sum(!mono_up & !mono_dn), as.character(sum(!mono_up & !mono_dn)),
    "benchmark_extent_grain_summary.csv")
rk <- stats::cor(gw$ratio.ndvi30, gw$ratio.ndvi120, method = "spearman")
add("grain_ratio_rho_30_120", rk, f2(rk), "benchmark_extent_grain_summary.csv")
write.csv(gw, here("output", "tables", "grain_ratio_w15.csv"), row.names = FALSE)
x <- unlist(gw[gw$metric == "surfscape_sr1", c("ratio.ndvi30", "ratio.ndvi60", "ratio.ndvi120")])
add("sr1_grain_ratio_max", max(x), s2(max(x)), "grain_ratio_w15.csv")

## ---- redundancy and key subset ------------------------------------------------------
cn <- as.matrix(read.csv(here("output", "tables", "redundancy_spearman_ndvi.csv"), row.names = 1,
                         check.names = FALSE))
cs <- as.matrix(read.csv(here("output", "tables", "redundancy_spearman_sim.csv"), row.names = 1,
                         check.names = FALSE))
add("red_n_metrics", ncol(cn), as.character(ncol(cn)), "redundancy_spearman_ndvi.csv")
cl <- tab("redundancy_clusters.csv")
add("red_n_clusters", length(unique(cl$cluster)), as.character(length(unique(cl$cluster))),
    "redundancy_clusters.csv")
ex <- tab("redundancy_excluded.csv")
add("red_excluded", paste(ex$metric, collapse = ";"),
    if (nrow(ex)) paste(ex$metric, collapse = ", ") else "none", "redundancy_excluded.csv")
add("rho_sent_sq_ndvi", cn["surfscape_sent", "surfscape_sq"], f2(cn["surfscape_sent", "surfscape_sq"]),
    "redundancy_spearman_ndvi.csv")
add("rho_srao_sq_ndvi", cn["surfscape_srao", "surfscape_sq"], f2(cn["surfscape_srao", "surfscape_sq"]),
    "redundancy_spearman_ndvi.csv")
add("rho_sent_sq_sim", cs["surfscape_sent", "surfscape_sq"], f2(cs["surfscape_sent", "surfscape_sq"]),
    "redundancy_spearman_sim.csv")
for (pr in list(c("surfscape_smean", "surfscape_sq"), c("terra_moran", "glcmtex_prob_contrast"),
                c("terra_moran", "surfscape_sr1"))) {
  nm <- paste0("rho_", sub("^[^_]+_", "", pr[1]), "_", sub("^[^_]+_", "", pr[2]))
  add(paste0(nm, "_ndvi"), cn[pr[1], pr[2]], f2(cn[pr[1], pr[2]]), "redundancy_spearman_ndvi.csv")
  add(paste0(nm, "_sim"), cs[pr[1], pr[2]], f2(cs[pr[1], pr[2]]), "redundancy_spearman_sim.csv")
}
strict <- tab("subset_strict_rule.csv")
add("subset_strict", paste(strict$metric, collapse = ";"), paste(strict$metric, collapse = ", "),
    "subset_strict_rule.csv")
sel <- tab("subset_selected.csv")
add("subset_metrics", paste(sel$metric, collapse = ";"), paste(sel$metric, collapse = ", "),
    "subset_selected.csv")
add("subset_max_rho", max(c(sel$max_abs_rho_ndvi, sel$max_abs_rho_sim)),
    f2(max(c(sel$max_abs_rho_ndvi, sel$max_abs_rho_sim))), "subset_selected.csv")

## ---- evenness and room for configuration -----------------------------------------------
evn <- tab("evenness_configuration_room.csv")
for (m in evn$metric) {
  x <- evn[evn$metric == m, ]
  add(paste0("rho_even_sdperm_", m), x$rho_evenness_sdperm, f2(x$rho_evenness_sdperm),
      "evenness_configuration_room.csv")
  add(paste0("rho_even_ses_", m), x$rho_evenness_abs_ses, f2(x$rho_evenness_abs_ses),
      "evenness_configuration_room.csv")
  add(paste0("rho_sq_sdperm_", m), x$rho_sq_sdperm, f2(x$rho_sq_sdperm),
      "evenness_configuration_room.csv")
  add(paste0("even_n_complete_", m), x$n_complete, as.character(x$n_complete),
      "evenness_configuration_room.csv")
  add(paste0("even_n_undefined_", m), x$n_undefined, as.character(x$n_undefined),
      "evenness_configuration_room.csv")
}
add("rho_even_sq", evn$rho_evenness_sq[1], f2(evn$rho_evenness_sq[1]), "evenness_configuration_room.csv")

## ---- decomposition ------------------------------------------------------------------------
dv <- tab("decomposition_virtual_windows.csv")
add("virt_max_residual_excl_scv", max(dv$max_residual[dv$metric != "scv"]),
    sci(max(dv$max_residual[dv$metric != "scv"])), "decomposition_virtual_windows.csv")
rel_res <- dv$max_residual / pmax(dv$mad_total, 1e-300)
add("virt_max_rel_residual", max(rel_res), sci(max(rel_res)), "decomposition_virtual_windows.csv")
add("virt_max_rel_residual_excl_scv", max(rel_res[dv$metric != "scv"]),
    sci(max(rel_res[dv$metric != "scv"])), "decomposition_virtual_windows.csv")
add("virt_max_rel_residual_scv", max(rel_res[dv$metric == "scv"]),
    sci(max(rel_res[dv$metric == "scv"])), "decomposition_virtual_windows.csv")
z1 <- dv[dv$lambda == 0 & dv$kappa > 0, ]
add("virt_lambda0_max_config_local", max(c(z1$max_abs_config, z1$max_abs_local)),
    s2(max(c(z1$max_abs_config, z1$max_abs_local))), "decomposition_virtual_windows.csv")
z2 <- dv[dv$kappa == 0 & dv$lambda > 0, ]
add("virt_kappa0_max_global", max(z2$max_abs_global), s2(max(z2$max_abs_global)),
    "decomposition_virtual_windows.csv")
comp_b <- c("smean", "sa", "sq", "ssk", "sku", "sp", "sv", "sz", "s10z", "scv", "sent", "srao")
cm <- dv[dv$metric %in% comp_b, ]
cfg_r <- cm$max_abs_config / pmax(cm$max_abs_global, cm$max_abs_local)
cfg_rel <- max(cfg_r, na.rm = TRUE)
add("virt_comp_config_rel", cfg_rel, sci(cfg_rel), "decomposition_virtual_windows.csv")
add("virt_comp_config_rel_excl_scv", max(cfg_r[cm$metric != "scv"], na.rm = TRUE),
    sci(max(cfg_r[cm$metric != "scv"], na.rm = TRUE)), "decomposition_virtual_windows.csv")
add("virt_comp_config_rel_scv", max(cfg_r[cm$metric == "scv"], na.rm = TRUE),
    sci(max(cfg_r[cm$metric == "scv"], na.rm = TRUE)), "decomposition_virtual_windows.csv")
ag <- aggregate(cbind(mad_global, mad_local, mad_config, mad_total, rank_cor) ~ kappa + lambda +
                  metric, dv, mean)
for (m in c("sr1", "gcon", "sq", "sent", "slac")) for (k in 0:2) for (l in c(0, 0.5, 1)) {
  x <- ag[ag$metric == m & ag$kappa == k & ag$lambda == l, ]
  for (cp in c("mad_global", "mad_local", "mad_config", "mad_total")) {
    add(paste("virt", m, cp, k, l, sep = "__"), x[[cp]], s2(x[[cp]]),
        "decomposition_virtual_windows.csv (mean over landscapes)")
  }
}
glob_dec <- unique(ag$metric[vapply(seq_len(nrow(ag)), function(i) {
  if (ag$kappa[i] != 2) return(FALSE)
  a1 <- ag$mad_global[ag$metric == ag$metric[i] & ag$kappa == 1 & ag$lambda == ag$lambda[i]]
  ag$mad_global[i] < a1
}, logical(1))])
add("virt_global_not_increasing", paste(glob_dec, collapse = ";"),
    if (length(glob_dec)) paste(glob_dec, collapse = ", ") else "none",
    "decomposition_virtual_windows.csv (mad_global smaller at kappa 2 than at kappa 1)")
noncomp_b <- setdiff(unique(ag$metric), comp_b)
cfg_mono <- all(vapply(noncomp_b, function(m) vapply(0:2, function(k) {
  x <- ag[ag$metric == m & ag$kappa == k, ]
  all(diff(x$mad_config[order(x$lambda)]) > 0)
}, logical(1)) |> all(), logical(1)))
add("virt_config_increases_with_lambda", cfg_mono, if (cfg_mono) "yes" else "no",
    "decomposition_virtual_windows.csv")
for (l in c(0.25, 0.5, 0.75, 1)) {
  x <- ag[ag$metric == "sq" & ag$kappa == 0 & ag$lambda == l, ]
  add(paste0("virt_rankcor_", l), x$rank_cor, f2(x$rank_cor), "decomposition_virtual_windows.csv")
}
dl <- tab("decomposition_ndvi_landscape.csv")
dl$metric <- ifelse(dl$set == "builtin", paste0("surfscape_", dl$metric), dl$metric)
dl$configuration[abs(dl$configuration) <= 1e-9 * abs(dl$total)] <- 0
dl$share_config[dl$configuration == 0] <- 0
for (i in seq_len(nrow(dl))) {
  add(paste0("land_share_", dl$metric[i]), dl$share_config[i], f1(100 * dl$share_config[i]),
      "decomposition_ndvi_landscape.csv")
}
add("sr1_share_land", dl$share_config[dl$metric == "surfscape_sr1"],
    f1(100 * dl$share_config[dl$metric == "surfscape_sr1"]), "decomposition_ndvi_landscape.csv")
for (m in c("surfscape_sr1", "terra_moran")) for (k in c("value1", "value2")) {
  add(paste0(sub("^.*_", "", m), "_", c(value1 = "2013", value2 = "2023")[[k]]),
      dl[[k]][dl$metric == m], f3(dl[[k]][dl$metric == m]), "decomposition_ndvi_landscape.csv")
}
comp_land_zero <- all(dl$configuration[dl$metric %in% c(paste0("surfscape_", comp_b),
  comp_nd)] == 0)
add("land_comp_config_zero", comp_land_zero, if (comp_land_zero) "yes" else "no",
    "decomposition_ndvi_landscape.csv")
w15 <- tab("decomposition_ndvi_w15_cells.csv")
w15$share <- w15$mad_config / (w15$mad_global + w15$mad_local + w15$mad_config)
for (i in seq_len(nrow(w15))) {
  add(paste0("w15_share_", w15$metric[i]), w15$share[i], f1(100 * w15$share[i]),
      "decomposition_ndvi_w15_cells.csv")
  add(paste0("w15_localglobal_", w15$metric[i]), w15$mad_local[i] / w15$mad_global[i],
      s2(w15$mad_local[i] / w15$mad_global[i]), "decomposition_ndvi_w15_cells.csv")
}
le <- w15$metric[w15$mad_local <= w15$mad_global]
add("w15_local_le_global", paste(le, collapse = ";"), if (length(le)) paste(le, collapse = ", ") else "none",
    "decomposition_ndvi_w15_cells.csv")
add("w15_max_residual", max(w15$max_residual), sci(max(w15$max_residual)),
    "decomposition_ndvi_w15_cells.csv")
add("w15_gcor_na_pct", 100 * (1 - w15$n[w15$metric == "gcor"] / w15$n[w15$metric == "sq"]),
    f1(100 * (1 - w15$n[w15$metric == "gcor"] / w15$n[w15$metric == "sq"])),
    "decomposition_ndvi_w15_cells.csv")
sp <- tab("decomposition_ndvi_scale_profile.csv")
for (m in c("sr1", "gcon", "slac", "sdq", "gcor")) for (w in c(5, 9, 15, 31, 61)) {
  x <- sp[sp$metric == m & sp$w == w, ]
  add(paste("scale_share", m, w, sep = "__"), x$share_config, f1(100 * x$share_config),
      "decomposition_ndvi_scale_profile.csv")
  add(paste("scale_sign", m, w, sep = "__"), x$sign_disagree, f2(x$sign_disagree),
      "decomposition_ndvi_scale_profile.csv")
}
comp_scale_zero <- all(sp$mad_config[sp$metric %in% comp_b] <= 1e-12 * sp$mad_total[sp$metric %in% comp_b])
add("scale_comp_config_zero", comp_scale_zero, if (comp_scale_zero) "yes" else "no",
    "decomposition_ndvi_scale_profile.csv")
zm <- tab("zoom_area.csv")
add("zoom_top", zm$top_row, as.character(zm$top_row), "zoom_area.csv")
add("zoom_left", zm$left_col, as.character(zm$left_col), "zoom_area.csv")

out <- do.call(rbind, items)
stopifnot(!anyDuplicated(out$name))
write.csv(out, here("output", "tables", "numbers_for_text.csv"), row.names = FALSE)
message(nrow(out), " numbers written")
