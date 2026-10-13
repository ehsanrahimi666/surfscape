## Tables for the manuscript, written as CSV files and as Markdown tables that
## are inserted into the manuscript sources.

suppressPackageStartupMessages({
  library(here)
})

tab_dir <- here("output", "tables")
md_dir <- here("paper", "src", "tables")
dir.create(md_dir, showWarnings = FALSE, recursive = TRUE)

fmt <- function(x, d = 2) {
  out <- ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
  out[!is.na(x) & x == 0] <- "0"
  out
}
fmt_sig <- function(x, d = 2) {
  out <- ifelse(is.na(x), "", formatC(signif(x, d), format = "fg", digits = d, flag = "#"))
  out[!is.na(x) & x == 0] <- "0"
  sub("\\.$", "", out)
}
minus <- function(x) gsub("-", "−", x, fixed = TRUE, useBytes = TRUE)
## relative column widths (number of dashes in the separator line, which pandoc
## uses for tables wider than the text)
col_widths <- list(
  table2.md = c(16, 20, 9, 9, 9, 9, 14, 14), table3.md = c(18, 18, 15, 13, 8, 8, 20),
  table4.md = c(26, 9, 9, 11, 12, 13, 10), tableA1.md = c(24, 14, 32, 24, 6),
  tableA3a.md = c(24, rep(9, 8)), tableA3b.md = c(24, rep(9, 8)),
  tableA4.md = c(26, 12, 11, 11, 11, 14, 15), tableA5.md = c(28, 12, 22, 12, 26),
  tableC1.md = c(10, 24, 11, 11, 11, 11, 11, 11), tableC2.md = c(25, 15, 15, 15, 15, 15),
  tableC3.md = c(20, 6, 6, 12, 12, 12, 12, 14), tableC4.md = c(26, 9, 9, 11, 12, 13, 10),
  tableC5.md = c(20, 9, 9, 11, 11, 11, 13, 9, 7), tableD1.md = c(30, 70),
  tableD2.md = c(18, 10, 52, 12))
md_table <- function(df, file, caption = NULL) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  wd <- col_widths[[file]]
  if (is.null(wd)) wd <- rep(3, ncol(df))
  stopifnot(length(wd) == ncol(df))
  sep <- paste0("|", paste(strrep("-", wd), collapse = "|"), "|")
  body <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  writeLines(c(hdr, sep, body), file.path(md_dir, file))
}

catalogue <- read.csv(file.path(tab_dir, "metric_catalogue.csv"))
family_of <- function(m) {
  pkg <- sub("_.*$", "", m)
  nm <- sub("^[^_]+_", "", m)
  if (pkg == "surfscape") {
    if (nm %in% c("smean", "sa", "sq", "ssk", "sku", "sp", "sv", "sz", "s10z", "scv"))
      return("Amplitude and moments")
    if (nm %in% c("sent", "srao")) return("Information-theoretic")
    if (nm == "sr1") return("Autocorrelation")
    if (nm %in% c("sdq", "sdr")) return("Slope, surface area and summits")
    if (nm == "slac") return("Lacunarity")
    return("GLCM, equal-range grey levels")
  }
  if (pkg == "geodiv") {
    if (nm %in% c("sa", "sq", "ssk", "sku", "s10z", "smean", "sv", "sph"))
      return("Amplitude and moments")
    if (nm %in% c("sbi", "sci", "svi", "sk", "spk", "svk")) return("Bearing area")
    if (nm %in% c("scl20", "str20", "srw", "srwi", "shw"))
      return("Autocorrelation length and spectrum")
    return("Slope, surface area and summits")
  }
  if (pkg == "terra") return("Autocorrelation")
  if (pkg == "rasterdiv") return("Information-theoretic")
  if (grepl("^glcmtex_prob", m)) return("GLCM, equal-probability grey levels")
  "GLCM, equal-range grey levels"
}
fam_levels <- c("Amplitude and moments", "Bearing area", "Information-theoretic",
                "Autocorrelation", "Autocorrelation length and spectrum",
                "Slope, surface area and summits", "Lacunarity",
                "GLCM, equal-range grey levels", "GLCM, equal-probability grey levels")
catalogue$family <- factor(vapply(catalogue$implementation, family_of, ""), fam_levels)
write.csv(catalogue, file.path(tab_dir, "metric_catalogue_families.csv"), row.names = FALSE)

summ <- read.csv(file.path(tab_dir, "benchmark_main_summary.csv"))
summ$family <- factor(vapply(summ$metric, family_of, ""), fam_levels)
resp <- read.csv(file.path(tab_dir, "benchmark_main_responses.csv"))

## ---- Table 2: classification by family ---------------------------------------
count_str <- function(cats) {
  k <- table(factor(cats, c("compositional", "configurational", "mixed", "unresponsive")))
  paste0(k[["compositional"]], " / ", k[["configurational"]], " / ", k[["mixed"]],
         if (k[["unresponsive"]] > 0) paste0(" (+", k[["unresponsive"]], " unresp.)") else "")
}
fx <- summ[summ$protocol == "fixed", ]
own <- summ[summ$protocol == "own", ]
t2 <- do.call(rbind, lapply(fam_levels, function(f) {
  m <- catalogue$implementation[catalogue$family == f]
  pk <- table(sub("_.*$", "", m))
  pk_txt <- paste(paste0(names(pk), " ", as.integer(pk)), collapse = ", ")
  sims <- vapply(c("r03", "r10", "r30"), function(L)
    count_str(fx$category[fx$landscape == L & fx$metric %in% m]), "")
  nd <- fx[fx$landscape == "ndvi" & fx$metric %in% m, ]
  own_nd <- own[own$landscape == "ndvi" & own$metric %in% m &
                  own$metric %in% catalogue$implementation[catalogue$binned], ]
  data.frame(
    Family = f, Implementations = paste0(length(m), " (", pk_txt, ")"),
    `Simulated, range 3` = sims[["r03"]], `Simulated, range 10` = sims[["r10"]],
    `Simulated, range 30` = sims[["r30"]], `NDVI` = count_str(nd$category),
    `NDVI invariant to shift / stretch / reshape` = paste(sum(nd$shift_invariant),
      sum(nd$scale_invariant), sum(nd$reshape_invariant), sep = " / "),
    `NDVI, own-range protocol (binned implementations)` =
      if (nrow(own_nd)) count_str(own_nd$category) else "–",
    check.names = FALSE)
}))
write.csv(t2, file.path(tab_dir, "table2_classification.csv"), row.names = FALSE)
md_table(t2, "table2.md")

## ---- Appendix Table A3: every implementation ---------------------------------
short_cat <- c(compositional = "comp", configurational = "conf", mixed = "mixed",
               unresponsive = "none")
a3 <- do.call(rbind, lapply(catalogue$implementation, function(m) {
  row <- data.frame(Implementation = m, check.names = FALSE)
  for (L in c("r03", "r10", "r30", "ndvi")) {
    s <- fx[fx$landscape == L & fx$metric == m, ]
    lab <- c(r03 = "range 3", r10 = "range 10", r30 = "range 30", ndvi = "NDVI")[[L]]
    row[[paste0(lab, ": class")]] <- short_cat[s$category]
    row[[paste0(lab, ": perm.")]] <- fmt_sig(s$configuration_response)
    row[[paste0(lab, ": comp.")]] <- fmt_sig(s$composition_response)
    row[[paste0(lab, ": signif.")]] <- fmt(s$signif, 2)
  }
  row
}))
write.csv(a3, file.path(tab_dir, "tableA3_all_implementations.csv"), row.names = FALSE)
md_table(a3[, c(1, 2:9)], "tableA3a.md")
md_table(a3[, c(1, 10:17)], "tableA3b.md")

## ---- Appendix Table A4: NDVI responses to each test --------------------------
tests <- c("permute", "shift", "scale", "reshape", "target", "rearrange")
rn <- resp[resp$landscape == "ndvi" & resp$protocol == "fixed", ]
a4 <- do.call(rbind, lapply(catalogue$implementation, function(m) {
  r <- rn[rn$metric == m, ]
  vals <- vapply(tests, function(t) {
    x <- r[r$test == t, ]
    if (!nrow(x)) return("")
    if (isTRUE(x$invariant)) "0" else fmt_sig(x$response)
  }, "")
  data.frame(Implementation = m, t(vals), check.names = FALSE)
}))
names(a4) <- c("Implementation", "Permutation", "Shift", "Stretch", "Reshape",
               "Observed distribution", "Landscape rearrangement")
write.csv(a4, file.path(tab_dir, "tableA4_ndvi_responses.csv"), row.names = FALSE)
md_table(a4, "tableA4.md")

## ---- Appendix Table A5: binning protocols for binned metrics -----------------
binned <- catalogue$implementation[catalogue$binned]
a5 <- do.call(rbind, lapply(binned, function(m) {
  out <- data.frame(Implementation = m, check.names = FALSE)
  for (p in c("fixed", "own")) {
    s <- summ[summ$landscape == "ndvi" & summ$protocol == p & summ$metric == m, ]
    inv <- c(s$shift_invariant, s$scale_invariant, s$reshape_invariant, s$target_invariant)
    out[[paste0(p, ": class")]] <- short_cat[s$category]
    out[[paste0(p, ": invariant to")]] <- {
      lab <- c("shift", "stretch", "reshape", "observed")[inv]
      if (length(lab)) paste(lab, collapse = ", ") else "none"
    }
  }
  out
}))
write.csv(a5, file.path(tab_dir, "tableA5_protocols.csv"), row.names = FALSE)
md_table(a5, "tableA5.md")
print(t2)

## ---- Table 4 and Appendix Table C4: landscape-level decomposition of NDVI ---------
dl <- read.csv(file.path(tab_dir, "decomposition_ndvi_landscape.csv"))
dl$metric <- ifelse(dl$set == "builtin", paste0("surfscape_", dl$metric), dl$metric)
tol0 <- function(x, ref) ifelse(abs(x) <= 1e-9 * pmax(abs(ref), 1e-300), 0, x)
dl$configuration <- tol0(dl$configuration, dl$total)
dl$share_config <- ifelse(dl$configuration == 0, 0, dl$share_config)
fmt_val <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    if (v == 0) return("0")
    if (abs(v) >= 1e-3 && abs(v) < 1e5) {
      out <- trimws(formatC(signif(v, 3), format = "fg", digits = 3))
    } else {
      e <- floor(log10(abs(v)))
      m <- signif(v / 10^e, 3)
      if (abs(m) >= 10) {
        m <- m / 10
        e <- e + 1
      }
      out <- paste0(trimws(formatC(m, format = "fg", digits = 3)), " × 10^", e, "^")
    }
    minus(out)
  }, "", USE.NAMES = FALSE)
}
mk_row <- function(r) data.frame(
  Metric = r$metric, `2013` = fmt_val(r$value1), `2023` = fmt_val(r$value2),
  `Total change` = fmt_val(r$total), Compositional = fmt_val(r$composition),
  Configurational = fmt_val(r$configuration),
  `Configuration (%)` = formatC(100 * r$share_config, format = "f", digits = 1),
  check.names = FALSE)
c3 <- do.call(rbind, lapply(seq_len(nrow(dl)), function(i) cbind(Set = dl$set[i], mk_row(dl[i, ]))))
write.csv(c3, file.path(tab_dir, "tableC4_landscape_decomposition_all.csv"), row.names = FALSE)
md_table(c3[, -1], "tableC4.md")
main_rows <- c("surfscape_smean", "surfscape_sq", "surfscape_ssk", "surfscape_sku",
               "surfscape_sent", "surfscape_srao", "surfscape_sr1", "terra_moran",
               "terra_geary", "surfscape_sdq", "surfscape_slac", "surfscape_gcon",
               "surfscape_gcor", "glcmtexown_prob_contrast", "glcmtexown_prob_correlation",
               "geodiv_sds", "geodiv_scl20")
t4 <- do.call(rbind, lapply(main_rows, function(m) mk_row(dl[dl$metric == m, ][1, ])))
write.csv(t4, file.path(tab_dir, "table4_landscape_decomposition.csv"), row.names = FALSE)
md_table(t4, "table4.md")

## ---- Appendix Table C3: controlled test of the decomposition ------------------------
dv <- read.csv(file.path(tab_dir, "decomposition_virtual_windows.csv"))
agg <- aggregate(cbind(rank_cor, mad_total, mad_global, mad_local, mad_config) ~ kappa + lambda +
                   metric, data = dv, FUN = mean)
c2 <- agg[agg$metric %in% c("sq", "sent", "sr1", "gcon", "slac"), ]
c2 <- c2[order(match(c2$metric, c("sq", "sent", "sr1", "gcon", "slac")), c2$kappa, c2$lambda), ]
c2_out <- data.frame(Metric = paste0("surfscape_", c2$metric), "κ" = c2$kappa, "λ" = c2$lambda,
                     `Rank correlation` = formatC(c2$rank_cor, format = "f", digits = 2),
                     Total = fmt_val(c2$mad_total), Global = fmt_val(c2$mad_global),
                     Local = fmt_val(c2$mad_local),
                     Configurational = fmt_val(ifelse(c2$mad_config <= 1e-9 * c2$mad_total, 0,
                                                      c2$mad_config)),
                     check.names = FALSE)
write.csv(c2_out, file.path(tab_dir, "tableC3_decomposition_test.csv"), row.names = FALSE)
md_table(c2_out, "tableC3.md")

## ---- Appendix Tables C1 and C2: extent and grain ------------------------------
eg <- read.csv(file.path(tab_dir, "benchmark_extent_grain_summary.csv"))
eg <- eg[eg$protocol == "fixed", ]
c1_metrics <- c("surfscape_sq", "surfscape_sent", "surfscape_sr1", "terra_moran",
                "surfscape_gcon", "glcmtex_range_correlation", "surfscape_slac", "surfscape_sdr")
lev <- c(w5 = "5", w9 = "9", w15 = "15", w31 = "31", w61 = "61", landscape = "L")
c1 <- do.call(rbind, lapply(c("r03", "r10", "r30", "ndvi30"), function(g) {
  do.call(rbind, lapply(c1_metrics, function(m) {
    s <- eg[eg$group == g & eg$metric == m, ]
    s <- s[match(names(lev), s$level), ]
    cells <- ifelse(s$category == "compositional", "inv.",
                    paste0(fmt_sig(s$configuration_response), " (", fmt(s$signif, 2), ")"))
    cells[s$level == "landscape" & s$category != "compositional"] <-
      paste0("sig. ", fmt(s$signif[s$level == "landscape" & s$category != "compositional"], 1))
    data.frame(Landscape = c(r03 = "range 3", r10 = "range 10", r30 = "range 30",
                             ndvi30 = "NDVI")[[g]], Metric = m, t(stats::setNames(cells, lev)),
               check.names = FALSE)
  }))
}))
write.csv(c1, file.path(tab_dir, "tableC1_extent.csv"), row.names = FALSE)
md_table(c1, "tableC1.md")

gr <- do.call(rbind, lapply(c1_metrics, function(m) {
  pick <- function(g, l) {
    s <- eg[eg$group == g & eg$level == l & eg$metric == m, ]
    if (!nrow(s)) return("")
    if (s$category == "compositional") "inv." else
      paste0(fmt_sig(s$configuration_response), " (", fmt(s$signif, 2), ")")
  }
  data.frame(Metric = m,
             `30 m, 15 cells` = pick("ndvi30", "w15"), `60 m, 15 cells` = pick("ndvi60", "w15"),
             `120 m, 15 cells` = pick("ndvi120", "w15"), `30 m, 31 cells` = pick("ndvi30", "w31"),
             `120 m, 7 cells` = pick("ndvi120", "w7"),
             check.names = FALSE)
}))
write.csv(gr, file.path(tab_dir, "tableC2_grain.csv"), row.names = FALSE)
md_table(gr, "tableC2.md")

## ---- Appendix Table A1: metric implementations ---------------------------------
a1 <- data.frame(Implementation = catalogue$implementation,
                 Package = catalogue$package,
                 Definition = catalogue$definition,
                 Settings = catalogue$settings,
                 Binned = ifelse(catalogue$binned, "yes", "no"),
                 check.names = FALSE)
md_table(a1, "tableA1.md")

## ---- Appendix Table A2: NDVI at the three grains -------------------------------
aws <- read.csv(file.path(tab_dir, "analysis_window_summary.csv"))
aws$year <- ifelse(grepl("2013", aws$file), "2013", "2023")
aws <- aws[order(aws$res, aws$year), ]
a2 <- data.frame(Grain = paste(aws$res, "m"), Year = aws$year,
                 Cells = paste(aws$nrow, "×", aws$ncol),
                 Mean = formatC(aws$mean, format = "f", digits = 3),
                 SD = formatC(aws$sd, format = "f", digits = 3),
                 Skewness = formatC(aws$skewness, format = "f", digits = 2),
                 Minimum = minus(formatC(aws$min, format = "f", digits = 3)),
                 Maximum = formatC(aws$max, format = "f", digits = 3),
                 check.names = FALSE)
write.csv(a2, file.path(tab_dir, "tableA2_ndvi_grains.csv"), row.names = FALSE)
md_table(a2, "tableA2.md")

## ---- Table 3: key subset of metrics ---------------------------------------------
sel <- read.csv(file.path(tab_dir, "subset_selected.csv"))
role_lab <- c(location = "Composition: location", spread = "Composition: spread",
              shape = "Composition: shape",
              diversity = "Compositional heterogeneity: diversity of values",
              configuration_autocorrelation = "Configuration: autocorrelation",
              configuration_texture = "Configuration: texture")
bin_set <- catalogue$implementation[catalogue$binned]
own_cat <- ifelse(sel$metric %in% bin_set,
                  paste0(sel$category_ndvi_own, " / ", sel$category_r10_own), "–")
red_txt <- function(a, b) {
  x <- c(if (nzchar(a)) paste0(gsub(";", ", ", a), " (NDVI)"),
         if (nzchar(b)) paste0(gsub(";", ", ", b), " (simulated)"))
  if (length(x)) paste(x, collapse = "; ") else "none"
}
t3 <- data.frame(Role = role_lab[sel$role], Implementation = sel$metric,
                 `Category, fixed range (NDVI / simulated, range 10)` =
                   paste0(sel$category_ndvi, " / ", sel$category_r10),
                 `Category, own range (NDVI / simulated, range 10)` = own_cat,
                 `Largest abs. rho, NDVI` = formatC(sel$max_abs_rho_ndvi, format = "f", digits = 2),
                 `Largest abs. rho, simulated` = formatC(sel$max_abs_rho_sim, format = "f", digits = 2),
                 `Redundant (abs. rho >= 0.7) with` = mapply(red_txt,
                   ifelse(is.na(sel$redundant_with_ndvi), "", sel$redundant_with_ndvi),
                   ifelse(is.na(sel$redundant_with_sim), "", sel$redundant_with_sim)),
                 check.names = FALSE)
slog <- read.csv(file.path(tab_dir, "subset_selection_log.csv"))
empty <- setdiff(names(role_lab), sel$role)
for (r in empty) {
  cands <- slog$candidate[slog$role == r]
  t3 <- rbind(t3, data.frame(Role = role_lab[[r]], Implementation = "none selected",
    `Category, fixed range (NDVI / simulated, range 10)` = "",
    `Category, own range (NDVI / simulated, range 10)` = "",
    `Largest abs. rho, NDVI` = "", `Largest abs. rho, simulated` = "",
    `Redundant (abs. rho >= 0.7) with` = paste0("candidates ", paste(cands, collapse = ", "),
      " redundant in both landscape types"), check.names = FALSE))
}
t3 <- t3[order(match(t3$Role, role_lab)), ]
write.csv(t3, file.path(tab_dir, "table3_subset.csv"), row.names = FALSE)
md_table(t3, "table3.md")

## ---- Appendix Table C5: decomposition across window sizes ------------------------
sp <- read.csv(file.path(tab_dir, "decomposition_ndvi_scale_profile.csv"))
noncomp <- c("sr1", "sdq", "sdr", "slac", "gcon", "gdis", "ghom", "gasm", "gent", "gcor")
sp <- sp[sp$metric %in% noncomp, ]
sp <- sp[order(match(sp$metric, noncomp), sp$w), ]
c5 <- data.frame(Metric = paste0("surfscape_", sp$metric), `Window (cells)` = sp$w,
                 Cells = sp$n,
                 Total = fmt_val(sp$mad_total), Global = fmt_val(sp$mad_global),
                 Local = fmt_val(sp$mad_local), Configurational = fmt_val(sp$mad_config),
                 `Share (%)` = formatC(100 * sp$share_config, format = "f", digits = 1),
                 Sign = formatC(sp$sign_disagree, format = "f", digits = 2),
                 check.names = FALSE)
write.csv(c5, file.path(tab_dir, "tableC5_scale_profile.csv"), row.names = FALSE)
md_table(c5, "tableC5.md")

## ---- Appendix Table D1: functions of the package ----------------------------------
d1 <- data.frame(
  Function = c("ss_benchmark()", "ss_benchmark_combine()", "ss_perturb()", "ss_simulate()",
               "ss_metrics()", "ss_metric_list()", "ss_rank_transfer()", "ss_decompose()",
               "ss_decompose_landscape()", "ss_scale_profile()", "ss_extract()"),
  Purpose = c("responses of any metric to the composition and configuration tests, and its classification",
              "pools benchmark results computed on separate windows or landscapes",
              "applies one of the controlled changes to a whole surface",
              "simulates Gaussian random fields with a given practical range",
              "moving-window maps of the built-in metrics",
              "catalogue of the built-in metrics",
              "one surface's spatial ordering with another surface's values",
              "three-component decomposition of change in moving windows",
              "two-component decomposition of change in landscape-level metrics",
              "the decomposition across window sizes",
              "decomposed change in buffers around points"))
md_table(d1, "tableD1.md")

## ---- Appendix Table D2: errors and warnings of other packages ------------------------
lg <- read.csv(here("logs", "benchmark_main_external_log.csv"))
lg <- lg[order(lg$metric, lg$type), ]
d2 <- data.frame(Function = lg$metric, Type = lg$type, Message = lg$message,
                 Occurrences = lg$n, check.names = FALSE)
md_table(d2, "tableD2.md")

## ---- Appendix Table D3: software versions ----------------------------------------------
pk <- c("surfscape", "terra", "geodiv", "GLCMTextures", "glcm", "rasterdiv", "Rcpp", "here")
d3 <- data.frame(Software = c("R", pk),
                 Version = c(paste(R.version$major, R.version$minor, sep = "."),
                             vapply(pk, function(p) utils::packageDescription(p)$Version, "")),
                 check.names = FALSE)
write.csv(d3, file.path(tab_dir, "tableD3_versions.csv"), row.names = FALSE)
md_table(d3, "tableD3.md")
