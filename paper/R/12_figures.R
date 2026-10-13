## Figures for the revised manuscript, drawn from the saved outputs.

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(surfscape)
})

fig_dir <- here("output", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
div <- hcl.colors(101, "Blue-Red 3")
ndvi_col <- hcl.colors(100, "YlGn", rev = TRUE)
cat_col <- c(compositional = "#1b9e77", configurational = "#d95f02", mixed = "#7570b3",
             unresponsive = "grey60")
png_open <- function(name, w, h) png(file.path(fig_dir, name), width = w, height = h, res = 300)
options(scipen = 10)

args <- commandArgs(trailingOnly = TRUE)
which_figs <- if (length(args)) args else c("fig1", "fig2", "fig3", "fig4", "fig5",
                                            "fig6", "fig7", "fig8")

## ---- Figure 2: study area -------------------------------------------------
if ("fig2" %in% which_figs) {
  r13 <- rast(here("data", "raw", "NDVI-2013.tif"))
  r23 <- rast(here("data", "raw", "NDVI-2023.tif"))
  aw <- read.csv(here("output", "tables", "analysis_window.csv"))
  box <- vect(ext(aw$xmin, aw$xmax, aw$ymin, aw$ymax), crs = crs(r13))
  d <- r23 - r13
  rng <- c(-0.1, 0.6)
  dl <- 0.3
  png_open("Fig2_study_area.png", 2400, 760)
  par(mfrow = c(1, 3), mar = c(1.8, 2.6, 1.4, 4.4), mgp = c(1.2, 0.3, 0), tcl = -0.2,
      cex.axis = 0.6)
  plot(clamp(r13, rng[1], rng[2]), range = rng, col = ndvi_col, main = "(a) NDVI 2013",
       cex.main = 0.95, plg = list(cex = 0.7))
  lines(box, lwd = 1.5)
  plot(clamp(r23, rng[1], rng[2]), range = rng, col = ndvi_col, main = "(b) NDVI 2023",
       cex.main = 0.95, plg = list(cex = 0.7))
  lines(box, lwd = 1.5)
  plot(clamp(d, -dl, dl), range = c(-dl, dl), col = div, main = "(c) 2023 minus 2013",
       cex.main = 0.95, plg = list(cex = 0.7))
  lines(box, lwd = 1.5)
  invisible(dev.off())
}

## ---- Figure 1: the tests ----------------------------------------------------
if ("fig1" %in% which_figs) {
  sim <- read.csv(here("output", "tables", "simulated_landscapes.csv"))
  s <- sim[sim$range == 10 & sim$rep == 1, ]
  z <- rast(here(s$file_landscape))
  a <- rast(here(s$file_arrangement))
  pan <- list("(a) original" = z,
              "(b) shift" = ss_perturb(z, "shift", 0.5),
              "(c) stretch" = ss_perturb(z, "scale", 1.5),
              "(d) reshape" = ss_perturb(z, "reshape", 2),
              "(e) landscape-wide\nrearrangement" = ss_perturb(z, "rearrange", arrangement = a))
  rng <- range(vapply(pan, function(r) range(values(r)[, 1]), numeric(2)))
  m <- as.matrix(z, wide = TRUE)
  win <- m[94:108, 94:108]
  set.seed(2)
  pw <- matrix(sample(as.vector(win)), 15)
  wr <- range(win)
  png_open("Fig1_tests.png", 2400, 1050)
  layout(matrix(c(1:5, 6, 7, 8, 8, 8), 2, byrow = TRUE), heights = c(1, 0.9))
  par(mar = c(0.4, 0.4, 2.4, 0.4))
  for (nm in names(pan)) {
    image(t(as.matrix(pan[[nm]], wide = TRUE))[, 201:1], col = hcl.colors(100, "viridis"),
          zlim = rng, axes = FALSE, asp = 1, main = nm, cex.main = 0.85)
    if (nm == "(a) original") rect(93.5 / 201, 1 - 108.5 / 201, 108.5 / 201, 1 - 93.5 / 201,
                                   border = "white", lwd = 1.5)
  }
  image(t(win)[, 15:1], col = hcl.colors(100, "viridis"), zlim = wr, axes = FALSE, asp = 1,
        main = "(f) 15 \u00d7 15 window", cex.main = 0.85)
  image(t(pw)[, 15:1], col = hcl.colors(100, "viridis"), zlim = wr, axes = FALSE, asp = 1,
        main = "(g) same values, permuted", cex.main = 0.85)
  par(mar = c(2.8, 3, 2.4, 1), mgp = c(1.6, 0.4, 0), tcl = -0.2)
  v <- values(z)[, 1]
  dens <- lapply(pan[1:4], function(r) density(values(r)[, 1]))
  plot(NA, xlim = range(vapply(dens, function(d) range(d$x), numeric(2))),
       ylim = c(0, max(vapply(dens, function(d) max(d$y), numeric(1)))),
       xlab = "value", ylab = "density", main = "(h) value distributions", cex.main = 0.85,
       cex.lab = 0.8, cex.axis = 0.7)
  lc <- c("black", "#e7298a", "#1b9e77", "#d95f02")
  for (i in seq_along(dens)) lines(dens[[i]], col = lc[i], lwd = 1.4)
  legend("topright", c("original (and permuted, rearranged)", "shift", "stretch", "reshape"),
         col = lc, lwd = 1.4, bty = "n", cex = 0.7)
  invisible(dev.off())
}

## ---- shared helpers for the benchmark figures --------------------------------
if (any(c("fig3", "fig4", "fig5", "fig6") %in% which_figs)) {
  catalogue <- read.csv(here("output", "tables", "metric_catalogue.csv"))
  short <- function(m) {
    m <- sub("^surfscape_", "ss: ", m)
    m <- sub("^geodiv_", "gd: ", m)
    m <- sub("^terra_", "terra: ", m)
    m <- sub("^glcmtex_range_", "GT-r: ", m)
    m <- sub("^glcmtex_prob_", "GT-p: ", m)
    m <- sub("^glcmpkg_", "glcm: ", m)
    sub("^rasterdiv_", "rd: ", m)
  }
}

## ---- Figure 3: responses of all implementations --------------------------------
if ("fig3" %in% which_figs) {
  resp <- read.csv(here("output", "tables", "benchmark_main_responses.csv"))
  summ <- read.csv(here("output", "tables", "benchmark_main_summary.csv"))
  resp <- resp[resp$protocol == "fixed", ]
  summ <- summ[summ$protocol == "fixed", ]
  mets <- catalogue$implementation
  land <- c(r10 = "Simulated, range 10 cells", ndvi = "NDVI 2013")
  tests <- list(r10 = c("permute", "shift", "scale", "reshape", "rearrange"),
                ndvi = c("permute", "shift", "scale", "reshape", "target", "rearrange"))
  tlab <- c(permute = "permutation", shift = "shift", scale = "stretch", reshape = "reshape",
            target = "observed\ndistribution", rearrange = "landscape\nrearrangement")
  brks <- seq(-3, 1.5, length.out = 101)
  pal <- hcl.colors(100, "Viridis")
  png_open("Fig3_benchmark_heatmap.png", 2400, 3300)
  layout(matrix(c(1, 2, 3, 4, 5, 5), 3, 2, byrow = TRUE), widths = c(1, 1),
         heights = c(1, 0.02, 0.06))
  for (L in names(land)) {
    par(mar = c(5.5, if (L == "r10") 8.5 else 1, 2.5, 3.2))
    tt <- tests[[L]]
    M <- matrix(NA_real_, length(mets), length(tt), dimnames = list(mets, tt))
    I <- M
    for (j in seq_along(tt)) {
      r <- resp[resp$landscape == L & resp$test == tt[j], ]
      i <- match(mets, r$metric)
      M[, j] <- log10(pmax(r$response[i], 1e-3))
      I[, j] <- r$invariant[i]
    }
    Z <- pmin(pmax(M, -3), 1.5)
    Z[I == TRUE] <- NA
    image(seq_along(tt), seq_along(mets), t(Z)[, length(mets):1, drop = FALSE], col = pal,
          breaks = brks, axes = FALSE, xlab = "", ylab = "", main = land[[L]], cex.main = 1)
    inv <- which(t(I)[, length(mets):1] == TRUE, arr.ind = TRUE)
    if (nrow(inv)) rect(inv[, 1] - 0.5, inv[, 2] - 0.5, inv[, 1] + 0.5, inv[, 2] + 0.5,
                        col = "grey92", border = NA)
    abline(h = seq(0.5, length(mets) + 0.5), col = "white", lwd = 0.3)
    axis(1, at = seq_along(tt), labels = tlab[tt], las = 2, cex.axis = 0.75, tick = FALSE)
    if (L == "r10") axis(2, at = rev(seq_along(mets)), labels = short(mets), las = 1,
                         cex.axis = 0.5, tick = FALSE, line = -0.6)
    cats <- summ$category[match(mets, summ$metric[summ$landscape == L])]
    cats <- summ[summ$landscape == L, ]$category[match(mets, summ[summ$landscape == L, ]$metric)]
    rect(length(tt) + 0.6, rev(seq_along(mets)) - 0.5, length(tt) + 1.1,
         rev(seq_along(mets)) + 0.5, col = cat_col[cats], border = NA, xpd = NA)
    groups <- cumsum(table(factor(sub("_.*$", "", mets), levels = unique(sub("_.*$", "", mets)))))
    abline(h = length(mets) - groups[-length(groups)] + 0.5, lwd = 1.2)
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  plot.new()
  par(mar = c(2.5, 10, 0.5, 10))
  plot.new()
  plot.window(c(-3, 1.5), c(0, 1))
  rect(brks[-101], 0.35, brks[-1], 1, col = pal, border = NA)
  axis(1, at = c(-3, -2, -1, 0, 1), labels = c("<= 0.001", "0.01", "0.1", "1", "10"),
       cex.axis = 0.7, mgp = c(1, 0.3, 0))
  rect(1.7, 0.35, 1.95, 1, col = "grey92", border = "grey60", xpd = NA)
  text(2.0, 0.67, "invariant", adj = 0, cex = 0.7, xpd = NA)
  rect(1.7, -0.45, 1.95, 0.2, col = "white", border = "grey60", xpd = NA)
  text(2.0, -0.12, "undefined", adj = 0, cex = 0.7, xpd = NA)
  legend(-4.35, 1.25, legend = names(cat_col)[1:3], fill = cat_col[1:3], border = NA,
         bty = "n", cex = 0.7, xpd = NA, title = "category", title.adj = 0)
  mtext("standardised response (SD among windows)", side = 1, line = 1.4, cex = 0.6)
  invisible(dev.off())
}

## ---- Figure 4: composition against configuration response ------------------
if ("fig4" %in% which_figs) {
  summ <- read.csv(here("output", "tables", "benchmark_main_summary.csv"))
  summ <- summ[summ$protocol == "fixed", ]
  pk_pch <- c(surfscape = 21, geodiv = 22, terra = 24, glcmtex = 23, glcmpkg = 25,
              rasterdiv = 4)
  lab_pos <- c(surfscape_sr1 = 4, terra_moran = 1, surfscape_gcon = 2,
               glcmtex_range_mean = 4, surfscape_slac = 4, surfscape_sdq = 1,
               geodiv_scl20 = 4, geodiv_sdq6 = 2)
  png_open("Fig4_sensitivity_plane.png", 2400, 1150)
  par(mfrow = c(1, 2), mar = c(3.4, 4.2, 2, 0.8), mgp = c(2.4, 0.5, 0), tcl = -0.2)
  for (L in c("r10", "ndvi")) {
    s <- summ[summ$landscape == L, ]
    k <- sub("_.*$", "", s$metric)
    comp_tests <- intersect(c("shift_invariant", "scale_invariant", "reshape_invariant",
                              "target_invariant"), names(s))
    comp_tests <- comp_tests[colSums(!is.na(s[comp_tests])) > 0]
    comp_inv <- apply(as.matrix(s[comp_tests]), 1, all)
    xv <- ifelse(comp_inv, 0, s$composition_response)
    yv <- ifelse(s$permute_invariant, 0, s$configuration_response)
    keep <- is.finite(xv) & is.finite(yv)
    s <- s[keep, ]; k <- k[keep]; xv <- xv[keep]; yv <- yv[keep]
    pos <- c(xv, yv)
    fl <- 10^floor(log10(min(pos[pos > 0])) - 0.5)
    px <- ifelse(xv == 0, fl, xv)
    py <- ifelse(yv == 0, fl, yv)
    plot(px, py, log = "xy", type = "n", axes = FALSE,
         xlab = "Largest response to a composition test (SD)",
         ylab = "Response to spatial permutation (SD)",
         main = c(r10 = "(a) Simulated landscapes, range 10 cells", ndvi = "(b) NDVI 2013")[L],
         cex.main = 0.9, cex.lab = 0.85)
    box()
    for (ax in 1:2) {
      at <- axTicks(ax, log = TRUE)
      at <- at[at > fl * 1.5]
      axis(ax, at = at, labels = format(at, drop0trailing = TRUE, scientific = FALSE),
           cex.axis = 0.7, las = if (ax == 2) 1 else 0)
    }
    abline(h = fl, v = fl, lty = 3, col = "grey60")
    usr <- 10^par("usr")
    if (any(py == fl)) text(usr[1], fl, "invariant", adj = c(-0.05, -0.5), cex = 0.6,
                            col = "grey30")
    if (any(px == fl)) text(fl, usr[4], "invariant", srt = 90, adj = c(1.05, -0.5), cex = 0.6,
                            col = "grey30")
    points(px, py, pch = pk_pch[k], bg = adjustcolor(cat_col[s$category], 0.85),
           col = ifelse(pk_pch[k] == 4, cat_col[s$category], "grey20"), cex = 1.05, lwd = 0.7)
    sel <- s$metric %in% names(lab_pos)
    text(px[sel], py[sel], short(s$metric[sel]), pos = lab_pos[s$metric[sel]], cex = 0.55,
         xpd = NA, offset = 0.35)
    if (L == "r10") {
      legend("bottomleft", legend = names(cat_col)[1:3], pt.bg = cat_col[1:3], pch = 21,
             bty = "n", cex = 0.65, title = "category", title.cex = 0.7, inset = c(0.02, 0.08))
    } else {
      legend("bottomleft", legend = c("surfscape", "geodiv", "terra", "GLCMTextures", "glcm",
                                      "rasterdiv"),
             pch = pk_pch, pt.bg = "white", bty = "n", cex = 0.65, title = "package",
             title.cex = 0.7, inset = c(0.02, 0.08))
    }
  }
  invisible(dev.off())
}

## ---- Figure 6: redundancy among metrics ---------------------------------------
if ("fig6" %in% which_figs) {
  hc <- readRDS(here("output", "rds", "redundancy_hclust_ndvi.rds"))
  cl <- read.csv(here("output", "tables", "redundancy_clusters.csv"))
  sel <- read.csv(here("output", "tables", "subset_selected.csv"))
  labs <- hc$labels[hc$order]
  cats <- cl$category_ndvi[match(labs, cl$metric)]
  n <- length(labs)
  png_open("Fig6_redundancy.png", 2400, 1150)
  par(mar = c(5.6, 3.6, 1.2, 0.6), mgp = c(2, 0.5, 0), tcl = -0.2)
  plot(hc, labels = FALSE, hang = -1, main = "", xlab = "", sub = "",
       ylab = expression(1 - "|Spearman correlation|"), cex.axis = 0.7, cex.lab = 0.8,
       axes = TRUE)
  abline(h = 0.3, lty = 2, col = "grey40")
  y0 <- par("usr")[3]
  is_sel <- labs %in% sel$metric
  text(seq_len(n), y0 - 0.02, short(labs), srt = 90, adj = 1, xpd = NA,
       cex = ifelse(is_sel, 0.5, 0.42), font = ifelse(is_sel, 2, 1), col = cat_col[cats])
  points(which(is_sel), rep(y0 + 0.015, sum(is_sel)), pch = 17, cex = 0.7, xpd = NA)
  present <- intersect(names(cat_col), cats)
  leg <- c(as.list(c(present, "key subset")), list(quote("|" * rho * "| = 0.7")))
  legend("topright", legend = as.expression(leg),
         col = c(cat_col[present], "black", "grey40"),
         pch = c(rep(15, length(present)), 17, NA), lty = c(rep(NA, length(present) + 1), 2),
         bty = "n", cex = 0.65)
  invisible(dev.off())
}

## ---- Figure 7: controlled test of the decomposition ---------------------------
if ("fig7" %in% which_figs) {
  dv <- read.csv(here("output", "tables", "decomposition_virtual_windows.csv"))
  mets <- c(sr1 = "lag-1 autocorrelation", gcon = "GLCM contrast", sq = "Sq")
  comps <- c(mad_global = "global composition", mad_local = "local composition",
             mad_config = "configuration")
  kcol <- c("0" = "grey35", "1" = "#e6ab02", "2" = "#e7298a")
  agg <- aggregate(cbind(mad_global, mad_local, mad_config, mad_total) ~ kappa + lambda + metric,
                   data = dv, FUN = mean)
  png_open("Fig7_decomposition_test.png", 2400, 2000)
  par(mfrow = c(3, 3), mar = c(3, 3.4, 1.8, 0.6), oma = c(0.5, 1.5, 0.5, 0), mgp = c(1.8, 0.45, 0),
      tcl = -0.2)
  for (m in names(mets)) {
    a <- agg[agg$metric == m, ]
    yl <- c(0, max(unlist(a[, names(comps)])) * 1.05)
    for (cp in names(comps)) {
      plot(NA, xlim = c(0, 1), ylim = yl, xlab = expression("change in arrangement (" * lambda * ")"),
           ylab = "mean absolute component", cex.lab = 0.8, cex.axis = 0.7,
           main = paste0(mets[[m]], ": ", comps[[cp]]), cex.main = 0.8)
      for (k in names(kcol)) {
        s <- a[a$kappa == as.numeric(k), ]
        s <- s[order(s$lambda), ]
        lines(s$lambda, s[[cp]], type = "b", pch = 19, cex = 0.7, col = kcol[[k]], lwd = 1.3)
      }
      if (m == "sr1" && cp == "mad_global") {
        legend("topleft", legend = c(expression("composition change, " * kappa == 0),
                                     expression("composition change, " * kappa == 1),
                                     expression("composition change, " * kappa == 2)),
               col = kcol, pch = 19, lwd = 1.3, bty = "n", cex = 0.7)
      }
    }
  }
  invisible(dev.off())
}

## ---- Figure 8: zoomed decomposition of the NDVI change ------------------------
if ("fig8" %in% which_figs) {
  z13 <- rast(here("data", "derived", "ndvi2013_aw30.tif"))
  z23 <- rast(here("data", "derived", "ndvi2023_aw30.tif"))
  d <- rast(here("output", "rasters", "decomposition_ndvi_w15.tif"))
  zm <- read.csv(here("output", "tables", "zoom_area.csv"))
  e <- ext(zm$xmin, zm$xmax, zm$ymin, zm$ymax)
  a13 <- crop(z13, e)
  a23 <- crop(z23, e)
  dz <- crop(d, e)
  nd_rng <- range(c(values(a13), values(a23)), na.rm = TRUE)
  sym <- function(r) {
    mx <- max(abs(quantile(values(r), c(0.01, 0.99), na.rm = TRUE)))
    if (mx == 0) mx <- 1e-12
    mx
  }
  panel <- function(r, rng, col, main) {
    plot(r, range = rng, col = col, axes = FALSE, main = main, cex.main = 0.75, mar = NA,
         plg = list(cex = 0.55), legend = TRUE, type = "continuous")
  }
  png_open("Fig8_zoom_decomposition.png", 2400, 1080)
  par(mfrow = c(2, 5), mar = c(0.4, 0.4, 1.6, 3.2), oma = c(0, 0, 0, 0))
  m1 <- sym(dz[[c("sr1_total", "sr1_global", "sr1_local", "sr1_config")]])
  m2 <- sym(dz[[c("sq_total", "sq_global", "sq_local", "sq_config")]])
  panel(a13, nd_rng, ndvi_col, "(a) NDVI 2013")
  panel(clamp(dz[["sr1_total"]], -m1, m1), c(-m1, m1) * 1.001, div, "(b) lag-1 autocorr.: total")
  panel(clamp(dz[["sr1_global"]], -m1, m1), c(-m1, m1) * 1.001, div, "(c) global composition")
  panel(clamp(dz[["sr1_local"]], -m1, m1), c(-m1, m1) * 1.001, div, "(d) local composition")
  panel(clamp(dz[["sr1_config"]], -m1, m1), c(-m1, m1) * 1.001, div, "(e) configuration")
  panel(a23, nd_rng, ndvi_col, "(f) NDVI 2023")
  panel(clamp(dz[["sq_total"]], -m2, m2), c(-m2, m2) * 1.001, div, "(g) Sq: total")
  panel(clamp(dz[["sq_global"]], -m2, m2), c(-m2, m2) * 1.001, div, "(h) global composition")
  panel(clamp(dz[["sq_local"]], -m2, m2), c(-m2, m2) * 1.001, div, "(i) local composition")
  panel(clamp(dz[["sq_config"]], -m2, m2), c(-m2, m2) * 1.001, div, "(j) configuration")
  invisible(dev.off())
}

## ---- Figure 5: extent and grain ---------------------------------------------
if ("fig5" %in% which_figs) {
  eg <- read.csv(here("output", "tables", "benchmark_extent_grain_summary.csv"))
  er <- read.csv(here("output", "tables", "benchmark_extent_grain_responses.csv"))
  eg <- eg[eg$protocol == "fixed", ]
  er <- er[er$protocol == "fixed", ]
  xs <- c(w5 = 5, w9 = 9, w15 = 15, w31 = 31, w61 = 61, landscape = 160)
  gcol <- c(r03 = "#a6cee3", r10 = "#1f78b4", r30 = "#08306b", ndvi30 = "#e31a1c")
  glab <- c(r03 = "simulated, range 3", r10 = "simulated, range 10",
            r30 = "simulated, range 30", ndvi30 = "NDVI 30 m")
  mixed_m <- c(surfscape_sr1 = "lag-1 autocorrelation", terra_moran = "Moran's I (terra)",
               surfscape_gcon = "GLCM contrast", glcmtex_range_correlation = "GLCM correlation (GT)",
               surfscape_slac = "lacunarity", surfscape_sdr = "Sdr")
  comp_m <- c(surfscape_sq = "Sq", surfscape_ssk = "Ssk", surfscape_sent = "Shannon entropy")
  png_open("Fig5_extent_grain.png", 2400, 2300)
  layout(matrix(c(1:6, 7:9, 10:12), 4, 3, byrow = TRUE))
  par(mar = c(3, 3.3, 1.8, 0.5), mgp = c(1.8, 0.45, 0), tcl = -0.2)
  for (m in names(mixed_m)) {
    plot(NA, xlim = c(4.5, 180), ylim = c(0, 1.05), log = "x", xaxt = "n",
         xlab = "window size (cells)", ylab = "proportion of significant windows",
         main = paste0("(", letters[match(m, names(mixed_m))], ") ", mixed_m[[m]]), cex.main = 0.85,
         cex.lab = 0.75, cex.axis = 0.7)
    axis(1, at = xs, labels = c(5, 9, 15, 31, 61, "L"), cex.axis = 0.7)
    for (g in names(gcol)) {
      s <- eg[eg$group == g & eg$metric == m, ]
      s <- s[match(names(xs), s$level), ]
      lines(xs, s$signif, type = "b", pch = 19, cex = 0.7, col = gcol[[g]], lwd = 1.2)
    }
    if (m == "surfscape_sr1") legend("bottomright", legend = glab, col = gcol, pch = 19, lwd = 1.2,
                                     bty = "n", cex = 0.65)
  }
  for (m in names(mixed_m)[1:3]) {
    plot(NA, xlim = c(4.5, 70), ylim = range(c(0.05, 50)), log = "xy", xaxt = "n", yaxt = "n",
         xlab = "window size (cells)", ylab = "response to permutation (SD)",
         main = paste0("(", letters[6 + match(m, names(mixed_m))], ") ", mixed_m[[m]]), cex.main = 0.85,
         cex.lab = 0.75, cex.axis = 0.7)
    for (g in names(gcol)) {
      s <- eg[eg$group == g & eg$metric == m & eg$level != "landscape", ]
      s <- s[match(names(xs)[1:5], s$level), ]
      lines(xs[1:5], s$configuration_response, type = "b", pch = 19, cex = 0.7, col = gcol[[g]],
            lwd = 1.2)
    }
    axis(1, at = xs[1:5], labels = c(5, 9, 15, 31, 61), cex.axis = 0.7)
    yt <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
    axis(2, at = yt, labels = format(yt, drop0trailing = TRUE), cex.axis = 0.7)
  }
  for (m in names(comp_m)) {
    plot(NA, xlim = c(4.5, 180), ylim = c(0, 1.6), log = "x", xaxt = "n",
         xlab = "window size (cells)", ylab = "response to landscape rearrangement (SD)",
         main = paste0("(", letters[9 + match(m, names(comp_m))], ") ", comp_m[[m]]), cex.main = 0.85,
         cex.lab = 0.75, cex.axis = 0.7)
    axis(1, at = xs, labels = c(5, 9, 15, 31, 61, "L"), cex.axis = 0.7)
    for (g in names(gcol)) {
      s <- er[er$group == g & er$metric == m & er$test == "rearrange", ]
      s <- s[match(names(xs), s$level), ]
      y <- ifelse(s$invariant, 0, s$response)
      lines(xs, y, type = "b", pch = 19, cex = 0.7, col = gcol[[g]], lwd = 1.2)
    }
  }
  invisible(dev.off())
}
