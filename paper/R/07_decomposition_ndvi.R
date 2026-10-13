## Decomposition of the 2013-2023 change in the NDVI analysis window: the whole
## window as one landscape, moving windows of 15 cells, and a series of window
## sizes. Run as: Rscript R/07_decomposition_ndvi.R [landscape|windows|scale]
## (all parts when no argument is given). Metrics are processed in groups to
## limit memory use; each metric is computed independently, so the results do
## not depend on the grouping.

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(surfscape)
})
source(here("R", "external_metrics.R"))
dir.create(here("output", "rasters"), showWarnings = FALSE)
terraOptions(todisk = TRUE, datatype = "FLT8S", progress = 0)

parts <- commandArgs(trailingOnly = TRUE)
if (!length(parts)) parts <- c("landscape", "windows", "scale")

z13 <- rast(here("data", "derived", "ndvi2013_aw30.tif"))
z23 <- rast(here("data", "derived", "ndvi2023_aw30.tif"))
metrics <- ss_metric_list()$metric
groups <- split(metrics, ceiling(seq_along(metrics) / 11))
mm <- cbind(minmax(z13), minmax(z23))
rg <- c(min(mm[1, ]), max(mm[2, ]))
v13 <- values(z13)[, 1]
v23 <- values(z23)[, 1]

## ---- landscape level ------------------------------------------------------
if ("landscape" %in% parts) {
  t0 <- proc.time()[3]
  land <- list(
    builtin = ss_decompose_landscape(z13, z23, metrics = metrics),
    terra = ss_decompose_landscape(z13, z23, fun = terra_metrics),
    glcmtex_pooled = ss_decompose_landscape(z13, z23, fun = function(r)
      glcm_texture_metrics(as.matrix(r, wide = TRUE), surface = c(v13, v23))),
    glcmtex_own = ss_decompose_landscape(z13, z23, fun = function(r)
      glcm_texture_metrics(as.matrix(r, wide = TRUE), surface = values(r)[, 1])),
    geodiv = ss_decompose_landscape(z13, z23, fun = function(r) {
      u <- unbinned_metrics(r)
      u[grepl("^geodiv_", names(u))]
    }))
  land$glcmtex_own$metric <- sub("^glcmtex_", "glcmtexown_", land$glcmtex_own$metric)
  land_tab <- do.call(rbind, lapply(names(land), function(s) cbind(set = s, land[[s]])))
  write.csv(land_tab, here("output", "tables", "decomposition_ndvi_landscape.csv"),
            row.names = FALSE)
  write.csv(ext_log(), here("logs", "decomposition_ndvi_landscape_external_log.csv"),
            row.names = FALSE)
  message("landscape level done in ", round(proc.time()[3] - t0), " s")
}

## ---- 15-cell windows: all metrics -----------------------------------------
if ("windows" %in% parts) {
  t0 <- proc.time()[3]
  keep <- c("sq", "sr1", "gcon", "slac", "smean")
  rows <- list()
  kept <- list()
  for (g in groups) {
    d <- ss_decompose(z13, z23, metrics = g, w = 15, range = rg)
    dmean <- if ("smean" %in% g) NULL else
      values(ss_decompose(z13, z23, metrics = "smean", w = 15, range = rg,
                          components = "total"))[, 1]
    dv <- values(d)
    if (is.null(dmean)) dmean <- dv[, "smean_total"]
    for (m in g) {
      cm <- dv[, paste0(m, c("_total", "_global", "_local", "_config"))]
      ok <- stats::complete.cases(cm) & is.finite(dmean)
      cm <- cm[ok, , drop = FALSE]
      cr <- function(a) if (stats::sd(a) == 0) NA_real_ else stats::cor(a, dmean[ok])
      rows[[length(rows) + 1L]] <- data.frame(
        metric = m, n = nrow(cm),
        max_residual = max(abs(cm[, 1] - rowSums(cm[, 2:4]))),
        mad_total = mean(abs(cm[, 1])), mad_global = mean(abs(cm[, 2])),
        mad_local = mean(abs(cm[, 3])), mad_config = mean(abs(cm[, 4])),
        max_abs_config = max(abs(cm[, 4])),
        r_total_dmean = cr(cm[, 1]), r_global_dmean = cr(cm[, 2]),
        r_local_dmean = cr(cm[, 3]), r_config_dmean = cr(cm[, 4]))
    }
    for (m in intersect(keep, g)) {
      kept[[m]] <- d[[paste0(m, c("_total", "_global", "_local", "_config"))]] * 1
    }
    rm(d, dv)
    gc()
  }
  write.csv(do.call(rbind, rows), here("output", "tables", "decomposition_ndvi_w15_cells.csv"),
            row.names = FALSE)
  dk <- rast(unname(kept[keep]))
  names(dk) <- unlist(lapply(kept[keep], names))
  writeRaster(dk, here("output", "rasters", "decomposition_ndvi_w15.tif"), overwrite = TRUE,
              datatype = "FLT8S")

  ## zoom area for the figure: the 120 x 120 cell block, on a grid of block
  ## origins spaced 60 cells apart, with the largest mean absolute
  ## configurational component of lag-1 autocorrelation
  cf <- as.matrix(dk[["sr1_config"]], wide = TRUE)
  starts <- seq(8L, nrow(cf) - 7L - 119L, by = 60L)
  best <- c(NA, NA, -Inf)
  for (r0 in starts) for (c0 in starts) {
    val <- mean(abs(cf[r0:(r0 + 119L), c0:(c0 + 119L)]), na.rm = TRUE)
    if (val > best[3]) best <- c(r0, c0, val)
  }
  zoom <- data.frame(top_row = best[1], left_col = best[2], mean_abs_sr1_config = best[3],
                     xmin = xFromCol(z13, best[2]) - 15,
                     xmax = xFromCol(z13, best[2] + 119L) + 15,
                     ymin = yFromRow(z13, best[1] + 119L) - 15,
                     ymax = yFromRow(z13, best[1]) + 15)
  write.csv(zoom, here("output", "tables", "zoom_area.csv"), row.names = FALSE)
  message("15-cell windows done in ", round(proc.time()[3] - t0), " s")
}

## ---- series of window sizes -----------------------------------------------
if ("scale" %in% parts) {
  t0 <- proc.time()[3]
  sp <- list()
  for (w in c(5L, 9L, 15L, 31L, 61L)) {
    for (g in groups) {
      sp[[length(sp) + 1L]] <- ss_scale_profile(z13, z23, metrics = g, w_seq = w,
                                                range = rg, progress = FALSE)
      gc()
    }
    message("window ", w, " done after ", round(proc.time()[3] - t0), " s")
  }
  sp <- do.call(rbind, sp)
  sp <- sp[order(sp$w, match(sp$metric, metrics)), ]
  write.csv(sp, here("output", "tables", "decomposition_ndvi_scale_profile.csv"),
            row.names = FALSE)
}
