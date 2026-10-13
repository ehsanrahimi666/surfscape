## Catalogue of the 75 metric implementations (Appendix Table A1).

suppressPackageStartupMessages({
  library(here)
  library(surfscape)
})

bi <- ss_metric_list()
builtin <- data.frame(
  implementation = paste0("surfscape_", bi$metric), package = "surfscape 0.3.0",
  call = paste0("ss_metrics(..., metrics = \"", bi$metric, "\")"),
  definition = bi$definition,
  settings = ifelse(bi$metric %in% c("sent", "srao"), "32 bins",
             ifelse(bi$metric == "slac", "3 × 3 gliding box",
             ifelse(grepl("^g", bi$metric), "16 grey levels, 4 directions pooled",
             ifelse(bi$metric %in% c("sdq", "sdr"), "central differences, cell size 30",
             ifelse(bi$metric == "sr1", "rook neighbours", ""))))),
  binned = bi$binned)

gd <- rbind(
  c("sa", "mean absolute deviation of values from their mean"),
  c("sq", "standard deviation of values (sample)"),
  c("ssk", "skewness of values"),
  c("sku", "kurtosis of values"),
  c("s10z", "mean height of the five highest local maxima plus mean depth of the five lowest local minima"),
  c("smean", "mean of the positive values (mean peak height)"),
  c("sv", "absolute value of the lowest value (maximum valley depth)"),
  c("sph", "absolute value of the highest value (maximum peak height)"),
  c("sbi", "surface bearing index: Sq divided by the height at 5% of the bearing area"),
  c("sci", "core fluid retention index from the bearing area curve"),
  c("svi", "valley fluid retention index from the bearing area curve"),
  c("sk", "core roughness depth from the bearing area curve"),
  c("spk", "reduced peak height from the bearing area curve"),
  c("svk", "reduced valley depth from the bearing area curve"),
  c("sdq", "root mean square slope, two-point method"),
  c("sdq6", "area root mean square slope, seven-point method"),
  c("sdr", "surface area ratio"),
  c("sds", "summit density, number of local peaks per unit area"),
  c("ssc", "mean summit curvature of local maxima"),
  c("sfd", "fractal dimension, triangular prism method"),
  c("scl20", "smallest distance at which the areal autocorrelation function falls to 0.2"),
  c("str20", "texture aspect ratio at autocorrelation 0.2"),
  c("srw", "dominant radial wavelength of the Fourier spectrum"),
  c("srwi", "radial wavelength index"),
  c("shw", "mean half wavelength of the radial spectrum"))
geodiv <- data.frame(
  implementation = paste0("geodiv_", gd[, 1]), package = "geodiv 1.1.1",
  call = paste0(sub("^(scl20)$", "scl", sub("^(str20)$", "stxr",
               sub("^(srwi|shw)$", "srw", gd[, 1]))), "(window)"),
  definition = gd[, 2], settings = "window passed as SpatRaster, cell size 30",
  binned = FALSE)

tx <- c(contrast = "contrast", dissimilarity = "dissimilarity", homogeneity = "homogeneity",
        ASM = "angular second moment", entropy = "entropy", mean = "GLCM mean",
        variance = "GLCM variance", correlation = "correlation")
glcmtex <- do.call(rbind, lapply(c("range", "prob"), function(q) data.frame(
  implementation = paste0("glcmtex_", q, "_", names(tx)), package = "GLCMTextures 0.6.3",
  call = "glcm_textures(quantized window, w = NULL, n_levels = 16, quant_method = \"none\")",
  definition = paste("GLCM", tx),
  settings = paste0("16 grey levels, ", if (q == "range") "equal-range" else "equal-probability",
                    " quantization, shifts (1,0), (1,1), (0,1), (-1,1)"),
  binned = TRUE)))
gp <- c(contrast = "contrast", dissimilarity = "dissimilarity", homogeneity = "homogeneity",
        second_moment = "second moment", entropy = "entropy", mean = "mean",
        variance = "variance", correlation = "correlation")
glcmpkg <- data.frame(
  implementation = paste0("glcmpkg_", names(gp)), package = "glcm 1.6.6",
  call = "glcm(window, n_grey = 16, window = c(13, 13), shift = 4 directions)",
  definition = paste("GLCM", gp), settings = "16 grey levels, equal-range quantization, statistic averaged over shifts (0,1), (1,1), (1,0), (1,-1)",
  binned = TRUE)
others <- data.frame(
  implementation = c("terra_moran", "terra_geary", "rasterdiv_shannon", "rasterdiv_rao"),
  package = c("terra 1.9-50", "terra 1.9-50", "rasterdiv 0.3.8", "rasterdiv 0.3.8"),
  call = c("autocor(window, method = \"moran\", global = TRUE)",
           "autocor(window, method = \"geary\", global = TRUE)",
           "Shannon(window, window = 15, na.tolerance = 0)",
           "paRao(window, window = 15, alpha = 1, na.tolerance = 0)"),
  definition = c("global Moran's I", "global Geary's C", "Shannon entropy of integer classes",
                 "parametric Rao's quadratic entropy, mean absolute difference between values"),
  settings = c("queen neighbours", "queen neighbours", "values rescaled to 8-bit integers",
               "values rescaled to 8-bit integers"),
  binned = c(FALSE, FALSE, TRUE, TRUE))

cat_tab <- rbind(builtin, geodiv, others[1:2, ], glcmtex, glcmpkg, others[3:4, ])
stopifnot(nrow(cat_tab) == 75L, !anyDuplicated(cat_tab$implementation))
write.csv(cat_tab, here("output", "tables", "metric_catalogue.csv"), row.names = FALSE)
