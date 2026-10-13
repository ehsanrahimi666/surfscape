#' Catalogue of continuous landscape metrics
#'
#' Lists the metrics computed by the \pkg{surfscape} window kernel with two
#' properties that follow from their definitions: whether a metric depends on
#' the positions of values inside the landscape window (`arrangement`), and
#' whether it bins values on a value range (`binned`), which makes its
#' response to changes in composition depend on how that range is set.
#'
#' @details
#' A metric with `arrangement = FALSE` is a function of the set of values in
#' the window only: rearranging the values inside the window leaves it
#' unchanged, so it describes composition and cannot respond to configuration.
#' This holds for the amplitude metrics of surface metrology (McGarigal et al.
#' 2009) and for the information-theoretic measures, as noted by Frazier
#' (2016) and Kedron et al. (2018). A metric with `arrangement = TRUE` can
#' respond to configuration, but that does not make it insensitive to
#' composition; [ss_benchmark()] measures both responses for any metric.
#'
#' @return A `data.frame` with columns `metric`, `id`, `family`,
#'   `arrangement`, `binned` and `definition`.
#' @references
#' McGarigal K, Tagil S, Cushman SA (2009) Surface metrics: an alternative to
#' patch metrics for the quantification of landscape structure.
#' *Landscape Ecology* 24:433-450.
#'
#' Frazier AE (2016) Surface metrics: scaling relationships and downscaling
#' behavior. *Landscape Ecology* 31:351-363.
#'
#' Kedron PJ, Frazier AE, Ovando-Montejo GA, Wang J (2018) Surface metrics for
#' landscape ecology: a comparison of landscape models across ecoregions and
#' scales. *Landscape Ecology* 33:1489-1504.
#'
#' Haralick RM, Shanmugam K, Dinstein I (1973) Textural features for image
#' classification. *IEEE Transactions on Systems, Man, and Cybernetics*
#' 3:610-621.
#'
#' Allain C, Cloitre M (1991) Characterizing the lacunarity of random and
#' deterministic fractal sets. *Physical Review A* 44:3552-3558.
#'
#' Plotnick RE, Gardner RH, Hargrove WW, Prestegaard K, Perlmutter M (1996)
#' Lacunarity analysis: a general technique for the analysis of spatial
#' patterns. *Physical Review E* 53:5461-5468.
#' @examples
#' ss_metric_list()
#' subset(ss_metric_list(), arrangement)
#' @export
ss_metric_list <- function() {
  data.frame(
    metric = c("smean", "sa", "sq", "ssk", "sku", "sp", "sv", "sz",
               "s10z", "scv", "sent", "srao",
               "sr1", "sdq", "sdr", "slac",
               "gcon", "gdis", "ghom", "gasm", "gent", "gcor"),
    id = c(1:12, 20:23, 30:35),
    family = c(rep("amplitude", 10L), rep("information", 2L),
               "autocorrelation", rep("gradient", 2L), "lacunarity",
               rep("texture", 6L)),
    arrangement = c(rep(FALSE, 12L), rep(TRUE, 10L)),
    binned = c(rep(FALSE, 10L), TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
               rep(TRUE, 6L)),
    definition = c(
      "window mean",
      "average roughness, mean(|z - zbar|)",
      "root mean square roughness (standard deviation)",
      "skewness of the value distribution",
      "kurtosis of the value distribution",
      "maximum peak height, max(z) - zbar",
      "maximum valley depth, zbar - min(z)",
      "range, max(z) - min(z)",
      "ten-point height, mean of five highest minus five lowest",
      "coefficient of variation",
      "Shannon entropy of values binned on the value range",
      "Rao's quadratic entropy of values binned on the value range",
      "lag-1 spatial autocorrelation, 1 - semivariance at lag 1 / variance (rook neighbours)",
      "root mean square surface slope",
      "developed interfacial area ratio",
      "grey-scale gliding-box lacunarity (box mass measured from the lower end of the value range)",
      "GLCM contrast",
      "GLCM dissimilarity",
      "GLCM homogeneity",
      "GLCM angular second moment (energy)",
      "GLCM entropy",
      "GLCM correlation"
    ),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.ss_metric_id <- function(metrics) {
  tab <- ss_metric_list()
  bad <- setdiff(metrics, tab$metric)
  if (length(bad)) {
    stop("Unknown metric(s): ", paste(bad, collapse = ", "),
         ". See ss_metric_list().", call. = FALSE)
  }
  as.integer(tab$id[match(metrics, tab$metric)])
}

#' Moving-window metrics for continuous landscape surfaces
#'
#' Computes continuous landscape metrics in a moving window. All metrics are
#' computed from the same square landscape window around each cell, so that
#' metrics that depend on the arrangement of values (lag-1 autocorrelation,
#' slope, lacunarity, grey-level co-occurrence textures) and metrics that
#' depend only on the values present (amplitude and information-theoretic
#' metrics) describe exactly the same area.
#'
#' @param x A single-layer `SpatRaster`, or anything [terra::rast()] accepts.
#' @param metrics Character vector of metric names; see [ss_metric_list()].
#' @param w Window size: a single odd integer, or `c(nrow, ncol)`, in cells.
#' @param range Numeric length-2 value range used to bin `sent`, `srao` and the
#'   grey levels of the GLCM metrics. Defaults to the range of `x`. **Supply the
#'   same `range` for both dates when comparing metrics through time**,
#'   otherwise the binning differs between dates and the comparison is not
#'   meaningful; [ss_decompose()] does this automatically.
#' @param nbins Number of bins for `sent` and `srao`.
#' @param nlevels Number of grey levels for the GLCM metrics.
#' @param lacbox Gliding-box side length for `slac`.
#' @param filename,overwrite,chunk_rows,wopt Output and blocking options.
#'
#' @return A `SpatRaster` with one layer per requested metric, in the order
#'   given by `metrics`.
#'
#' @section Relationship to other packages:
#' \pkg{geodiv} implements a wider set of surface metrics, including the
#' bearing-area, spectral and correlation-length families, \pkg{GLCMTextures}
#' implements grey-level co-occurrence textures with more options, and
#' [terra::autocor()] computes Moran's I and Geary's C. These should be
#' preferred when those metrics are the goal, and any of them can be tested
#' with [ss_benchmark()]. The metrics are reimplemented here only because
#' [ss_decompose()] evaluates them on artificial windows created during the
#' computation (for example, the values of one date in the arrangement of the
#' other), which are not available as raster layers that could be passed to
#' another package.
#'
#' @examples
#' f <- system.file("extdata", "ndvi_2013.tif", package = "surfscape")
#' if (nzchar(f)) {
#'   r <- terra::rast(f)
#'   m <- ss_metrics(r, c("sq", "sr1", "gcon"), w = 9)
#'   names(m)
#' }
#' @export
ss_metrics <- function(x,
                       metrics = c("sq", "ssk", "sr1", "gcon"),
                       w = 9,
                       range = NULL,
                       nbins = 32L,
                       nlevels = 16L,
                       lacbox = 3L,
                       filename = "",
                       overwrite = FALSE,
                       chunk_rows = NULL,
                       wopt = list()) {

  x <- .ss_layer(x)
  w <- .ss_window(w)
  metrics <- unique(metrics)
  ids <- .ss_metric_id(metrics)
  pars <- .ss_pars(x, nbins, range, nlevels, lacbox)
  halo <- (w[1] - 1L) %/% 2L

  fun <- function(vlist, nr, nc)
    ss_focal_metrics_cpp(vlist[[1]], nr, nc, w[1], w[2], ids, pars)

  .ss_block_apply(list(x), halo, fun, metrics, filename = filename,
                  overwrite = overwrite, chunk_rows = chunk_rows, wopt = wopt)
}
