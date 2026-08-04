#' Catalogue of continuous landscape metrics
#'
#' Lists the metrics computed by the \pkg{surfscape} window kernel together with
#' the property the decomposition framework depends on: whether a metric
#' responds to the *arrangement* of values inside the analysis neighbourhood, or
#' only to the *multiset* of values it contains.
#'
#' @details
#' A metric belongs to the **composition** family if permuting the values within
#' a neighbourhood leaves it unchanged. Composition metrics can only ever
#' measure what values are present, never how they are organised; no spatial
#' rearrangement, however extreme, will move them. This is a property of the
#' definition rather than of the data, and it can be verified for any metric,
#' including metrics from other packages, with [ss_sensitivity()].
#'
#' Metrics belonging to the **configuration** family respond to arrangement and
#' are therefore the only ones capable of detecting a change in spatial
#' structure at a fixed value distribution.
#'
#' The amplitude family of surface metrics (Sa, Sq, Ssk, Sku, S10z) is aspatial,
#' as noted in passing by McGarigal et al. (2009), Frazier and Wang (2016) and
#' Kedron et al. (2018); so are the information-theoretic measures (Shannon
#' entropy, Rao's quadratic entropy), because they are computed from the value
#' distribution alone. The consequences of this for change detection are what
#' [ss_decompose()] addresses.
#'
#' @return A `data.frame` with columns `metric`, `id`, `family`, `arrangement` and
#'   `definition`.
#' @references
#' McGarigal K, Tagil S, Cushman SA (2009) Surface metrics: an alternative to
#' patch metrics for the quantification of landscape structure.
#' *Landscape Ecology* 24:433-450.
#'
#' Haralick RM, Shanmugam K, Dinstein I (1973) Textural features for image
#' classification. *IEEE Transactions on Systems, Man, and Cybernetics*
#' 3:610-621.
#'
#' Allain C, Cloitre M (1991) Characterizing the lacunarity of random and
#' deterministic fractal sets. *Physical Review A* 44:3552-3558.
#' @examples
#' ss_metric_list()
#' subset(ss_metric_list(), arrangement)     # configuration-sensitive only
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
      "Shannon entropy of binned values",
      "Rao's quadratic entropy of binned values",
      "Moran's I within the neighbourhood (rook adjacency)",
      "root mean square surface slope",
      "developed interfacial area ratio",
      "grey-scale gliding-box lacunarity",
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
#' Computes continuous landscape metrics in a moving window. Every metric is
#' derived from a single gather of each neighbourhood, and the neighbourhood is
#' held as a two-dimensional patch so that configuration-sensitive metrics
#' (within-window autocorrelation, slope, lacunarity, grey-level co-occurrence
#' textures) are available alongside the aspatial amplitude and
#' information-theoretic metrics.
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
#' bearing-area, spectral and correlation-length families, and \pkg{GLCMTextures}
#' implements grey-level co-occurrence textures with more options; both should
#' be preferred when those specific metrics are the goal, and both can be passed
#' to [ss_sensitivity()] as user functions. The kernel here exists because
#' [ss_decompose()] must evaluate the same metric on counterfactual
#' neighbourhoods that never exist as a raster and therefore cannot be handed to
#' an external package.
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
