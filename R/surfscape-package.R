#' surfscape: decomposing structural change in continuous landscapes
#'
#' Landscape ecology distinguishes composition from configuration, but the
#' continuous (gradient) representation of landscapes lost that distinction:
#' most metrics applied to continuous surfaces are functions of the value
#' distribution alone and are mathematically incapable of responding to how
#' those values are arranged in space. \pkg{surfscape} makes the distinction
#' explicit again.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{ss_decompose}}}{split change between two dates into
#'     distributional, redistributional and configurational components}
#'   \item{\code{\link{ss_sensitivity}}}{audit whether a metric can see spatial
#'     arrangement, and how much it does}
#'   \item{\code{\link{ss_scale_profile}}}{how the components scale}
#'   \item{\code{\link{ss_extract}}}{decomposed change in buffers around sites}
#'   \item{\code{\link{ss_metrics}}}{the underlying moving-window engine}
#'   \item{\code{\link{ss_rank_transfer}}}{counterfactual surface construction}
#' }
#'
#' @useDynLib surfscape, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
