#' surfscape: composition and configuration in continuous landscapes
#'
#' Tools to test whether a continuous landscape metric responds to the values
#' present in a landscape window (composition), to their spatial arrangement
#' (configuration), or to both, and to separate the compositional and
#' configurational parts of an observed change between two dates.
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{ss_benchmark}}}{responses of any metric to controlled
#'     changes of composition and configuration, and its classification}
#'   \item{\code{\link{ss_perturb}}}{the controlled changes used by the
#'     benchmark}
#'   \item{\code{\link{ss_simulate}}}{continuous neutral landscapes with known
#'     spatial structure}
#'   \item{\code{\link{ss_decompose}}}{split the change between two dates
#'     into global compositional, local compositional and configurational
#'     components}
#'   \item{\code{\link{ss_scale_profile}}}{the decomposition across window
#'     sizes}
#'   \item{\code{\link{ss_extract}}}{decomposed change in buffers around
#'     sites}
#'   \item{\code{\link{ss_metrics}}}{moving-window metric maps}
#'   \item{\code{\link{ss_rank_transfer}}}{counterfactual surfaces}
#' }
#'
#' @useDynLib surfscape, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
