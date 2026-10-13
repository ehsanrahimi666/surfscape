## Wrappers that evaluate metrics implemented in other R packages on one
## landscape window, for use as `fun` in surfscape::ss_benchmark().
##
## Failures of individual metrics and warnings are not discarded: they are
## counted and their messages written to the log returned by ext_log().

.ext_env <- new.env()
.ext_env$log <- list()
.ext_env$qcache <- list(key = NULL, br = NULL)

ext_log_reset <- function() .ext_env$log <- list()

ext_log <- function() {
  if (!length(.ext_env$log)) {
    return(data.frame(metric = character(), type = character(),
                      message = character(), n = integer()))
  }
  d <- do.call(rbind, .ext_env$log)
  a <- aggregate(n ~ metric + type + message, data = transform(d, n = 1L), FUN = sum)
  a[order(a$metric, a$type), ]
}

.log_event <- function(metric, type, message) {
  .ext_env$log[[length(.ext_env$log) + 1L]] <-
    data.frame(metric = metric, type = type, message = message)
}

## Evaluate `expr` for metric `nm`. Errors give NA and are logged; warnings
## whose message matches `known` are logged and muffled, others are logged and
## propagated.
.guard <- function(nm, expr, n = 1L, known = character()) {
  withCallingHandlers(
    tryCatch(expr, error = function(e) {
      .log_event(nm, "error", conditionMessage(e))
      rep(NA_real_, n)
    }),
    warning = function(w) {
      msg <- conditionMessage(w)
      .log_event(nm, "warning", msg)
      if (any(vapply(known, grepl, logical(1), x = msg, fixed = TRUE))) {
        invokeRestart("muffleWarning")
      }
    })
}

## ---- geodiv and terra: metrics that do not bin values ---------------------

geodiv_single <- c("sa", "sq", "ssk", "sku", "s10z", "smean", "sv", "sph",
                   "sbi", "sci", "svi", "sk", "spk", "svk",
                   "sdq", "sdq6", "sdr", "sds", "ssc")

unbinned_metrics <- function(r) {
  gd <- asNamespace("geodiv")
  out <- vapply(geodiv_single, function(fn)
    as.numeric(.guard(paste0("geodiv_", fn), get(fn, gd)(r)))[1], numeric(1))
  names(out) <- paste0("geodiv_", geodiv_single)
  sfd <- .guard("geodiv_sfd", geodiv::sfd(r, silent = TRUE))
  scl <- .guard("geodiv_scl", unlist(geodiv::scl(r)), n = 4L)
  stx <- .guard("geodiv_stxr", unlist(geodiv::stxr(r)), n = 2L)
  srw <- .guard("geodiv_srw", unlist(geodiv::srw(r)), n = 3L)
  moran <- .guard("terra_moran", terra::autocor(r, method = "moran", global = TRUE))
  geary <- .guard("terra_geary", terra::autocor(r, method = "geary", global = TRUE))
  c(out,
    geodiv_sfd = as.numeric(sfd)[1],
    geodiv_scl20 = as.numeric(scl)[1],
    geodiv_str20 = as.numeric(stx)[1],
    geodiv_srw = as.numeric(srw)[1],
    geodiv_srwi = as.numeric(srw)[2],
    geodiv_shw = as.numeric(srw)[3],
    terra_moran = as.numeric(moran)[1],
    terra_geary = as.numeric(geary)[1])
}

## terra only, for the extent and grain analyses
terra_metrics <- function(r) {
  c(terra_moran = as.numeric(.guard("terra_moran",
      terra::autocor(r, method = "moran", global = TRUE)))[1],
    terra_geary = as.numeric(.guard("terra_geary",
      terra::autocor(r, method = "geary", global = TRUE)))[1])
}

## ---- binned metrics: GLCMTextures, rasterdiv, glcm ------------------------

glcm_levels <- 16L

## grey levels 0..L-1 as in GLCMTextures::quantize_raster(): breaks with the
## outer limits set to -Inf and Inf, intervals closed on the left
.quantize <- function(m, br) {
  matrix(findInterval(m, br[2:(length(br) - 1L)]), nrow(m))
}

.prob_breaks <- function(surface, L) {
  key <- c(length(surface), sum(surface), surface[1], surface[length(surface)])
  if (!identical(key, .ext_env$qcache$key)) {
    .ext_env$qcache <- list(key = key, br = stats::quantile(
      surface, probs = seq(0, 1, length.out = L + 1L), type = 8, names = FALSE))
  }
  .ext_env$qcache$br
}

.glcmtex <- function(q, prefix) {
  g <- .guard(prefix, GLCMTextures::glcm_textures(q, w = NULL, n_levels = glcm_levels,
                                                    quant_method = "none"), n = 8L)
  stats::setNames(as.numeric(g), paste0(prefix, "_", sub("^glcm_", "", c(
    "glcm_contrast", "glcm_dissimilarity", "glcm_homogeneity", "glcm_ASM",
    "glcm_entropy", "glcm_mean", "glcm_variance", "glcm_correlation"))))
}

glcm_texture_metrics <- function(m, surface) {
  rng <- range(surface)
  q_range <- .quantize(m, seq(rng[1], rng[2], length.out = glcm_levels + 1L))
  q_prob <- .quantize(m, .prob_breaks(surface, glcm_levels))
  c(.glcmtex(q_range, "glcmtex_range"), .glcmtex(q_prob, "glcmtex_prob"))
}

rasterdiv_known <- c("NAs introduced by coercion", "Simplify=0")

## GLCMTextures (equal-range and equal-probability grey levels) and glcm
texture_metrics <- function(m, surface) {
  rng <- range(surface)
  ctr <- c((nrow(m) + 1L) %/% 2L, (ncol(m) + 1L) %/% 2L)
  ## glcm::glcm pairs a base window with shifted copies of it; a base window
  ## two cells smaller keeps every pair inside the landscape window
  st <- c("contrast", "dissimilarity", "homogeneity", "second_moment",
          "entropy", "mean", "variance", "correlation")
  gp <- .guard("glcmpkg", glcm::glcm(m, n_grey = glcm_levels, window = dim(m) - 2L,
                                      shift = list(c(0, 1), c(1, 1), c(1, 0), c(1, -1)),
                                      statistics = st, min_x = rng[1], max_x = rng[2],
                                      na_opt = "any")[ctr[1], ctr[2], ],
               n = length(st))
  c(glcm_texture_metrics(m, surface),
    stats::setNames(as.numeric(gp), paste0("glcmpkg_", st)))
}

## rasterdiv on values rescaled to 8-bit integers over the binning range.
## The window is passed as a matrix of its own size; with na.tolerance = 0 only
## the complete (central) moving window is computed.
rasterdiv_metrics <- function(m, surface) {
  rng <- range(surface)
  x8 <- matrix(as.integer(round(255 * (m - rng[1]) / (rng[2] - rng[1]))), nrow(m))
  ctr <- c((nrow(m) + 1L) %/% 2L, (ncol(m) + 1L) %/% 2L)
  sh <- .guard("rasterdiv_shannon", suppressMessages(
    rasterdiv::Shannon(x8, window = nrow(m), np = 1, na.tolerance = 0))[ctr[1], ctr[2]])
  rao <- .guard("rasterdiv_rao", suppressMessages(
    rasterdiv::paRao(x8, window = nrow(m), alpha = 1, np = 1, progBar = FALSE,
                     na.tolerance = 0))[[1]][[1]][ctr[1], ctr[2]],
    known = rasterdiv_known)
  c(rasterdiv_shannon = as.numeric(sh), rasterdiv_rao = as.numeric(rao))
}
