## Tied values and the observed-distribution (target) test. Rank transfer
## breaks ties by cell order, so cells with equal 2013 values can receive
## different 2023 values and cells with different 2013 values can receive equal
## ones. This counts such cases and the cells whose equal-probability grey
## level (16 levels, as in the benchmark) differs between each 2013 surface and
## its target-test counterpart.

suppressPackageStartupMessages({
  library(terra)
  library(here)
  library(surfscape)
})
source(here("R", "external_metrics.R"))

grains <- c(aw30 = 30L, aw60m = 60L, aw120m = 120L)
res <- do.call(rbind, lapply(names(grains), function(sfx) {
  z <- rast(here("data", "derived", sprintf("ndvi2013_%s.tif", sfx)))
  z23 <- rast(here("data", "derived", sprintf("ndvi2023_%s.tif", sfx)))
  v <- values(z)[, 1]
  v23 <- values(z23)[, 1]
  vt <- values(ss_perturb(z, "target", target = z23))[, 1]
  o <- order(v)
  g0 <- .quantize(matrix(v), .prob_breaks(v, glcm_levels))
  g1 <- .quantize(matrix(vt), .prob_breaks(vt, glcm_levels))
  data.frame(grain = grains[[sfx]], cells = length(v),
             distinct_2013 = length(unique(v)), distinct_2023 = length(unique(v23)),
             cells_tied_2013 = sum(duplicated(v) | duplicated(v, fromLast = TRUE)),
             target_has_2023_values = identical(sort(vt), sort(v23)),
             target_nondecreasing = !is.unsorted(vt[o]),
             ties_split = sum(diff(v[o]) == 0 & diff(vt[o]) != 0),
             ties_created = sum(diff(v[o]) != 0 & diff(vt[o]) == 0),
             cells_grey_level_changed = sum(g0 != g1))
}))
write.csv(res, here("output", "tables", "target_ties_check.csv"), row.names = FALSE)
print(res)
