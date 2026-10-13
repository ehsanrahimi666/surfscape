## Number of distinct 8-bit values in the NDVI benchmark windows, relevant to
## the class limit of rasterdiv::Shannon() (it tabulates values with
## summary(), which lists at most 100 classes).

suppressPackageStartupMessages({
  library(terra)
  library(here)
})
z <- rast(here("data", "derived", "ndvi2013_aw30.tif"))
ndw <- read.csv(here("output", "tables", "ndvi_windows_w15.csv"))
v <- values(z)[, 1]
rng <- range(v)
m <- as.matrix(z, wide = TRUE)
n_distinct <- vapply(seq_len(nrow(ndw)), function(i) {
  w <- m[(ndw$row[i] - 7):(ndw$row[i] + 7), (ndw$col[i] - 7):(ndw$col[i] + 7)]
  length(unique(as.integer(round(255 * (w - rng[1]) / (rng[2] - rng[1])))))
}, integer(1))
res <- data.frame(n_windows = length(n_distinct), max_distinct = max(n_distinct),
                  median_distinct = stats::median(n_distinct),
                  n_over_100 = sum(n_distinct > 100))
write.csv(res, here("output", "tables", "rasterdiv_class_check.csv"), row.names = FALSE)
print(res)
