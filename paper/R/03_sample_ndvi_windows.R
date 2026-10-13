library(terra)
library(here)

z <- rast(here("data", "derived", "ndvi2013_aw30.tif"))
w <- 15L
h <- (w - 1L) %/% 2L
set.seed(20131L)
rows <- (h + 1L):(nrow(z) - h)
cols <- (h + 1L):(ncol(z) - h)
cand <- as.vector(outer((rows - 1L) * ncol(z), cols, "+"))
cells <- sort(sample(cand, 100L))
rc <- rowColFromCell(z, cells)
xy <- xyFromCell(z, cells)
win <- data.frame(window = seq_along(cells), cell = cells, row = rc[, 1], col = rc[, 2],
                  x = xy[, 1], y = xy[, 2], chunk = rep(1:10, each = 10L))
write.csv(win, here("output", "tables", "ndvi_windows_w15.csv"), row.names = FALSE)
