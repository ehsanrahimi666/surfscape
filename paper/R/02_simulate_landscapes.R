library(terra)
library(here)
library(surfscape)

ranges <- c(3, 10, 30)
nrep <- 10L
size <- 201L

grid <- expand.grid(rep = seq_len(nrep), range = ranges)
grid$id <- sprintf("r%02d_rep%02d", grid$range, grid$rep)
grid$seed_landscape <- 100000L + 100L * grid$range + grid$rep
grid$seed_arrangement <- 200000L + 100L * grid$range + grid$rep
grid$file_landscape <- file.path("data", "derived", "sim", paste0("grf_", grid$id, ".tif"))
grid$file_arrangement <- file.path("data", "derived", "sim", paste0("arr_", grid$id, ".tif"))

lag1 <- function(r) {
  m <- as.matrix(r, wide = TRUE)
  cor(c(as.vector(m[, -1]), as.vector(m[-1, ])),
      c(as.vector(m[, -ncol(m)]), as.vector(m[-nrow(m), ])))
}
grid$lag1_landscape <- NA_real_
grid$lag1_expected <- exp(-3 / grid$range)
for (i in seq_len(nrow(grid))) {
  z <- ss_simulate(size, size, range = grid$range[i], res = 30, seed = grid$seed_landscape[i])
  a <- ss_simulate(size, size, range = grid$range[i], res = 30, seed = grid$seed_arrangement[i])
  writeRaster(z, here(grid$file_landscape[i]), overwrite = TRUE, datatype = "FLT8S")
  writeRaster(a, here(grid$file_arrangement[i]), overwrite = TRUE, datatype = "FLT8S")
  grid$lag1_landscape[i] <- lag1(z)
}
write.csv(grid, here("output", "tables", "simulated_landscapes.csv"), row.names = FALSE)
print(aggregate(cbind(lag1_landscape, lag1_expected) ~ range, grid, mean))
