library(terra)
library(here)

r13 <- rast(here("data", "raw", "NDVI-2013.tif"))
r23 <- rast(here("data", "raw", "NDVI-2023.tif"))
stopifnot(compareGeom(r13, r23))

k <- 700L
valid <- as.matrix(!is.na(r13) & !is.na(r23), wide = TRUE)
storage.mode(valid) <- "integer"
S <- rbind(0, cbind(0, t(apply(apply(valid, 2, cumsum), 1, cumsum))))
nr <- nrow(valid)
nc <- ncol(valid)
tot <- S[(k + 1):(nr + 1), (k + 1):(nc + 1)] - S[1:(nr - k + 1), (k + 1):(nc + 1)] -
  S[(k + 1):(nr + 1), 1:(nc - k + 1)] + S[1:(nr - k + 1), 1:(nc - k + 1)]
feasible <- which(tot == k * k, arr.ind = TRUE)
centre <- colMeans(feasible)
d2 <- (feasible[, 1] - centre[1])^2 + (feasible[, 2] - centre[2])^2
cand <- feasible[d2 == min(d2), , drop = FALSE]
cand <- cand[order(cand[, 1], cand[, 2]), , drop = FALSE]
top <- cand[1, 1]
left <- cand[1, 2]

aw13 <- r13[top:(top + k - 1L), left:(left + k - 1L), drop = FALSE]
aw23 <- r23[top:(top + k - 1L), left:(left + k - 1L), drop = FALSE]
names(aw13) <- "ndvi_2013"
names(aw23) <- "ndvi_2023"
stopifnot(!anyNA(values(aw13)), !anyNA(values(aw23)))
writeRaster(aw13, here("data", "derived", "ndvi2013_aw30.tif"), overwrite = TRUE,
            datatype = "FLT4S")
writeRaster(aw23, here("data", "derived", "ndvi2023_aw30.tif"), overwrite = TRUE,
            datatype = "FLT4S")

for (f in c(2L, 4L)) {
  for (yr in c("2013", "2023")) {
    a <- if (yr == "2013") aw13 else aw23
    g <- aggregate(a, fact = f, fun = "mean")
    writeRaster(g, here("data", "derived", sprintf("ndvi%s_aw%dm.tif", yr, 30L * f)),
                overwrite = TRUE, datatype = "FLT8S")
  }
}

ctr <- vect(cbind((xmin(aw13) + xmax(aw13)) / 2, (ymin(aw13) + ymax(aw13)) / 2),
            crs = crs(aw13))
ll <- crds(project(ctr, "EPSG:4326"))
lat <- ll[1, 2] * pi / 180
dlon <- (ll[1, 1] - 63) * pi / 180
k_tm <- 0.9996 / sqrt(1 - (cos(lat) * sin(dlon))^2)

win <- data.frame(
  rule = "fully valid 700 x 700 window closest to the centroid of all fully valid positions",
  n_feasible_positions = nrow(feasible), top_row = top, left_col = left,
  xmin = xmin(aw13), xmax = xmax(aw13), ymin = ymin(aw13), ymax = ymax(aw13),
  crs = crs(aw13, describe = TRUE)$name, epsg = crs(aw13, describe = TRUE)$code,
  centre_lon = ll[1, 1], centre_lat = ll[1, 2],
  utm_zone_of_centre = floor((ll[1, 1] + 180) / 6) + 1,
  scale_factor_zone41 = k_tm)
write.csv(win, here("output", "tables", "analysis_window.csv"), row.names = FALSE)

psd <- function(v) sqrt(mean((v - mean(v))^2))
files <- list.files(here("data", "derived"), pattern = "^ndvi.*\\.tif$", full.names = TRUE)
summ <- do.call(rbind, lapply(files, function(f) {
  r <- rast(f)
  v <- values(r)[, 1]
  data.frame(file = basename(f), nrow = nrow(r), ncol = ncol(r), res = res(r)[1],
             n = sum(!is.na(v)), n_na = sum(is.na(v)), min = min(v), max = max(v),
             mean = mean(v), sd = psd(v),
             skewness = mean((v - mean(v))^3) / psd(v)^3)
}))
write.csv(summ, here("output", "tables", "analysis_window_summary.csv"), row.names = FALSE)
print(win)
print(summ)
