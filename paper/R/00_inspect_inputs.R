library(terra)

files <- c(ndvi_2013 = "data/raw/NDVI-2013.tif", ndvi_2023 = "data/raw/NDVI-2023.tif")
info <- lapply(names(files), function(nm) {
  r <- rast(files[[nm]])
  v <- values(r, mat = FALSE)
  data.frame(layer = nm, nrow = nrow(r), ncol = ncol(r), nlyr = nlyr(r),
             res_x = res(r)[1], res_y = res(r)[2],
             crs_epsg = crs(r, describe = TRUE)$code,
             crs_name = crs(r, describe = TRUE)$name,
             datatype = datatype(r), n_cells = ncell(r),
             n_na = sum(is.na(v)), n_valid = sum(!is.na(v)),
             min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE),
             mean = mean(v, na.rm = TRUE), sd = sd(v, na.rm = TRUE))
})
info <- do.call(rbind, info)
write.csv(info, "output/tables/input_raster_summary.csv", row.names = FALSE)
print(info)

r1 <- rast(files[[1]]); r2 <- rast(files[[2]])
cat("identical geometry:", compareGeom(r1, r2, stopOnError = FALSE), "\n")
both <- !is.na(r1) & !is.na(r2)
cat("cells valid in both:", global(both, "sum")[1, 1], "\n")
cat("metadata tags 2013:\n"); print(metags(r1))
cat("time 2013:", as.character(time(r1)), "\n")
ctr <- vect(cbind(mean(c(xmin(r1), xmax(r1))), mean(c(ymin(r1), ymax(r1)))), crs = crs(r1))
cat("scene centre lon/lat:", round(crds(project(ctr, "EPSG:4326")), 4), "\n")
