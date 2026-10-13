## Software versions, checksums of the raw inputs and the lineage of the
## outputs used in the manuscript.

suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(surfscape)
  library(geodiv)
  library(GLCMTextures)
  library(glcm)
  library(rasterdiv)
})

writeLines(capture.output(sessionInfo()), here("logs", "sessionInfo.txt"))

raw <- list.files(here("data", "raw"), full.names = TRUE)
sha256 <- function(f) {
  out <- system2("sha256sum", shQuote(f), stdout = TRUE)
  sub(" .*$", "", out)
}
inputs <- data.frame(file = basename(raw), bytes = file.size(raw),
                     sha256 = vapply(raw, sha256, ""), role = "input")
write.csv(inputs, here("logs", "input_checksums_end.csv"), row.names = FALSE)

lineage <- read.csv(text = "output,script,inputs
output/tables/input_raster_summary.csv,R/00_inspect_inputs.R,data/raw/NDVI-2013.tif; data/raw/NDVI-2023.tif
data/derived/ndvi2013_aw30.tif,R/01_prepare_data.R,data/raw/NDVI-2013.tif; data/raw/NDVI-2023.tif
data/derived/ndvi2023_aw30.tif,R/01_prepare_data.R,data/raw/NDVI-2013.tif; data/raw/NDVI-2023.tif
output/tables/analysis_window.csv,R/01_prepare_data.R,data/raw/NDVI-2013.tif; data/raw/NDVI-2023.tif
output/tables/analysis_window_summary.csv,R/01_prepare_data.R,data/derived/ndvi*_aw*.tif
output/tables/simulated_landscapes.csv,R/02_simulate_landscapes.R,seeds in script
output/tables/ndvi_windows_w15.csv,R/03_sample_ndvi_windows.R,data/derived/ndvi2013_aw30.tif
output/tables/metric_catalogue.csv,R/11_metric_catalogue.R,package metadata
output/tables/benchmark_main_summary.csv,R/08_combine_benchmarks.R,output/rds/units/ndvi_chunk*.rds; output/rds/units/sim_*.rds (R/04_benchmark_units.R)
output/tables/benchmark_main_responses.csv,R/08_combine_benchmarks.R,output/rds/units/ndvi_chunk*.rds; output/rds/units/sim_*.rds (R/04_benchmark_units.R)
output/tables/benchmark_extent_grain_summary.csv,R/10_combine_extent_grain.R,output/rds/units/ext_*.rds; output/rds/units/grain*.rds (R/05_extent_grain_units.R)
output/tables/benchmark_extent_grain_responses.csv,R/10_combine_extent_grain.R,output/rds/units/ext_*.rds; output/rds/units/grain*.rds (R/05_extent_grain_units.R)
output/tables/benchmark_landscape_level_values.csv,R/10_combine_extent_grain.R,output/rds/units/ext_*.rds (R/05_extent_grain_units.R)
output/tables/benchmark_own_target_windows.csv,R/10_combine_extent_grain.R,output/rds/units/ext_*.rds; output/rds/units/grain*.rds (R/05_extent_grain_units.R)
output/tables/benchmark_geodiv_centred_summary.csv,R/09_geodiv_centred.R,data/derived/ndvi20*_aw30.tif; output/tables/ndvi_windows_w15.csv
output/tables/decomposition_virtual_windows.csv,R/06_decomposition_virtual.R,data/derived/sim/*.tif
output/tables/decomposition_virtual_landscape.csv,R/06_decomposition_virtual.R,data/derived/sim/*.tif
output/tables/decomposition_ndvi_landscape.csv,R/07_decomposition_ndvi.R,data/derived/ndvi20*_aw30.tif
output/tables/decomposition_ndvi_w15_cells.csv,R/07_decomposition_ndvi.R,data/derived/ndvi20*_aw30.tif
output/tables/decomposition_ndvi_scale_profile.csv,R/07_decomposition_ndvi.R,data/derived/ndvi20*_aw30.tif
output/rasters/decomposition_ndvi_w15.tif,R/07_decomposition_ndvi.R,data/derived/ndvi20*_aw30.tif
output/tables/zoom_area.csv,R/07_decomposition_ndvi.R,output/rasters/decomposition_ndvi_w15.tif
output/tables/redundancy_clusters.csv,R/13_redundancy.R,output/rds/benchmark_main_observed.rds; output/tables/benchmark_main_summary.csv
output/tables/subset_selected.csv,R/13_redundancy.R,output/rds/benchmark_main_observed.rds; output/tables/benchmark_main_summary.csv
output/tables/evenness_configuration_room.csv,R/15_evenness.R,output/rds/units/ndvi_chunk*.rds
output/tables/rasterdiv_class_check.csv,R/16_rasterdiv_check.R,data/derived/ndvi2013_aw30.tif; output/tables/ndvi_windows_w15.csv
output/tables/target_ties_check.csv,R/19_target_ties.R,data/derived/ndvi20*_aw*.tif
output/tables/numbers_for_text.csv,R/17_numbers.R,output/tables/*.csv
output/tables/grain_ratio_w15.csv,R/17_numbers.R,output/tables/benchmark_extent_grain_summary.csv
output/tables/extent_response_parts.csv,R/17_numbers.R,output/tables/benchmark_extent_grain_responses.csv
")
lineage$rows <- vapply(lineage$output, function(f) {
  p <- here(f)
  if (!file.exists(p)) return(NA_integer_)
  if (grepl("\\.csv$", f)) nrow(read.csv(p)) else NA_integer_
}, integer(1))
lineage$sha256 <- vapply(lineage$output, function(f) {
  p <- here(f)
  if (file.exists(p)) sha256(p) else NA_character_
}, "")
lineage$modified <- vapply(lineage$output, function(f) {
  p <- here(f)
  if (file.exists(p)) format(file.mtime(p), "%Y-%m-%d %H:%M") else NA_character_
}, "")
write.csv(lineage, here("logs", "data_lineage.csv"), row.names = FALSE)
print(inputs)
