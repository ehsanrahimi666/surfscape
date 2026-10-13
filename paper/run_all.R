## Runs the complete analysis from the raw NDVI files. Each step is a separate
## R process, so that the steps can also be run one at a time. The benchmark
## steps (04, 05) accept a worker number and the number of workers, so that
## they can be split over several processes; here they run as one worker.
## The complete run takes about 20 hours of computing time on one core.

steps <- list(
  c("R/00_inspect_inputs.R"),
  c("R/01_prepare_data.R"),
  c("R/02_simulate_landscapes.R"),
  c("R/03_sample_ndvi_windows.R"),
  c("R/11_metric_catalogue.R"),
  c("R/04_benchmark_units.R", "1", "1"),
  c("R/05_extent_grain_units.R", "1", "1"),
  c("R/06_decomposition_virtual.R"),
  c("R/07_decomposition_ndvi.R"),
  c("R/08_combine_benchmarks.R"),
  c("R/09_geodiv_centred.R"),
  c("R/10_combine_extent_grain.R"),
  c("R/13_redundancy.R"),
  c("R/15_evenness.R"),
  c("R/16_rasterdiv_check.R"),
  c("R/19_target_ties.R"),
  c("R/14_tables.R"),
  c("R/12_figures.R"),
  c("R/17_numbers.R"),
  c("R/18_session_info.R")
)
for (s in steps) {
  message(format(Sys.time()), " running ", paste(s, collapse = " "))
  status <- system2("Rscript", s)
  if (status != 0) stop("step failed: ", paste(s, collapse = " "))
}
