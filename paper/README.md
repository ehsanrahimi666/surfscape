# Analyses for "Can continuous landscape metrics detect spatial configuration?"

R scripts that reproduce the analyses, tables and figures of the manuscript
from the NDVI surfaces of western Isfahan (Iran) for 2013 and 2023.

## Requirements

- R (>= 4.1; the analyses were run with R 4.6.1)
- R packages: surfscape 0.3.0, terra, geodiv, GLCMTextures, glcm, rasterdiv,
  here. The versions used are listed in `logs/sessionInfo.txt` after a run.
- The NDVI surfaces `NDVI-2013.tif` and `NDVI-2023.tif` (2077 x 2755 cells,
  30 m, WGS 84 / UTM zone 41N) in `data/raw/`. They are not distributed with
  the scripts; see the Data availability statement of the article.

## Running

From this folder:

```
Rscript run_all.R
```

runs every step in order. The complete run takes about 20 hours on one core.
The two benchmark steps can be split over several processes, for example

```
Rscript R/04_benchmark_units.R 1 2 &
Rscript R/04_benchmark_units.R 2 2 &
```

Each unit of work is saved in `output/rds/units/` and skipped if it exists,
so an interrupted run can be resumed.

## Steps

| script | output |
|---|---|
| `00_inspect_inputs.R` | structure of the raw rasters |
| `01_prepare_data.R` | 700 x 700 cell analysis window and its aggregation to 60 m and 120 m |
| `02_simulate_landscapes.R` | 30 Gaussian random fields (ranges 3, 10 and 30 cells) and their paired fields |
| `03_sample_ndvi_windows.R` | locations of the 100 NDVI benchmark windows |
| `11_metric_catalogue.R` | catalogue of the 75 metric implementations |
| `04_benchmark_units.R` | benchmark of all implementations, 15-cell windows |
| `05_extent_grain_units.R` | benchmark at window sizes of 5-61 cells, at the landscape level and at 60 m and 120 m |
| `06_decomposition_virtual.R` | controlled test of the decomposition |
| `07_decomposition_ndvi.R` | decomposition of the 2013-2023 NDVI change |
| `08_combine_benchmarks.R`, `10_combine_extent_grain.R` | pooled benchmark tables |
| `09_geodiv_centred.R` | geodiv autocorrelation and spectral metrics on centred windows |
| `13_redundancy.R` | correlations among metrics and the key subset |
| `15_evenness.R` | room for configuration and evenness of values |
| `16_rasterdiv_check.R` | number of distinct values per window (rasterdiv class limit) |
| `19_target_ties.R` | tied NDVI values and the observed-distribution test |
| `14_tables.R`, `12_figures.R` | tables and figures |
| `17_numbers.R` | every number quoted in the text, with its source file |
| `18_session_info.R` | software versions, input checksums and output lineage |

`external_metrics.R` contains the functions that evaluate the metrics of
geodiv, terra, GLCMTextures, glcm and rasterdiv on one landscape window.
Errors and warnings raised by these functions are recorded, not suppressed.

All random steps use fixed seeds, and the sampled window locations are
written to `output/tables/`.
