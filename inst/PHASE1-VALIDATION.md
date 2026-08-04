# surfscape 0.1.0 — Phase 1 validation and benchmarks

Environment: R 4.3.3, terra 1.7-65, Rcpp 1.0.12, g++ 13.3 (`-O2`), single core.
`R CMD check`: 0 errors, 0 package warnings, 37/37 tests passing.

## 1. Numerical correctness of the window engine

Every window metric was compared against an independent R implementation that
extracts each window explicitly and evaluates the definition directly
(37 × 41 Gaussian raster, 5 × 5 window). Maximum absolute differences:

| Metric | Max abs. difference |
|---|---|
| `smean` | 3.3e-16 |
| `sa` | 4.4e-16 |
| `sq` | 4.4e-16 |
| `ssk` | 1.0e-15 |
| `sku` | 5.3e-15 |
| `sp`, `sv` | 4.4e-16 |
| `sz` | 0 |
| `scv` | 1.1e-10 (ratio metric; amplified when the window mean is near zero) |

`s10z` and `sent` were checked against separate reference implementations to the
same tolerance.

**Gradient metrics.** On a plane `z = 0.7x`, `sdq` returns 0.7 and `sdr`
returns `sqrt(1 + 0.7^2) - 1`, both to 1e-8.

**Block invariance.** Results with `chunk_rows = 1`, `chunk_rows = 3` and
terra's automatic blocking are bit-identical (`expect_identical` on the raw
values), including for metrics that depend on global quantities.

**Missing data.** With `na_rm = TRUE` the number of missing output cells equals
the number of missing input cells; with `na_rm = FALSE` it increases as
expected (200 → 1450 for 200 scattered `NA`s and a 5 × 5 window).

## 2. Correctness of the local association statistics

| Statistic | Reference | Max abs. difference |
|---|---|---|
| Gi\* statistic | closed form | 2.4e-12 |
| Gi\* z-score | Ord & Getis (1995) | 8.9e-16 |
| Gi statistic (focal excluded) | closed form | 0 |
| Gi z-score | Ord & Getis (1995) | 2.6e-14 |
| Local Moran's I | Anselin (1995) | 1.2e-14 |

## 3. Calibration of inference (the substantive result)

200 independent Gaussian 40 × 40 fields, queen neighbourhood, interior cells
only (n = 180,000 evaluations). Realised rejection rates against nominal
levels:

| Statistic / inference | sd(z) | α = 0.10 | α = 0.05 | α = 0.01 | α = 0.001 |
|---|---|---|---|---|---|
| Gi\*, analytical | 0.997 | 0.0989 | 0.0491 | 0.0096 | 0.00100 |
| Local Moran, analytical | 0.994 | 0.0937 | 0.0642 | 0.0305 | 0.01364 |
| Local Moran, Monte-Carlo (999) | — | 0.1006 | 0.0498 | 0.0096 | — |

Two conclusions:

1. **The variance formulae are right** — z-scores have unit standard deviation
   for both statistics, confirming the implemented randomisation moments.
2. **The normal approximation is right for Gi\* and wrong for local Moran.**
   Local Moran's I is a product of a value and its spatial lag, so its
   randomisation distribution has much heavier tails than a Gaussian. At
   α = 0.01 the analytical test is three times too liberal, and at α = 0.001 it
   is fourteen times too liberal. Because a per-cell test on a raster is
   followed by a multiple-testing correction, the effective threshold is pushed
   into exactly the region where the approximation fails worst. Published local
   Moran cluster maps of NDVI and similar surfaces that rely on
   normal-approximation p-values will therefore over-declare significance.

The Monte-Carlo alternative is exactly calibrated: Kolmogorov–Smirnov test of
the p-values against a uniform distribution gives D = 0.0107, p = 0.89.

**Why it is cheap.** Once the focal value is held fixed, `I_i = z_i L_i / m2`
is a fixed multiple of the spatial lag `L_i`, so the two-sided randomisation
p-value of `I_i` equals that of `L_i`. The null distribution of `L_i` depends
only on the number of contributing neighbours, so it is simulated once per
distinct neighbour count (a handful of values for any window) instead of once
per cell. This is O(K · nsim · k) rather than the O(n · nsim · k) of a per-cell
conditional permutation.

## 4. Performance

500 × 500 raster, four metrics (`sa`, `sq`, `ssk`, `sku`), single core. The
comparison is against one `terra::focal()` pass per metric with an R callback,
which is the structure used by per-metric moving-window implementations; it is
a proxy for `geodiv::focal_metrics()`, which could not be installed in this
environment.

| Window | `ss_metrics()` | focal + R callback | speed-up |
|---|---|---|---|
| 5 × 5 | 0.12 s | 16.33 s | 140× |
| 11 × 11 | 0.26 s | 26.94 s | 104× |
| 21 × 21 | 0.71 s | 86.37 s | 122× |

2000 × 2000 raster (4 × 10^6 cells):

| Task | Time |
|---|---|
| Eight window metrics, 11 × 11 | 10.2 s |
| Gi\* with z, p, FDR and cluster map, 3 × 3 | 5.2 s |
| Local Moran with 999 Monte-Carlo randomisations, 3 × 3 | 17.2 s |

For reference, the study area in Rahimi & Jung (2025, *Land* 14:244) contained
327,104 cells, roughly one twelfth of the raster benchmarked here.

## 5. Known limitations of this release

* `s10z` uses order statistics rather than local peak detection; this is
  documented and deliberate, but it means values are not directly comparable
  with `geodiv::s10z()`.
* Local Geary's c is returned without inference.
* Monte-Carlo inference currently covers local Moran with binary weights only.
* The window engine is O(w²) per cell. An incremental column-sum update would
  reduce this to O(w) for the moment-based metrics; deferred because the
  current implementation is already fast enough for the target rasters and
  correctness was the priority.
* Multiple-testing correction holds the full p-value vector in memory.
