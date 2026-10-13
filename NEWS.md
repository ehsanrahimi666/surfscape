# surfscape 0.3.0

* `ss_benchmark()`: tests any continuous landscape metric, built-in or from
  another package, against controlled changes of composition (uniform shift,
  uniform stretch, rank-preserving change of distribution shape, imposed
  distribution of another surface) and configuration (random spatial
  permutation inside windows, landscape-wide rearrangement), and classifies it
  as compositional, configurational or mixed. Includes a per-window
  permutation test, standardised responses, fixed- and own-range binning
  protocols, and evaluation of the whole landscape as one window.
* `ss_benchmark_combine()`: pools results computed on separate windows or
  replicate landscapes.
* `ss_perturb()`: the controlled changes used by the benchmark.
* `ss_simulate()`: Gaussian random field landscapes by circulant embedding,
  optionally with any value distribution.
* `ss_decompose_landscape()`: compositional and configurational components of
  the change in landscape-level metrics, for built-in metrics or any metric
  function.
* `ss_scale_profile()`: `share_config` is now the mean absolute
  configurational component divided by the sum of the mean absolute values of
  the three components; correlations are `NA` when a component is constant.
* `ss_metric_list()`: new `binned` column; corrected description of `sr1`
  (lag-1 autocorrelation from the semivariance between rook neighbours).
* Documentation of `ss_decompose()` corrected: the limiting cases are stated
  for the components to which they apply.
* Removed `ss_sensitivity()` and `ss_sensitivity_fun()`, superseded by
  `ss_benchmark()`.
* `paper/` (not part of the installed package): R scripts that reproduce the
  analyses, tables and figures of the accompanying article.

# surfscape 0.2.0

* `ss_decompose()`, `ss_scale_profile()`, `ss_rank_transfer()`,
  `ss_extract()` and 22 built-in metrics computed by a shared C++ window
  kernel.
* Example data: Landsat NDVI for 2013 and 2023, western Isfahan, Iran.

# surfscape 0.1.0

* Initial moving-window surface metrics.
