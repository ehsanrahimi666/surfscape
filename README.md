# surfscape

**Continuous landscape metrics cannot see spatial configuration.**
`surfscape` proves which ones can, and decomposes structural change into the
part that comes from the value distribution and the part that comes from
spatial arrangement.

## The problem

Landscape ecology separates *composition* from *configuration*. When landscapes
are represented as continuous surfaces rather than patch mosaics, that
distinction was never carried across. Most metrics used on continuous
surfaces are functions of the **multiset** of values in the analysis
neighbourhood: rearranging every pixel inside the window leaves them
numerically unchanged.

Audited on real Landsat NDVI (Isfahan, Iran; 15 x 15 windows), across 47
metrics from three packages:

| Verdict | Metrics |
|---|---|
| **Configuration-blind** (CSI = 0, deviation ~1e-15) | `sa`, `sq`, `ssk`, `sku`, `smean`, `sp`, `sv`, `sz`, `scv`, Shannon entropy, **Rao's Q**, and the whole surface-bearing family `sbi`, `sci`, `svi`, `sk`, `spk`, `svk` |
| **Formally spatial, nearly blind in practice** | `glcm_mean` (CSI 0.0005), `glcm_variance` (0.002), `glcm_entropy` (0.007), `glcm_ASM` (0.012), `glcm_homogeneity` (0.065), `sdq` (0.066), lacunarity (0.065) |
| **Genuinely configurational** | within-window Moran (1.00), `glcm_correlation` (0.99), `ssc` (0.997), `sds` (0.963), `sdr` (0.71), `glcm_contrast` (0.28) |

The audit runs on the **published implementations** of `geodiv` and
`GLCMTextures`, not on reimplementations.

## The decomposition

For two dates, `ss_decompose()` splits the change in any metric into three
components that sum **exactly** to the total:

* `global` -- the landscape-wide value distribution shifted
* `local` -- values were redistributed among neighbourhoods
* `config` -- values were rearranged *within* neighbourhoods

For an aspatial metric the `config` component is identically zero, at every
scale, by construction. On the Isfahan landscape, 87% of the decade-scale
change in spatial autocorrelation was configurational -- and the standard
amplitude and entropy metrics detected 0.000% of it.

## Install

```r
# install.packages("remotes")
remotes::install_github("ehsanrahimi666/surfscape")
```

## Use

```r
library(terra); library(surfscape)

z1 <- rast(system.file("extdata", "ndvi_2013.tif", package = "surfscape"))
z2 <- rast(system.file("extdata", "ndvi_2023.tif", package = "surfscape"))

# which metrics can see arrangement, and how much do they in this landscape?
ss_sensitivity(z1, w = 15, n = 400, nperm = 49, seed = 1)

# audit a metric from any other package on the same footing
ss_sensitivity_fun(z1, function(m) c(sd = sd(as.numeric(m))), w = 15)

# decompose a decade of change
d <- ss_decompose(z1, z2, metrics = c("sq", "sent", "sr1", "gcon"), w = 15)
range(values(d[["sq_config"]]), na.rm = TRUE)   # exactly zero

# how the components respond to scale
ss_scale_profile(z1, z2, metrics = c("sq", "sr1"), w_seq = c(5, 11, 21, 31))

# decomposed change in buffers around field sites
# `pts` may be a data.frame with x/y columns, an sf object, or a SpatVector
set.seed(1)
pts <- data.frame(id = 1:6,
                  x  = runif(6, xmin(z1) + 600, xmax(z1) - 600),
                  y  = runif(6, ymin(z1) + 600, ymax(z1) - 600))
ss_extract(d, pts, radii = c(150, 300, 600))
```

## Functions

| Function | Purpose |
|---|---|
| `ss_sensitivity()` | permutation audit of the built-in metrics |
| `ss_sensitivity_fun()` | the same audit for any user or third-party metric |
| `ss_decompose()` | three-way decomposition of change |
| `ss_scale_profile()` | how the components respond to window size |
| `ss_rank_transfer()` | counterfactual surface: one date's arrangement, the other's values |
| `ss_metrics()` | moving-window metric maps |
| `ss_extract()` | decomposed change in buffers around points |
| `ss_metric_list()` | metric catalogue with the arrangement flag |

## Scope

`surfscape` does not reimplement other packages. Use `landscapemetrics` for
categorical metrics, `geodiv` for the full surface-metric suite,
`GLCMTextures` for textures, `rasterdiv` for diversity, `multiScaleR` for
scale-of-effect optimisation. All of them can be fed to `ss_sensitivity_fun()`.
The kernel here exists because the decomposition must evaluate metrics on
counterfactual neighbourhoods that never exist as a raster.

## Documentation

* `vignette("surfscape")` -- worked tutorial
* `manual/surfscape-manual.pdf` -- reference manual
* `manual/surfscape-tutorial.pdf` -- tutorial as PDF

## License

GPL-3.
