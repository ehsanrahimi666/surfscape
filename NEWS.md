# surfscape 0.2.0

The decomposition release.

* `ss_decompose()`: three-way additive decomposition of the change in any
  continuous landscape metric into global compositional, local compositional
  and within-neighbourhood configurational components. Additivity is exact.
* `ss_sensitivity()`: permutation audit establishing whether a metric can
  respond to spatial arrangement, and how much it does in a given landscape.
* `ss_sensitivity_fun()`: the same audit for arbitrary metric functions, so
  implementations from `geodiv`, `GLCMTextures`, `rasterdiv` or the user can be
  tested on identical footing.
* `ss_scale_profile()`: scale response of the decomposition.
* `ss_rank_transfer()`: counterfactual surfaces used by the decomposition.
* `ss_extract()`: decomposed change in buffers around field sites.
* `ss_metrics()`: 22 built-in metrics including within-window Moran's I,
  grey-scale gliding-box lacunarity, Rao's quadratic entropy and six GLCM
  textures, all evaluated from a single window kernel.
* Example data: Landsat NDVI for 2013 and 2023, western Isfahan, Iran.

# surfscape 0.1.0

* Initial engine: moving-window surface metrics and raster-native local
  indicators of spatial association with inference.
