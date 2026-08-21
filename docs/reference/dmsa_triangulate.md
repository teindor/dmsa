# Three complementary lenses on the same family

Three complementary lenses on the same family

## Usage

``` r
dmsa_triangulate(
  M,
  data,
  rhs,
  term,
  units,
  alignment,
  block = NULL,
  B = 1999,
  winsor = 3,
  method = c("expected", "fixed"),
  correction = c("maxT", "minP"),
  weighting = c("combined", "reliability", "flat"),
  w_floor = 1.5,
  seed = NULL,
  ri_group = NULL
)
```

## Arguments

- M:

  Numeric matrix, samples x probes, methylation on the analysis scale
  (probes are the responses).

- data:

  data.frame of predictors, rows aligned to `M`.

- rhs:

  Right-hand-side terms, e.g. `c("anx","avo","sex_c","chip_f")` or with
  an interaction `c("anx*ace", ...)`.

- term:

  Design column to test, as `model.matrix` names it.

- units:

  Character vector, one entry per column of `M`, naming the unit that
  column belongs to (typically the gene).

- alignment:

  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result covering `M`'s columns.

- block:

  Permutation block labels, one per row.

- B:

  Permutations.

- winsor:

  MADs for winsorisation; `NULL` to disable.

- method:

  Pooling method for
  [`dmsa_test()`](https://teindor.github.io/dmsa/reference/dmsa_test.md).

- correction:

  `"maxT"` (default; keeps the information advantage of a unit with many
  concordant probes) or `"minP"`.

- weighting:

  Character. Probe weighting engine within a unit: `"combined"`
  (default) fuses the flat and reliability statistics on one shared
  permutation stream, `"flat"` weights every usable aligned probe
  equally, `"reliability"` weights each probe by its item-rest
  correlation with the rest of its unit. Weights are computed from
  methylation alone, so the permutation null is unaffected.

- w_floor:

  Numeric. Lower bound applied to reliability weights before
  normalisation, so a single poorly behaved probe cannot be driven to
  zero influence. Ignored when `weighting = "flat"`.

- seed:

  Optional integer.

- ri_group:

  Optional grouping factor for a one-way random intercept, one entry per
  row of `M` (typically the chip). When supplied with more than one
  level, the aligned response and every design column are quasi-demeaned
  by a REML-estimated shrinkage factor before any permutation - the GLS
  solution for `(1 | group)`, at a fraction of the cost of a mixed fit
  per probe. `NULL` (default) leaves the design untouched. The transform
  is applied identically under permutation, so it moves power, not
  type-I error.

## Value

data.frame, one row per unit: the raw and family-adjusted p under each
lens, the ACAT omnibus, the direction, and how many probes agreed.

## Examples

``` r
set.seed(1)
n <- 100
## half the probes raise expression when methylated, half lower it: on the
## methylation scale they cancel, on the expression scale they agree
d <- rep(c(1, -1), each = 3)
f <- rnorm(n)
M <- sapply(d, function(dj) dj * .6 * f + sqrt(1 - .6^2) * rnorm(n))
colnames(M) <- paste0("cg", 1:6)
dat <- data.frame(y = .5 * f + rnorm(n), cv = rnorm(n))
al <- dmsa_align(data.frame(cpg = colnames(M), d = d,
                            p_plus = ifelse(d > 0, .9, .1)),
                 genes = rep("OXTR", 6), level = "gene")

res <- dmsa_triangulate(M, dat, rhs = c("y", "cv"), term = "y",
                        units = rep("OXTR", 6), alignment = al, B = 99, seed = 1)
res[, c("unit", "n_probes", "concordance", "direction", "p_omnibus")]
#>   unit n_probes concordance direction p_omnibus
#> 1 OXTR        6           1         1      0.01
```
