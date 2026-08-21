# Outcome-free reliability weights for aligned probes

Item-rest correlation of each aligned probe with the rest of its unit,
computed from the (standardised, already winsorised) methylation matrix
only. Used internally by
[`dmsa_triangulate`](https://teindor.github.io/dmsa/reference/dmsa_triangulate.md),
[`dmsa_scores`](https://teindor.github.io/dmsa/reference/dmsa_scores.md)
and
[`dmsa_change`](https://teindor.github.io/dmsa/reference/dmsa_change.md)
when `weighting = "reliability"`.

## Usage

``` r
dmsa_relweights(
  Z,
  units,
  mlt,
  w_floor = 1.5,
  weighting = c("reliability", "flat")
)
```

## Arguments

- Z:

  Numeric matrix, samples x probes, standardised methylation (the same
  matrix the lenses see; outcome-free).

- units:

  Character/label vector, one per column of `Z`.

- mlt:

  Numeric per-probe aligned multiplier (0 for unusable probes).

- w_floor:

  Minimum leading eigenvalue of the aligned probes' correlation matrix
  for weighting to engage (independent probes give eigenvalues near 1;
  the default 1.5 asks for a shared factor stronger than sampling
  noise). Below it the unit is weighted flat.

- weighting:

  `"reliability"` (default) or `"flat"` (returns all ones).

## Value

Numeric vector of per-probe weights, length `ncol(Z)`, mean 1 within
each weighted unit and 1 everywhere a guard fired.

## Examples

``` r
set.seed(2)
n <- 120
f <- rnorm(n)
## three probes report the gene's shared axis well, three barely at all
lam <- c(.8, .8, .8, .1, .1, .1)
Z <- scale(sapply(lam, function(l) l * f + sqrt(1 - l^2) * rnorm(n)))

## weights come from methylation alone, so the permutation null is untouched
round(dmsa_relweights(Z, units = rep("NR3C1", 6), mlt = rep(1, 6)), 2)
#> [1] 1.93 1.79 1.71 0.07 0.27 0.23
```
