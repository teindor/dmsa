# Expand a categorical focal exposure into one directional design per contrast

DMSA pools signed per-probe effects, so a focal factor with `k > 2`
levels cannot be tested as one directional term. This builds the `k - 1`
single-contrast designs, each of which *is* directional, and states the
multiplicity correction that applies across them.

## Usage

``` r
dmsa_contrasts(
  data,
  design,
  reference = NULL,
  correct = c("holm", "BH", "bonferroni", "none")
)
```

## Arguments

- data:

  The analysis data frame.

- design:

  A `dmsa_design` whose `focal_test` names a factor (or character, or a
  numeric column with few unique values) with `k > 2` levels.

- reference:

  Level to contrast against. Defaults to the first level, which for a
  factor is R's own reference.

- correct:

  Multiplicity correction across the contrasts: `"holm"` (default),
  `"BH"`, `"bonferroni"`, or `"none"`. Recorded so the choice is visible
  rather than assumed.

## Value

An object of class `dmsa_contrasts`: a list of `dmsa_design` objects in
`$designs`, the data frame with the contrast columns added in `$data`,
and the correction in `$correct`. Feed each design to
[`dmsa_fit()`](https://teindor.github.io/dmsa/reference/dmsa_fit.md) and
pass the resulting p-values to
[`dmsa_contrast_adjust()`](https://teindor.github.io/dmsa/reference/dmsa_contrast_adjust.md).

## Examples

``` r
set.seed(1)
d <- data.frame(style = factor(sample(c("secure", "anxious", "avoidant"), 40,
                                      TRUE)),
                age = rnorm(40), cID = rep(1:20, each = 2))
des <- dmsa_design(focal = "style", fixed = "age", exchangeable = "cID")
dmsa_contrasts(d, des)
#> Warning: level(s) with fewer than 20 observations: anxious, avoidant, secure. A contrast on a thin cell is unstable and its permutation null is coarse.
#> DMSA categorical focal exposure: style
#>   levels:    anxious, avoidant, secure
#>   reference: anxious   (n = 15)
#>   2 directional contrast(s), corrected by holm:
#>     1. style__avoidant_vs_anxiousavoidant   (n = 10 vs 15)
#>     2. style__secure_vs_anxioussecure   (n = 15 vs 15)
#>   Each contrast keeps a sign. A single omnibus test over 2 df would not.
```
