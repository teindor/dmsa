# Two-wave change-score DMSA: time x methylation x factor

Tests Delta-S x exposure (x mod2) per unit through the composite
machinery
([`dmsa_model`](https://teindor.github.io/dmsa/reference/dmsa_model.md),
Freedman-Lane, family blocks), with maxT across the units of the
level-local family on one shared permutation stream.

## Usage

``` r
dmsa_change(
  M_before,
  M_after,
  data,
  outcome = NULL,
  exposure,
  mod2 = NULL,
  units,
  alignment,
  covariates = NULL,
  baseline = TRUE,
  role = c("predictor", "outcome"),
  block = NULL,
  B = 1999,
  correction = c("maxT", "minP"),
  weighting = c("combined", "reliability", "flat"),
  w_floor = 1.5,
  seed = 1
)
```

## Arguments

- M_before, M_after:

  numeric matrices, same probes in the same order (rows = the same
  people in the same order).

- data:

  data.frame with outcome, exposure, moderators and covariates.

- outcome:

  outcome column (role = "predictor").

- exposure:

  the moderating factor E; with two waves this is the \`mod\` of the
  time x S x E three-way (Delta carries the time axis).

- mod2:

  optional second moderator: tests Delta-S x E x mod2 (moderators
  mean-centered, all lower-order terms in the model, highest-order
  tested).

- units:

  unit label per probe column (genes, modules, or one system).

- alignment:

  alignment as in
  [`dmsa_triangulate`](https://teindor.github.io/dmsa/reference/dmsa_triangulate.md).

- covariates:

  character vector of covariate columns. Supply the CHANGE-SCALE cell
  composition here; a message reminds you if empty.

- baseline:

  include the unit's before-wave score as a covariate.

- role:

  "predictor": outcome ~ dS \* E (\* mod2) + S0 + covariates. "outcome":
  dS ~ E (\* mod2) + covariates (exposure predicts the change).

- block:

  exchangeable blocks (e.g. couple id).

- B, correction, seed:

  as elsewhere; correction applies across the units

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
  zero influence. Ignored when `weighting = "flat"`. of the family (maxT
  default).

## Value

data.frame of class "dmsa_change": unit, n, n_probes, b, t, p, p_adj;
attributes term, role, formula.

## Examples

``` r
set.seed(5)
n <- 80; g <- rep(c("g1", "g2"), each = 4); P <- length(g)
d <- rep(c(1, -1), length.out = P); E <- rnorm(n); y <- rnorm(n)
M0 <- matrix(rnorm(n * P), n, P)
M1 <- 0.6 * M0 + matrix(rnorm(n * P), n, P)
# plant a change in g1 that tracks y x E through the alignment
M1[, g == "g1"] <- M1[, g == "g1"] + outer(0.5 * y * E, d[g == "g1"])
al <- dmsa_align(data.frame(cpg = paste0("cg", 1:P), d = d,
                            p_plus = ifelse(d > 0, 0.9, 0.1)), genes = g)
dat <- data.frame(y = y, E = E, dEpi = rnorm(n))
dmsa_change(M0, M1, dat, outcome = "y", exposure = "E", units = g,
            alignment = al, covariates = "dEpi", B = 99, seed = 1)
#> DMSA change-score test - dS x E (time x S x E, 2-wave identity) 
#>   family of 2 unit(s), maxT across the family
#>   g1                           n   80  probes   4  b +0.507  t +6.15  p 0.0100  adj 0.0200  *
#>   g2                           n   80  probes   4  b -0.173  t -1.37  p 0.1600  adj 0.2900
```
