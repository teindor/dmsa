# Aligned set scores, and the alternatives worth comparing against

The aligned projection is
[`dmsa_score()`](https://teindor.github.io/dmsa/reference/dmsa_score.md)'s,
restated for a matrix input and with optional winsorisation. The other
three are what competing methods implicitly use, returned so that the
comparison can be run rather than argued.

## Usage

``` r
dmsa_scores(
  M,
  alignment,
  method = c("expected", "fixed"),
  winsor = 3,
  flavours = c("aligned", "mean", "pc1", "blind"),
  weighting = c("combined", "reliability", "flat"),
  w_floor = 1.5
)
```

## Arguments

- M:

  Numeric matrix, samples x probes, on the analysis scale (M-values).
  Columns are standardised internally.

- alignment:

  A
  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result for the same probes in the same order, or a numeric vector of
  multipliers.

- method:

  `"expected"` uses `2*p_s_plus - 1`; `"fixed"` uses `s`.

- winsor:

  Winsorise each probe at this many MADs first; `NULL` disables. 3 is
  the default because it won the positive-control estimator grid -
  outliers otherwise dominate the composite.

- flavours:

  Which projections to return. `aligned` is DMSA's;

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
  zero influence. Ignored when `weighting = "flat"`. `mean` is the
  unweighted mean z that "averaging the gene set" means; `pc1` is the
  regional principal component (DMRPC-style), oriented by the alignment
  so its sign stays interpretable; `blind` is mean `|z|`, the sign-blind
  control.

## Value

data.frame of standardised scores, one column per flavour.

## See also

[`dmsa_score`](https://teindor.github.io/dmsa/reference/dmsa_score.md)
for the single aligned vector on a data.frame,
[`dmsa_model`](https://teindor.github.io/dmsa/reference/dmsa_model.md)
for the upper-level model.

## Examples

``` r
if (FALSE) { # \dontrun{
al <- dmsa_align(dir_table, genes = gene, level = "gene")
S  <- dmsa_scores(M[, probes], al)
dat$oxytocin <- S$aligned
} # }
```
