# Test an aligned tone score against a non-Gaussian outcome

The subject-level arm of DMSA for outcomes that are not continuous:
binary (logistic), ordered-categorical (proportional odds),
unordered-categorical (multinomial), or counts (Poisson / negative
binomial). Optionally moderated.

## Usage

``` r
dmsa_outcome(
  data,
  outcome,
  score,
  family = c("gaussian", "binomial", "ordinal", "multinomial", "poisson"),
  moderator = NULL,
  covariates = character(),
  block = NULL,
  center = TRUE,
  B = 999L,
  seed = NULL
)
```

## Arguments

- data:

  data.frame with one row per subject.

- outcome:

  Name of the outcome column. For `"binomial"` it may be 0/1, logical,
  or a two-level factor; for `"ordinal"` an ordered factor or an
  integer-coded severity; for `"multinomial"` a factor.

- score:

  Name of the tone-score column (from
  [`dmsa_score()`](https://teindor.github.io/dmsa/reference/dmsa_score.md)).

- family:

  One of `"gaussian"`, `"binomial"`, `"ordinal"`, `"multinomial"`,
  `"poisson"`.

- moderator:

  Optional moderator column; when given, the tested term is
  `score:moderator`.

- covariates:

  Character vector of covariate columns.

- block:

  Optional column naming the exchangeable unit for permutation
  (families, repeated measures).

- center:

  Center the moderator (default TRUE).

- B:

  Permutations. 0 skips them and returns model-based inference only.

- seed:

  Optional integer.

## Value

An object of class `dmsa_outcome`.

## Details

**Direction.** A binary, count or ordinal outcome yields ONE coefficient
for the tone score, so DMSA's directional claim survives intact - the
sign of the log-odds (or log-rate, or latent-scale slope) is the
expression-aligned direction. A **multinomial** outcome yields K-1
coefficients and there is no single sign: the test becomes an omnibus
"tone predicts category" and the per-category log-odds are reported
separately. If the categories are in any sense ordered,
`family = "ordinal"` keeps the direction and is strictly preferable.

## Examples

``` r
set.seed(1)
n <- 120
d <- data.frame(S = rnorm(n), cv = rnorm(n),
                cid = rep(seq_len(n / 2), each = 2))
d$y <- rbinom(n, 1, plogis(0.9 * d$S))
## a binary outcome yields ONE coefficient, so the sign of the log-odds
## is still the expression-aligned direction
dmsa_outcome(d, outcome = "y", score = "S", family = "binomial",
             covariates = "cv", block = "cid", B = 99, seed = 2)
#> DMSA subject-level test (binomial)
#>   y ~ S   n = 120
#>   log-odds = +0.7094  se = 0.2377  |z| = 2.98  p(model) = 0.002845
#>   exp(coef) = 2.033
#>   permutation p = 0.01  (B = 99, predictor-side residual permutation)
```
