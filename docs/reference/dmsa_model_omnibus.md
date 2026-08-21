# Omnibus across the score projections

Refits the same model once per projection, substituting each into the
column named `score_col`, and combines the p-values with the aggregated
Cauchy test. ACAT is valid under arbitrary dependence, so hedging across
projections costs no multiplicity - unlike a Bonferroni over arms, which
is what the old `min(1, 2*min(p))` gate amounted to.

## Usage

``` r
dmsa_model_omnibus(
  scores,
  formula,
  data,
  term,
  block = NULL,
  score_col = "S",
  B = 1999,
  seed = NULL
)
```

## Arguments

- scores:

  data.frame from
  [`dmsa_scores()`](https://teindor.github.io/dmsa/reference/dmsa_scores.md).

- formula:

  Formula referring to the score by `score_col`.

- data, block, B, seed:

  As in
  [`dmsa_model()`](https://teindor.github.io/dmsa/reference/dmsa_model.md).

- term:

  Term to test.

- score_col:

  Column name the formula uses for the score.

## Value

list: `p_by_flavour`, the omnibus `p`, and the fits.

## Examples

``` r
set.seed(1)
n <- 80
M <- matrix(rnorm(n * 6), n, 6)
al <- dmsa_align(data.frame(cpg = paste0("cg", 1:6), d = c(1, 1, -1, 1, -1, 1),
                            p_plus = rep(.9, 6)),
                 genes = rep("OXTR", 6), level = "gene")
sc <- dmsa_scores(M, al)
d <- data.frame(age = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
d$y <- 0.6 * sc$aligned + rnorm(n)
## one p per projection, combined by ACAT - hedging costs no multiplicity
o <- dmsa_model_omnibus(sc, y ~ S + age, d, term = "S", block = d$cID,
                        B = 99, seed = 1)
round(c(o$p_by_flavour, omnibus = o$p), 3)
#> aligned    mean     pc1   blind omnibus 
#>   0.010   0.010   0.330   0.250   0.019 
```
