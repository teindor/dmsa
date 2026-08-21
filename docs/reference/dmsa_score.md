# Subject-level aligned tone score

Weighted mean of within-probe standardised M-values, each weighted by
its aligned multiplier (`s_j` for fixed-sign, `E[s_j]` for
expected-sign), then standardised. Positive = the set reads as higher
system activation tone.

## Usage

``` r
dmsa_score(
  data,
  probes,
  alignment,
  method = c("expected", "fixed"),
  beta_input = TRUE
)
```

## Arguments

- data:

  data.frame.

- probes:

  methylation columns, in the row order of `alignment`.

- alignment:

  a
  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result.

- method:

  `"expected"` or `"fixed"`.

- beta_input:

  TRUE if `probes` hold beta values.

## Value

numeric vector, one value per row of `data` (NA where the probes are
missing).

## Examples

``` r
set.seed(1)
n <- 80
d <- rep(c(1, -1), each = 3)      # methylation-to-expression direction
f <- rnorm(n)                     # the tone the gene actually tracks
probes <- paste0("cg", 1:6)
dat <- as.data.frame(sapply(d, function(dj) plogis(dj * .8 * f + rnorm(n))))
names(dat) <- probes

al <- dmsa_align(data.frame(cpg = probes, d = d, p_plus = ifelse(d > 0, .9, .1)),
                 genes = rep("FKBP5", 6), level = "gene")
## one standardised tone per subject, despite the probes opposing each other
dat$tone <- dmsa_score(dat, probes, al)
round(cor(dat$tone, f), 2)
#> [1] 0.89
```
