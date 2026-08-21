# Fit, pool and permute a DMSA set test under a declared design

Fit, pool and permute a DMSA set test under a declared design

## Usage

``` r
dmsa_fit(
  data,
  probes,
  alignment,
  design,
  method = c("expected", "fixed"),
  B = 999L,
  engine = c("auto", "lm", "lmer"),
  beta_input = TRUE,
  check = TRUE,
  seed = NULL
)
```

## Arguments

- data:

  data.frame with one row per observation.

- probes:

  Character vector of methylation columns, in the SAME order as the rows
  of `alignment`.

- alignment:

  A
  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result.

- design:

  A
  [`dmsa_design()`](https://teindor.github.io/dmsa/reference/dmsa_design.md)
  (see
  [`alpha_design()`](https://teindor.github.io/dmsa/reference/alpha_design.md)).

- method:

  `"expected"` (certainty-weighted) or `"fixed"`.

- B:

  Permutations. `0` skips the permutation (descriptive only).

- engine:

  `"auto"` reports from `lmer` when the design declares random effects
  and lme4 is installed, otherwise `lm`. The permutation null is ALWAYS
  built with the matched `lm`, and the observed `lm` statistic is the
  one compared to it - mixing estimators between observed and null is
  the classic way to manufacture significance. The agreement between the
  two is returned.

- beta_input:

  If TRUE (default) `probes` hold beta values and are converted to M;
  set FALSE if they are already M.

- check:

  Run
  [`dmsa_check_design()`](https://teindor.github.io/dmsa/reference/dmsa_check_design.md)
  first (default TRUE).

- seed:

  Optional integer for reproducible permutation.

## Value

An object of class `dmsa_fit`.

## Examples

``` r
set.seed(1)
n <- 80; K <- 6; x <- rnorm(n); d <- sample(c(-1, 1), K, TRUE)
M <- matrix(rnorm(n * K), n, K) + outer(x, 0.5 * d)
dat <- data.frame(cID = rep(1:40, each = 2), chip = rep(1:10, each = 8),
                  x = x, sex_c = rep(c(-0.5, 0.5), 40))
probes <- paste0("cg", 1:K); dat[probes] <- 1 / (1 + 2^(-M))
al <- dmsa_align(data.frame(cpg = probes, d = d,
                            p_plus = ifelse(d > 0, 0.95, 0.05)),
                 genes = rep(c("g1", "g2"), each = 3))
des <- dmsa_design("x", "sex_c", random = c("chip", "cID"),
                   exchangeable = "cID")
dmsa_fit(dat, probes, al, des, B = 99, engine = "lm", seed = 1)
#> DMSA set test (expected-sign)
#>   focal x   n = 80   probes 6 used of 6
#>   estimate 0.4768  se 0.0529  z 9.01
#>   permutation p = 0.01  (B = 99, block = cID, z_null-matched = 9.01)
#>   engine lm
```
