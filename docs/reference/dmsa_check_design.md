# Validate a design against the data before fitting

Checks the things that silently ruin a methylation set analysis: a
covariate that is constant in this subsample, a covariate that proxies
the focal term, a grouping factor with too few levels to estimate, a
grouping factor aliased with the focal term, and declared independence
in data that plainly repeats.

## Usage

``` r
dmsa_check_design(design, data, strict = TRUE)
```

## Arguments

- design:

  A `dmsa_design`.

- data:

  A data.frame.

- strict:

  If TRUE (default) problems are errors; if FALSE they are warnings and
  the checks are still returned.

## Value

Invisibly, a list with `n`, `problems`, `notes` and the per-term
diagnostics.

## Examples

``` r
set.seed(1)
d <- data.frame(exposure = rnorm(80), sex_c = rep(c(-0.5, 0.5), 40),
                age = rnorm(80, 40, 5), Fib = 0,
                chip = rep(1:10, each = 8), cID = rep(1:40, each = 2))
des <- dmsa_design("exposure", c("sex_c", "age"),
                   random = c("chip", "cID"), exchangeable = "cID")
dmsa_check_design(des, d)$n
#> [1] 80
# Fib is constant in this subsample: absorbed silently by lm, caught here
bad <- dmsa_design("exposure", c("sex_c", "Fib"), random = "chip",
                   exchangeable = "cID")
dmsa_check_design(bad, d, strict = FALSE)$problems
#> Warning: design check failed:
#>   - constant in this subsample (drop it): Fib
#> [1] "constant in this subsample (drop it): Fib"
```
