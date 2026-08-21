# Local false discovery rate by central matching

Efron-style two-group model on z-values: the null is fitted to the
central mass, the alternative is a wider normal, and `pi0` comes from
the central proportion. Returns one local fdr per unit.

## Usage

``` r
dmsa_lfdr(z, prior_odds_multiplier = 1, min_n = 50L)
```

## Arguments

- z:

  numeric vector of unit-level z statistics (null approximately standard
  normal after calibration).

- prior_odds_multiplier:

  Multiplies the estimated `pi1`, bounded to `[0.01, 0.9]`. Values above
  1 encode evidence inherited from a parent.

- min_n:

  Below this many units the estimate is not attempted and `NULL` is
  returned, signalling the caller to fall back to a frequentist rule.

## Value

numeric vector of local fdr values, or NULL if not estimable.

## Examples

``` r
set.seed(1)
## 300 null genes and 30 carrying a real aligned effect
z <- c(rnorm(300), rnorm(30, mean = 4))
lfdr <- dmsa_lfdr(z)
round(quantile(lfdr, c(0, .5, 1)), 3)
#>    0%   50%  100% 
#> 0.000 0.999 0.999 
mean(lfdr[301:330]) < mean(lfdr[1:300])
#> [1] TRUE
## too few units to fit the two-group model: NULL asks the caller to
## fall back to a frequentist rule rather than guess
is.null(dmsa_lfdr(rnorm(10)))
#> [1] TRUE
```
