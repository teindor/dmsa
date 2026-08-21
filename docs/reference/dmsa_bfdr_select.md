# Select units by Bayesian FDR

Largest set whose MEAN local fdr does not exceed `q`.

## Usage

``` r
dmsa_bfdr_select(lfdr, q = 0.05)
```

## Arguments

- lfdr:

  vector of local fdr values

- q:

  target Bayesian FDR

## Value

integer indices of selected units

## Examples

``` r
lfdr <- c(rep(0.001, 10), rep(0.5, 10), rep(0.99, 10))
sel <- dmsa_bfdr_select(lfdr, q = 0.05)
sel
#>  [1]  1  2  3  4  5  6  7  8  9 10 11
mean(lfdr[sel])
#> [1] 0.04636364
# nothing is selected when no unit is credible enough to fit the budget
dmsa_bfdr_select(rep(0.9, 20), q = 0.05)
#> integer(0)
```
