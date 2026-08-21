# Permutation p-value for a DMSA statistic

Permutation p-value for a DMSA statistic

## Usage

``` r
dmsa_perm_pvalue(observed, null_stats)
```

## Arguments

- observed:

  The observed statistic (e.g. z from `dmsa_test`).

- null_stats:

  Vector of the same statistic computed on B family-wise-permuted
  datasets (shuffle the exposure BY FAMILY, refit the probe models with
  the SAME estimator, re-run `dmsa_test`).

## Value

two-sided permutation p with the +1 correction (resolution floor
1/(B+1)).

## Examples

``` r
set.seed(1)
null_z <- rnorm(999)          # z from B family-wise permuted datasets
dmsa_perm_pvalue(3.1, null_z)
#> [1] 0.002
## the +1 correction puts a floor of 1/(B+1) on any permutation p
dmsa_perm_pvalue(50, null_z)
#> [1] 0.001
```
