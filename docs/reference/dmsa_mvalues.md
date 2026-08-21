# M-value transform with clamping

M-value transform with clamping

## Usage

``` r
dmsa_mvalues(beta, eps = 1e-04)
```

## Arguments

- beta:

  numeric matrix or vector of beta values in (0,1)

- eps:

  clamp bound

## Value

A numeric matrix or vector of M-values, the same shape as `beta`, on the
logit scale.

## Examples

``` r
beta <- matrix(runif(24, 0.05, 0.95), nrow = 6,
               dimnames = list(NULL, paste0("cg", 1:4)))
round(head(dmsa_mvalues(beta), 3), 2)
#>        cg1   cg2   cg3   cg4
#> [1,] -2.24 -0.02 -0.88 -3.38
#> [2,] -3.29 -1.05  1.09  2.54
#> [3,] -3.51 -1.24  2.06  0.73
## clamping is what keeps a beta of exactly 0 or 1 from becoming +/-Inf
dmsa_mvalues(c(0, 0.5, 1))
#> [1] -13.28757   0.00000  13.28757
```
