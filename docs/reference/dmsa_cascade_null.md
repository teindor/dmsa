# Null statistics for cascade calibration

Simulates the pooled aligned statistic under the global null for each
unit shape present in the tree, so that
[`dmsa_cascade()`](https://teindor.github.io/dmsa/reference/dmsa_cascade.md)
can convert observed statistics into calibrated p-values. Uses the
supplied residual covariance structure by resampling the exposure within
blocks.

## Usage

``` r
dmsa_cascade_null(
  b,
  se,
  alignment,
  tree,
  exposure,
  M,
  block = NULL,
  B = 999L,
  method = c("expected", "fixed")
)
```

## Arguments

- b, se, alignment, tree:

  as in `dmsa_cascade`.

- exposure:

  numeric exposure vector used to fit `b`.

- M:

  matrix of probe values (rows = observations), used to refit under
  permutation.

- block:

  optional grouping vector defining exchangeable units.

- B:

  number of permutations.

- method:

  pooling method.

## Value

named list suitable for the `nulls` argument.

## Examples

``` r
set.seed(1)
tree <- expand.grid(probe = 1:3, gene = 1:2, system = 1:4)
tree <- tree[, c("system", "gene")]      # outermost level first
P <- nrow(tree); n <- 60; x <- rnorm(n)
M <- matrix(rnorm(n * P), n, P)
al <- dmsa_align(data.frame(cpg = paste0("p", 1:P), d = 1, p_plus = 0.95),
                 genes = paste0("g", tree$gene), level = "gene")
fit <- lm(M ~ x); b <- coef(fit)["x", ]
se <- sqrt(colSums(residuals(fit)^2) / (n - 2) / sum((x - mean(x))^2))
# B = 49 keeps the example quick; a real calibration wants B >= 999
nulls <- dmsa_cascade_null(b, se, al, tree, exposure = x, M = M, B = 49)
vapply(nulls, length, integer(1))
#>        system          gene system_sparse 
#>           196           392           196 
```
