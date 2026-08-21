# Gated hierarchical DMSA with level-specific inference

Gated hierarchical DMSA with level-specific inference

## Usage

``` r
dmsa_cascade(
  b,
  se,
  alignment,
  tree,
  nulls,
  engine = NULL,
  gate = c("both", "dense", "sparse"),
  q = 0.05,
  method = c("expected", "fixed"),
  borrow = 4
)
```

## Arguments

- b, se:

  Probe-level coefficient and standard error, one per probe.

- alignment:

  A
  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result covering the same probes in the same order (supplies the
  aligned multipliers).

- tree:

  data.frame with one row per probe and one column per level, ordered
  outermost first, e.g. `data.frame(system=, module=, gene=)`. Probes
  are the implicit innermost level.

- nulls:

  Named list of null statistic vectors for calibration, one entry per
  level named as in `tree`, plus `"sparse"` entries named
  `paste0(level, "_sparse")` where a sparse statistic is used. Generate
  these once per design with
  [`dmsa_cascade_null()`](https://teindor.github.io/dmsa/reference/dmsa_cascade_null.md).

- engine:

  Named character vector giving `"freq"` or `"eb"` per level (plus
  `"probe"`). Default is "freq" everywhere: in a gated tree the
  surviving unit counts are too small for an empirical-Bayes mixture to
  be estimable, and validation showed EB neutral-to-worse there. EB
  belongs on an ungated level with many units.

- gate:

  `"both"` (default), `"dense"` or `"sparse"` - the statistic used at
  frequentist levels. `"both"` takes the smaller of the two calibrated
  p-values with a factor-2 correction.

- q:

  Target error rate at each level (BH for `"freq"`, Bayesian FDR for
  `"eb"`).

- method:

  Pooling method passed to
  [`dmsa_test()`](https://teindor.github.io/dmsa/reference/dmsa_test.md).

- borrow:

  Upper bound on the prior-odds multiplier passed from a parent to its
  children at empirical-Bayes levels. 1 disables borrowing.

## Value

An object of class `dmsa_cascade`.

## Examples

``` r
set.seed(1)
tree <- expand.grid(probe = 1:3, gene = 1:2, system = 1:4)
tree <- tree[, c("system", "gene")]      # outermost level first
P <- nrow(tree); n <- 60; x <- rnorm(n)
M <- matrix(rnorm(n * P), n, P) + outer(x, 0.9 * (tree$system == 1))
al <- dmsa_align(data.frame(cpg = paste0("p", 1:P), d = 1, p_plus = 0.95),
                 genes = paste0("g", tree$gene), level = "gene")
fit <- lm(M ~ x); b <- coef(fit)["x", ]
se <- sqrt(colSums(residuals(fit)^2) / (n - 2) / sum((x - mean(x))^2))
nulls <- dmsa_cascade_null(b, se, al, tree, exposure = x, M = M, B = 49)
cs <- dmsa_cascade(b, se, al, tree, nulls)
cs$tables$system
#>    level   parent unit         z  sparse_z       stat selected n_probes engine
#> 1 system __root__    1 15.689438 12.451816 0.01015228     TRUE        6   freq
#> 2 system __root__    2 -1.072006  1.046073 0.69035533    FALSE        6   freq
#> 3 system __root__    3  2.959446  2.896919 0.08121827    FALSE        6   freq
#> 4 system __root__    4 -2.091301  2.839274 0.09137056    FALSE        6   freq
```
