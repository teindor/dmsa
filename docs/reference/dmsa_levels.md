# Gated hierarchical DMSA on the pooled aligned statistic

Gated hierarchical DMSA on the pooled aligned statistic

## Usage

``` r
dmsa_levels(
  M,
  data,
  rhs,
  term,
  map,
  alignment,
  roots,
  block = NULL,
  alpha = 0.05,
  B = 1999,
  winsor = 3,
  method = c("expected", "fixed"),
  gate = c("both", "sparse", "dense"),
  family_correction = c("maxT", "minP"),
  sparse_reach = "children",
  seed = NULL
)
```

## Arguments

- M:

  Numeric matrix, samples x probes, methylation on the analysis scale.
  The probes are the RESPONSES, as in
  [`dmsa_fit()`](https://teindor.github.io/dmsa/reference/dmsa_fit.md).

- data:

  data.frame of predictors, same row order as `M`.

- rhs:

  Character vector of right-hand-side terms, e.g.
  `c("anx","avo","sex_c","chip_f")` or with an interaction
  `c("anx*aceZ", ...)`.

- term:

  The design column to pool, as `model.matrix` names it - `"anx"` or
  `"anx:aceZ"`.

- map:

  data.frame, one row per column of `M`, level columns outermost first
  (e.g. system, module, gene, probe).

- alignment:

  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result for `M`'s columns in order.

- roots:

  Pre-registered units at the outermost level.

- block:

  Permutation block labels, one per row.

- alpha, B, winsor, seed:

  As in
  [`dmsa_gate`](https://teindor.github.io/dmsa/reference/dmsa_gate.md).

- method:

  Pooling method for
  [`dmsa_test()`](https://teindor.github.io/dmsa/reference/dmsa_test.md).

- gate:

  `"both"`, `"sparse"` or `"dense"`.

- family_correction:

  Character. Family-wise correction applied inside each level-local
  family: `"maxT"` (default) Westfall-Young step-down on the strength
  scale, or `"minP"` step-down on marginal p-values. The choice is
  consequential and should be fixed in advance; see the package vignette
  and `REPRODUCE.md`. sparse arm may look straight at.

- sparse_reach:

  `"children"` or the name of a deeper level a unit's

## Value

Object of class `dmsa_gate`.

## Examples

``` r
set.seed(1)
n <- 60; probes <- paste0("cg", 1:6)
map <- data.frame(gene = rep(c("OXTR", "AVP"), each = 3), probe = probes)
al <- dmsa_align(data.frame(cpg = probes, d = rep(c(1, -1), 3),
                            p_plus = rep(c(.9, .1), 3)),
                 genes = map$gene, level = "gene")
d <- data.frame(anx = rnorm(n), age = rnorm(n),
                cID = rep(seq_len(n / 2), each = 2))
## the probes are the RESPONSES: their aligned coefficients are pooled,
## so cross-probe concordance is what carries the evidence
M <- matrix(rnorm(n * 6), n, 6) + outer(0.3 * d$anx, 2 * al$p_s_plus - 1)
dmsa_levels(M, d, rhs = c("anx", "age"), term = "anx", map = map,
            alignment = al, roots = c("OXTR", "AVP"), block = d$cID,
            B = 99, seed = 1)
#> Gated DMSA - one score and one model per unit; correction within level and family only
#>   term 'anx'   alpha = 0.05   B = 99
#>   each gate = calibrated min(dense score model, minP over children)
#> 
#>   --- GENE ---
#>    OXTR                           fam (pre-registered)        ( 2)  own 0.0100  gate 0.0100  adj 0.0100  PASS
#>    AVP                            fam (pre-registered)        ( 2)  own 0.0100  gate 0.0100  adj 0.0100  PASS
#> 
#>   --- PROBE ---
#>    cg3                            fam OXTR                    ( 3)  own 0.0100  gate 0.0100  adj 0.0200  PASS
#>    cg5                            fam AVP                     ( 3)  own 0.0100  gate 0.0100  adj 0.0800  
#>    cg2                            fam OXTR                    ( 3)  own 0.0300  gate 0.0300  adj 0.0800  
#>    cg6                            fam AVP                     ( 3)  own 0.0800  gate 0.0800  adj 0.1500  
#>    cg1                            fam OXTR                    ( 3)  own 0.0800  gate 0.0800  adj 0.2900  
#>    cg4                            fam AVP                     ( 3)  own 0.8900  gate 0.8900  adj 1.0000  
#> 
```
