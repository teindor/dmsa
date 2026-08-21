# Gated hierarchical DMSA: one score, one model, level-local correction

Gated hierarchical DMSA: one score, one model, level-local correction

## Usage

``` r
dmsa_gate(
  data,
  M,
  map,
  alignment,
  formula,
  term,
  roots,
  block = NULL,
  alpha = 0.05,
  B = 1999,
  winsor = 3,
  gate = c("both", "sparse", "dense"),
  sparse_reach = "children",
  seed = NULL
)
```

## Arguments

- data:

  data.frame with outcome, moderators and covariates.

- M:

  Numeric matrix, samples x probes, rows aligned to `data`.

- map:

  data.frame, one row per column of `M`, with the level columns ordered
  outermost first, e.g. `system`, `module`, `gene`, `probe`. The
  innermost column is the probe level.

- alignment:

  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result for `M`'s columns in order, or a numeric multiplier vector.

- formula:

  Upper-level model formula, using the literal name `S` where the unit's
  score belongs.

- term:

  Coefficient to test, as `model.matrix` names it - e.g. `"S"` or
  `"S:ACE"`.

- roots:

  Units at the outermost level that were PRE-REGISTERED. Their count is
  the first family's size; one root means no correction at the top.

- block:

  Permutation block labels (family/cluster id), one per row.

- alpha:

  Level for every family, propagated to children that pass.

- B:

  Permutations.

- winsor:

  Passed to
  [`dmsa_scores()`](https://teindor.github.io/dmsa/reference/dmsa_scores.md).

- gate:

  Logical. If `TRUE` (default) a level is entered only after the level
  above it has produced a surviving unit, which is the level-local
  cascade; if `FALSE` every declared level is tested independently.

- sparse_reach:

  Integer. Maximum number of units carried forward from a level into the
  level below it, used to keep a wide family from expanding the cascade
  without bound.

- seed:

  Optional integer.

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
M <- matrix(rnorm(n * 6), n, 6, dimnames = list(NULL, probes))
d <- data.frame(age = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
d$y <- as.numeric(scale(M %*% (2 * al$p_s_plus - 1))) + rnorm(n)
## the two genes are the pre-registered family; a probe is corrected
## inside its own gene only, never against all six
dmsa_gate(d, M, map, al, y ~ S + age, term = "S", roots = c("OXTR", "AVP"),
          block = d$cID, B = 99, seed = 1)
#> Gated DMSA - one score and one model per unit; correction within level and family only
#>   term 'S'   alpha = 0.05   B = 99
#>   each gate = calibrated min(dense score model, minP over children)
#> 
#>   --- GENE ---
#>    OXTR                           fam (pre-registered)        ( 2)  own 0.0100  gate 0.0100  adj 0.0100  PASS
#>    AVP                            fam (pre-registered)        ( 2)  own 0.0100  gate 0.0100  adj 0.0100  PASS
#> 
#>   --- PROBE ---
#>    cg4                            fam AVP                     ( 3)  own 0.0100  gate 0.0100  adj 0.0100  PASS
#>    cg3                            fam OXTR                    ( 3)  own 0.0200  gate 0.0200  adj 0.0300  PASS
#>    cg5                            fam AVP                     ( 3)  own 0.0400  gate 0.0400  adj 0.0900  
#>    cg2                            fam OXTR                    ( 3)  own 0.0700  gate 0.0700  adj 0.1800  
#>    cg6                            fam AVP                     ( 3)  own 0.1300  gate 0.1300  adj 0.3100  
#>    cg1                            fam OXTR                    ( 3)  own 0.1600  gate 0.1600  adj 0.3900  
#> 
```
