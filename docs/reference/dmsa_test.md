# Fixed-sign or expected-sign DMSA pooled test

Pools probe-level coefficients after alignment. Two estimators:
`"fixed"` multiplies each b_j by the hard sign s_j; expected-sign
(`"expected"`) multiplies by E\[s_j\] = 2 P(s_j = +1) - 1, the
first-order latent-sign approximation - low-tier calls contribute in
proportion to their calibrated certainty, which is why no tier filtering
is needed upstream.

## Usage

``` r
dmsa_test(b, se, alignment, method = c("expected", "fixed"), w = NULL)
```

## Arguments

- b, se:

  Probe-level coefficient and SE from YOUR model (fitted on M-values
  standardised within probe; same estimator observed and null).

- alignment:

  A
  [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  result covering the same probes, or vectors via `s` / `p_s_plus`.

- method:

  `"fixed"` or `"expected"`.

- w:

  Optional per-probe reliability weight (see
  [`dmsa_relweights`](https://teindor.github.io/dmsa/reference/dmsa_relweights.md));
  scales each aligned multiplier. `NULL` (default) or all ones is the
  flat engine.

## Value

list: estimate (pooled aligned slope), se, z, p_normal, n_used,
n_excluded, table (per-probe contributions).

## Details

The z test here is descriptive; the reportable error control is the
family-wise permutation (`dmsa_perm_pvalue`) on the same statistic.

## Examples

``` r
## Two promoter CpGs in NR3C1: methylation silences the gene (d = -1), and
## GR is the HPA brake (w_g = -1), so both probes vote for HIGHER axis tone.
calls <- data.frame(cpg = c("cg01", "cg02"), d = c(-1, -1), p_plus = c(.1, .2))
al <- dmsa_align(calls, genes = c("NR3C1", "NR3C1"), level = "system",
                 polarity = data.frame(gene = "NR3C1", w_g = -1))

tt <- dmsa_test(b = c(0.28, 0.19), se = c(0.09, 0.08), al, method = "expected")
tt[c("estimate", "se", "z", "p_normal", "n_used")]
#> $estimate
#> [1] 0.336138
#> 
#> $se
#> [1] 0.08598279
#> 
#> $z
#> [1] 3.909364
#> 
#> $p_normal
#> [1] 9.253936e-05
#> 
#> $n_used
#> [1] 2
#> 
tt$table
#>   probe  gene    b   se s p_s_plus multiplier b_aligned used
#> 1  cg01 NR3C1 0.28 0.09 1      0.9        0.8     0.224 TRUE
#> 2  cg02 NR3C1 0.19 0.08 1      0.8        0.6     0.114 TRUE
```
