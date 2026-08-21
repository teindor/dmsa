# Reading direction from methylation: a first DMSA analysis

## The problem this package exists to solve

A CpG in a promoter and a CpG in a gene body carry opposite consequences
for expression when their methylation moves the same way. Conventional
set-level methods average them, so a gene whose probes move *coherently
on the expression scale* can cancel itself to nothing on the methylation
scale.

This is not a power problem that a larger sample fixes. At a 50/50
composition the cancellation is exact and a direction-blind test is null
by construction.

We can build exactly that situation and watch it happen.

``` r

set.seed(1)

n <- 300L                      # participants
K <- 12L                       # probes in the gene
f <- rnorm(n)                  # the latent process the gene actually tracks

# Half the probes are predicted to RAISE expression when methylated (d = +1),
# half to LOWER it (d = -1). On the expression scale they agree. On the
# methylation scale they point in opposite directions.
d <- rep(c(1, -1), each = K / 2)

M <- sapply(seq_len(K), function(j) d[j] * 0.55 * f + sqrt(1 - 0.55^2) * rnorm(n))
colnames(M) <- sprintf("cg%04d", seq_len(K))

dat <- data.frame(y = 0.35 * f + rnorm(n), cv = rnorm(n))
```

## What a direction-blind test sees

Average the probes and test the average. This is, in essence, what a
composite-score method does when it ignores direction.

``` r

blind <- summary(lm(rowMeans(M) ~ y + cv, data = dat))$coefficients["y", ]
round(blind, 4)
#>   Estimate Std. Error    t value   Pr(>|t|) 
#>     0.0037     0.0154     0.2397     0.8108
```

The estimate sits near zero and the *p*-value near 1. The signal did not
fail to exist; it was destroyed in the averaging step.

## The alignment step

[`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
takes the per-probe direction calls and turns them into a multiplier.
`d` is the methylation-to-expression direction, `p_plus` the confidence
that the direction is positive. At the gene level no polarity is
applied; at the system level each gene is additionally signed against
its system’s activation tone.

``` r

direction <- data.frame(
  cpg    = colnames(M),
  d      = d,
  p_plus = ifelse(d > 0, 0.9, 0.1)
)

al <- dmsa_align(direction, genes = rep("GENE1", K), level = "gene")
al[, c("probe", "d", "s", "p_s_plus", "usable")]
#>     probe  d  s p_s_plus usable
#> 1  cg0001  1  1      0.9   TRUE
#> 2  cg0002  1  1      0.9   TRUE
#> 3  cg0003  1  1      0.9   TRUE
#> 4  cg0004  1  1      0.9   TRUE
#> 5  cg0005  1  1      0.9   TRUE
#> 6  cg0006  1  1      0.9   TRUE
#> 7  cg0007 -1 -1      0.1   TRUE
#> 8  cg0008 -1 -1      0.1   TRUE
#> 9  cg0009 -1 -1      0.1   TRUE
#> 10 cg0010 -1 -1      0.1   TRUE
#> 11 cg0011 -1 -1      0.1   TRUE
#> 12 cg0012 -1 -1      0.1   TRUE
```

`s` is the sign each probe contributes once it is read on the expression
scale. Every probe now votes the same way, which is what the biology
said all along.

## Running the three lenses

[`dmsa_triangulate()`](https://teindor.github.io/dmsa/reference/dmsa_triangulate.md)
tests the unit through three pre-specified lenses on one shared
permutation stream: **coherence** (pooled aligned burden, which rewards
cross-probe sign concordance), **composite** (the unit’s aligned score
as a single regressor), and **diffuse** (a confidence-weighted quadratic
score). An ACAT omnibus combines them, and because ACAT is valid under
arbitrary dependence, hedging across the three costs no multiplicity.

``` r

res <- dmsa_triangulate(
  M         = M,
  data      = dat,
  rhs       = c("y", "cv"),
  term      = "y",
  units     = rep("GENE1", K),
  alignment = al,
  B         = 499,
  seed      = 1
)

res[, c("unit", "n_probes", "concordance", "direction",
        "p_coherence", "p_composite", "p_omnibus")]
#>    unit n_probes concordance direction p_coherence p_composite p_omnibus
#> 1 GENE1       12           1         1       0.002       0.002     0.002
```

The same data that read as null a moment ago now carry the signal.
Nothing about the measurements changed — only whether the aggregation
step was told which way each probe points. `concordance` of 1.00 says
every usable probe agreed once aligned, and `direction` gives the sign
of the effect on the expression scale.

Permutation *p*-values have granularity `1/(B + 1)`, so at `B = 499` the
smallest attainable value is 0.002 and the numbers above are at that
floor. Raise `B` for a finer resolution; the manuscript uses `B = 1999`.

## Which lens fires is itself a result

The three lenses fail in different ways, so their pattern is informative
rather than redundant. Coherence rewards a unit whose probes agree in
sign; composite rewards a unified shift that survives as one regressor;
diffuse rewards signal spread across many probes whichever way each
points. A finding that fires on coherence alone has a different shape
from one that fires on composite alone, and the package reports both
rather than collapsing them to a single number.

## What DMSA does not do

Directional alignment buys power against cancellation and costs power
elsewhere. Where the question is whether a *whole panel* is associated
with an outcome, and no directional structure is claimed, a quadratic
score such as `globaltest` detects associations that DMSA’s pooled score
does not — in the manuscript’s benchmark, that comparison is one DMSA
loses. Alignment is a bet that the direction calls carry information.
When they do not, the bet costs something.

## Moving to real data

Two things replace the simulation above:

- **Direction calls.** `d` and `p_plus` come from a curated
  methylation-to-expression map. The companion package `cpgdirection`
  (<https://github.com/teindor/cpgdirection>, archived at
  <https://doi.org/10.5281/zenodo.22024185>) supplies them for 3.1
  million CpG-gene pairs across four evidence classes, but any table
  with the same three columns works. `dmsa` does not require it.
- **The two-call interface.**
  [`dmsa_frame()`](https://teindor.github.io/dmsa/reference/dmsa_frame.md)
  declares the design once — outcome, covariates, permutation blocks,
  the level hierarchy — and
  [`dmsa_report()`](https://teindor.github.io/dmsa/reference/dmsa_report.md)
  runs it and emits figures, machine-readable tables at every level, and
  a written paragraph stating what was tested and what was controlled.
  See
  [`?dmsa_frame`](https://teindor.github.io/dmsa/reference/dmsa_frame.md)
  for the full argument set.

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: aarch64-apple-darwin23
#> Running under: macOS Tahoe 26.5.2
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib 
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
#> 
#> locale:
#> [1] he_IL.UTF-8/he_IL.UTF-8/he_IL.UTF-8/C/he_IL.UTF-8/he_IL.UTF-8
#> 
#> time zone: Asia/Jerusalem
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dmsa_0.99.2
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
#>  [5] xfun_0.60         cachem_1.1.0      knitr_1.51        htmltools_0.5.9  
#>  [9] rmarkdown_2.31    lifecycle_1.0.5   cli_3.6.6         sass_0.4.10      
#> [13] pkgdown_2.2.1     textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2
#> [17] compiler_4.6.1    tools_4.6.1       ragg_1.5.2        bslib_0.12.0     
#> [21] evaluate_1.0.5    yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0   
#> [25] rlang_1.3.0       fs_2.1.0          htmlwidgets_1.6.4
```
