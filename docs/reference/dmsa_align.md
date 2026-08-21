# Build the aligned sign s_j for a set of probes

Build the aligned sign s_j for a set of probes

## Usage

``` r
dmsa_align(
  direction,
  genes,
  level = c("gene", "system"),
  polarity = "alpha",
  system_id = NULL,
  missing_polarity = c("error", "zero", "drop")
)
```

## Arguments

- direction:

  A data.frame from
  [`cpgdirection::cpg_expression_direction()`](https://rdrr.io/pkg/cpgdirection/man/cpg_expression_direction.html)
  (needs columns `cpg_id` or `input`, `best_direction`, and if available
  `best_confidence`), OR a data.frame with columns `cpg`, `d` (+1/-1)
  and optionally `p_plus` = P(d = +1).

- genes:

  Character vector, same length as the probes in `direction`: the
  set-membership gene of each probe (the panel's annotation, e.g. from
  the column suffix).

- level:

  `"gene"` (alignment = d_j only) or `"system"` (alignment = d_j \*
  w_g).

- polarity:

  For `level = "system"`: where w_g comes from. `"alpha"` (default) uses
  the bundled Alpha polarity table. Or supply your own data.frame with
  columns `gene` and `w_g` (+1 / -1 / 0). User-supplied values OVERRIDE
  the bundled table where both cover a gene.

- system_id:

  Optional integer: restrict the bundled polarity lookup to one Alpha
  system (a gene can carry different roles in different systems).

- missing_polarity:

  What to do with a called probe whose gene has no polarity entry:
  `"error"` (default - stop and list the genes so the user can specify
  each, per protocol), `"zero"` (treat as off-axis, probe drops from the
  system pool but is reported), or `"drop"` (silently exclude).

## Value

data.frame: probe, gene, d, p_plus, w_g, s (= d \* w_g at system level,
d at gene level), p_s_plus = P(s = +1) (chained), usable (logical), and
reason for every non-usable probe.

## Examples

``` r
# d is the CpG -> expression sign, w_g the gene -> system-tone sign.
# A promoter CpG silences NR3C1 (d = -1) and GR is the HPA brake (w_g = -1),
# so the probe votes for HIGHER axis tone: s = d * w_g = +1.
g <- c("CRH", "NR3C1")
dcall <- data.frame(cpg = c("cg_crh", "cg_nr3c1"), d = c(1, -1),
                    p_plus = c(0.9, 0.1))
pol <- data.frame(gene = g, w_g = c(1, -1))
dmsa_align(dcall, genes = g, level = "gene")
#>      probe  d p_plus  gene w_g  s p_s_plus usable reason
#> 1   cg_crh  1    0.9   CRH  NA  1      0.9   TRUE       
#> 2 cg_nr3c1 -1    0.1 NR3C1  NA -1      0.1   TRUE       
dmsa_align(dcall, genes = g, level = "system", polarity = pol)
#>      probe  d p_plus  gene w_g s p_s_plus usable reason
#> 1   cg_crh  1    0.9   CRH   1 1      0.9   TRUE       
#> 2 cg_nr3c1 -1    0.1 NR3C1  -1 1      0.9   TRUE       
```
