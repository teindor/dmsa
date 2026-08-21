# Build the aligned sign s_j for a set of probes

Build the aligned sign s_j for a set of probes

## Usage

``` r
dmsa_align(
  direction,
  genes = NULL,
  level = c("gene", "system"),
  polarity = "alpha",
  system_id = NULL,
  missing_polarity = c("error", "zero", "drop"),
  tissue = "blood"
)
```

## Arguments

- direction:

  Either a plain character vector of CpG identifiers, in which case the
  direction calls and the gene mapping are taken from the map bundled
  with dmsa (see
  [`dmsa_directions`](https://teindor.github.io/dmsa/reference/dmsa_directions.md)),
  or a table of calls you supply yourself. A supplied table may be the
  output of
  [`cpgdirection::cpg_expression_direction()`](https://rdrr.io/pkg/cpgdirection/man/cpg_expression_direction.html)
  (needs `cpg_id` or `input`, `best_direction`, and if available
  `best_confidence`), or any data.frame with `cpg`, `d` (+1/-1) and
  optionally `p_plus` = P(d = +1).

- genes:

  The set-membership gene of each probe - which gene this probe is being
  read against. Required, and the same length as the probes, when
  `direction` is a table. Optional when `direction` is a character
  vector of probe IDs: leave it `NULL` to take every gene the bundled
  map maps each probe to, or supply one gene per probe to pin specific
  pairs.

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

- tissue:

  Which bundled direction layer to read when `direction` is a character
  vector of probe IDs. Only `"blood"` ships at present; other tissues
  live in `cpgdirection`.

## Value

data.frame: probe, gene, d, p_plus, w_g, s (= d \* w_g at system level,
d at gene level), p_s_plus = P(s = +1) (chained), usable (logical), and
reason for every non-usable probe.

## Examples

``` r
# The short way: probe IDs alone. Directions and gene mapping come from
# the bundled map, so nothing has to be obtained first.
al <- dmsa_align(c("cg00052046", "cg00176879", "cg00308631"))
#> dmsa: 3 of 3 probes carry a call in the bundled blood map (100%), 5 probe-gene pairs [map 2026.08]
al[, c("probe", "gene", "d", "s", "usable")]
#>        probe          gene  d  s usable
#> 1 cg00052046           AVP -1 -1   TRUE
#> 2 cg00052046 UBOX5;FASTKD5 -1 -1   TRUE
#> 3 cg00176879           AVP -1 -1   TRUE
#> 4 cg00308631           AVP -1 -1   TRUE
#> 5 cg00308631 UBOX5;FASTKD5 -1 -1   TRUE

# Which map produced this alignment travels with it.
attr(al, "direction_map")
#> [1] "dmsa bundled blood map v2026.08 (from cpgdirection 2.3.0, 10.5281/zenodo.22024185)"

# Pin each probe to one gene when the panel says which gene it belongs to.
dmsa_align(c("cg00052046", "cg00176879"), genes = c("AVP", "AVP"))$s
#> dmsa: 2 of 2 probes carry a call in the bundled blood map (100%), 2 probe-gene pairs [map 2026.08]
#> [1] -1 -1

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
