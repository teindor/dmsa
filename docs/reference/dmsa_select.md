# Select a set of units from a cascade

Choose one or more systems by short name. Everything below the system -
modules, genes, probes - defaults to `"full"`, meaning every unit the
cascade assigns to the chosen system(s). Narrow a level only when the
question is genuinely narrower, and note that narrowing changes the
declared family and therefore the multiplicity toll.

## Usage

``` r
dmsa_select(
  x = "alpha",
  systems = NULL,
  modules = "full",
  genes = "full",
  probes = "full",
  columns = NULL,
  data = NULL
)
```

## Arguments

- x:

  A `dmsa_sets`, or anything
  [`dmsa_sets()`](https://teindor.github.io/dmsa/reference/dmsa_sets.md)
  accepts.

- systems:

  Short names (`c("hpa", "oxytocin")`), system ids, or full system
  names. `NULL` selects every system in the cascade. Matching is
  case-insensitive and accepts a unique prefix of a short name.

- modules, genes, probes:

  `"full"` (default) or explicit ids/symbols to restrict to. `modules`
  takes module ids such as `"2.6"`, `genes` takes symbols, `probes`
  takes CpG ids or probe ids.

- columns:

  Which cascade column holds the data-matrix column names for this
  analysis - e.g. `"col_parent_T1"`. `"auto"` picks the `col_*` column
  with the most overlap against `data`, when `data` is supplied. `NULL`
  leaves it unresolved.

- data:

  Optional data.frame or matrix used by `columns = "auto"`.

## Value

An object of class `dmsa_selection`.

## Examples

``` r
# \donttest{
sel <- dmsa_select(systems = c("hpa", "oxytocin"))
sel
#> dmsa selection from: Project Alpha 2026c (module-audited) 
#>   2 system(s) selected | modules full | genes full | probes full
#>    hpa              HPA axis & glucocorticoid signalling      8 mod    49 gene    538 cpg
#>    oxytocin         Oxytocin, vasopressin & CICR              9 mod    90 gene   1690 cpg
#>   total: 17 modules, 139 genes, 2228 CpGs
#>   module evidence: 17 High, 0 Moderate
#>   polarity: 110 signed, 29 off-axis (curated 52, database 10, literature 13, heuristic 62, none 2)
#>     35 row(s) need a decision - dmsa_polarity_review()
dmsa_select(systems = "hpa", modules = c("2.6", "2.8"))
#> dmsa selection from: Project Alpha 2026c (module-audited) 
#>   1 system(s) selected | modules restricted | genes full | probes full
#>    hpa              HPA axis & glucocorticoid signalling      2 mod    19 gene    154 cpg
#>   total: 2 modules, 19 genes, 154 CpGs
#>   module evidence: 2 High, 0 Moderate
#>   polarity: 14 signed, 5 off-axis (curated 11, database 1, heuristic 6, none 1)
#>     7 row(s) need a decision - dmsa_polarity_review()
# }
```
