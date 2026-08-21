# Validate a candidate cascade before using it

Validate a candidate cascade before using it

## Usage

``` r
dmsa_sets_check(x, verbose = TRUE)
```

## Arguments

- x:

  Path, `data.frame`, or `dmsa_sets`.

- verbose:

  Print the report. Default `TRUE`.

## Value

Invisibly, a list with `ok` and one entry per check.

## Examples

``` r
sets <- data.frame(
  system_id = "1", system_short = "hpa", system = "HPA axis",
  module_id = "1.1", module = "Corticosteroid receptors",
  gene = c("NR3C1", "NR3C1", "FKBP5"),
  cpg = c("cg01234567", "cg07654321", "cg11223344"))

## a gene in two modules, or a duplicated row, would make the family
## ambiguous; missing module evidence is only a caveat
chk <- dmsa_sets_check(sets)
#> cascade: 3 rows | 1 systems | 1 modules | 2 genes | 3 CpGs
#>   [ok] required columns                         
#>   [ok] each module in one system                offenders: 
#>   [ok] each gene in one module                  offenders: 
#>   [ok] system_short unique                      duplicated: 
#>   [ok] module metadata constant within module   varying: 
#>   [ok] no duplicate module+gene+cpg rows        
#>   [WARN] module evidence annotated                no evidence_strength - results will print without an evidence banner
#>   -> usable; 1 caveat(s) to carry into how results read
chk$ok
#> [1] TRUE
```
