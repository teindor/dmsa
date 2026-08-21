# Validate a candidate polarity table

Validate a candidate polarity table

## Usage

``` r
dmsa_polarity_check(x = "alpha", sets = NULL, verbose = TRUE)
```

## Arguments

- x:

  Path, `data.frame`, or `dmsa_polarity`.

- sets:

  Optional cascade to check coverage against.

- verbose:

  Print the report. Default `TRUE`.

## Value

Invisibly, a list with `ok` and one entry per check.

## Examples

``` r
## An anchor gene is what defines "more activation" for its system, so an
## anchor that brakes its own system is a contradiction, not a small error.
cand <- data.frame(system_id = "HPA", gene = c("CRH", "POMC", "NR3C1"),
                   w_g = c(-1, 1, -1), anchor = c(TRUE, FALSE, FALSE),
                   w_g_source = "curated")
chk <- dmsa_polarity_check(cand)
#> polarity: 3 rows | 1 systems | 0% heuristic | 0 flagged
#>   [ok] w_g in [-1, 1]                               
#>   [ok] no duplicate gene within a system            
#>   [WARN] every anchor has w_g > 0                     offenders: CRH
#>   [ok] every system has at least one anchor         
#>   [ok] no system is entirely one-sided              
#>   [ok] no role contradicts its own sign             
#>   [ok] role vocabulary is closed                    
#>   [ok] no module signed toward its own sub-process  
#>   [ok] source stated for every row                  0 row(s) without w_g_source
#>   [ok] most signs rest on evidence, not heuristics  0% heuristic
#>   -> 1 check(s) raised a warning; none is fatal, but each changes how a system-level sign should be read
chk$ok
#> [1] FALSE
dmsa_polarity_check("alpha", verbose = FALSE)$ok
#> [1] FALSE
```
