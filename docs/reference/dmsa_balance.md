# Per-system coverage and polarity balance

The honest-reporting companion: how many probes each system can actually
align, on which side of the activation axis. A system whose brake genes
have no usable probes silently degrades to an activation-gene test -
this table is how that is caught and reported.

## Usage

``` r
dmsa_balance(alignment)
```

## Arguments

- alignment:

  result of dmsa_align(level = "system")

## Value

one-row data.frame per call

## Examples

``` r
dcall <- data.frame(cpg = paste0("p", 1:4), d = c(1, -1, 1, NA),
                    p_plus = c(0.9, 0.1, 0.8, NA))
genes <- c("CRH", "POMC", "AVPR2", "NR3C1")
pol <- data.frame(gene = genes, w_g = c(1, 1, 0, -1))
al <- dmsa_align(dcall, genes, level = "system", polarity = pol)
# the only brake gene abstained, so this "system" test is really an
# activation-gene test - which is what the balance table makes visible
dmsa_balance(al)
#>   probes called usable activation_probes brake_probes off_axis no_direction
#> 1      4      3      2                 2            0        1            1
#>   activation_genes brake_genes
#> 1                2           0
```
