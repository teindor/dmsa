# Alpha panel gene -\> system map (bundled)

Read the Project Alpha gene codebook shipped with the package: 549 genes
assigned to 30 biological systems.

## Usage

``` r
alpha_gene_systems()
```

## Value

data.frame with system_id, system, gene, n_kept, usability

## Examples

``` r
gs <- alpha_gene_systems()
dim(gs)
#> [1] 549  10
length(unique(gs$system_id))
#> [1] 30
head(gs[, c("system_id", "system", "gene", "n_kept", "usability")], 4)
#>   system_id                       system  gene n_kept usability
#> 1         1 Oxytocin, vasopressin & CICR TRPM2     16      good
#> 2         1 Oxytocin, vasopressin & CICR AVPR2     13      good
#> 3         1 Oxytocin, vasopressin & CICR   OXT     11      good
#> 4         1 Oxytocin, vasopressin & CICR  NOS1     10      good
```
