# Write a cascade template and the instructions for building one

Write a cascade template and the instructions for building one

## Usage

``` r
dmsa_sets_template(path = "dmsa_sets_template.csv", example = TRUE)
```

## Arguments

- path:

  Where to write the template CSV. The instructions are written
  alongside it as a `.md` file with the same stem.

- example:

  Include two filled example rows. Default `TRUE`.

## Value

The paths written, invisibly.

## Examples

``` r
p <- tempfile(fileext = ".csv")
dmsa_sets_template(p)
#> wrote /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpEgjJxj/file60485423308c.csv and /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpEgjJxj/file60485423308c.md

## the .md alongside it states the rules the validator enforces
basename(c(p, sub("\\.csv$", ".md", p)))
#> [1] "file60485423308c.csv" "file60485423308c.md" 
utils::read.csv(p)[, c("system_short", "module_id", "gene", "cpg")]
#>   system_short module_id  gene        cpg
#> 1          hpa       1.1 NR3C1 cg01234567
#> 2          hpa       1.1 FKBP5 cg07654321
```
