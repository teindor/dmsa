# Write a template script that fetches probe coordinates from Bioconductor

The heavyweight route, kept for people who already have the annotation
packages installed.
[`dmsa_probe_coords`](https://teindor.github.io/dmsa/reference/dmsa_probe_coords.md)
does the same job from a plain-text manifest with base R only, and
[`dmsa_plot_locus`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
needs no coordinates at all.

## Usage

``` r
dmsa_probe_annotation_template(
  path = "get_probe_coords.R",
  array = c("EPIC", "450K")
)
```

## Arguments

- path:

  Where to write the script.

- array:

  `"EPIC"` or `"450K"`.

## Value

`path`, invisibly.

## Details

Run the written script once; it produces a csv with
`probe, chr, pos, gene_region, island` that
[`dmsa_plot_locus()`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
consumes.

## See also

[`dmsa_probe_coords`](https://teindor.github.io/dmsa/reference/dmsa_probe_coords.md)

## Examples

``` r
f <- tempfile(fileext = ".R")
dmsa_probe_annotation_template(f, array = "EPIC")
#> wrote /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpEgjJxj/file6048ae4978f.R - run it once, then join probe_coords.csv to your probes. Lighter: dmsa_probe_coords(), or skip coordinates entirely.
## the script is not run here: it installs and queries the minfi annotation
## packages. dmsa_probe_coords() does the same job from a plain manifest.
head(readLines(f), 5)
#> [1] "# Probe coordinates for dmsa_plot_locus(). Run once."                                 
#> [2] "# Lighter alternatives: dmsa_probe_coords(), or no coordinates at all."               
#> [3] "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
#> [4] "BiocManager::install(c('minfi', 'IlluminaHumanMethylationEPICanno.ilm10b4.hg19'))"    
#> [5] "library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19); library(minfi)"               
```
