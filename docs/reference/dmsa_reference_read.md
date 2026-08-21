# Read a reference bundle from CSV files

The format written by the bundle builder: `systems.csv`, `polarity.csv`
and optionally `anchors.csv` in one directory, plus an optional
`manifest.dcf` carrying name, version and notes.

## Usage

``` r
dmsa_reference_read(dir, anchor_method = "user")
```

## Arguments

- dir:

  Directory holding the files.

- anchor_method:

  Passed to
  [`dmsa_reference()`](https://teindor.github.io/dmsa/reference/dmsa_reference.md)
  when no manifest is present.

## Value

A `dmsa_reference`.

## Examples

``` r
## A bundle on disk is systems.csv, polarity.csv, anchors.csv and an optional
## manifest.dcf carrying the provenance, all in one directory.
sys <- data.frame(gene = c("CRH", "POMC", "NR3C1"), system_id = "1",
                  system = "HPA axis")
pol <- data.frame(gene = sys$gene, system_id = "1", w_g = c(1, 1, -1))
d <- file.path(tempdir(), "hpa_bundle")
dmsa_reference_write(dmsa_reference(sys, pol, anchor_method = "curated",
                                    name = "hpa-toy", version = "1"), d)
list.files(d)
#> [1] "manifest.dcf" "polarity.csv" "systems.csv" 
dmsa_reference_read(d)
#> dmsa reference: hpa-toy (1)
#>   3 genes across 1 systems  (no module layer)
#>   polarity: 3 gene-system pairs, signed +-1 (2 activating, 1 braking, 0 off-axis)
#>   anchors: curated
#>   note: 
```
