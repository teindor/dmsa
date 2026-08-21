# Write a reference bundle to a directory

The inverse of
[`dmsa_reference_read()`](https://teindor.github.io/dmsa/reference/dmsa_reference_read.md):
useful for saving a bundle built from a one-file csv, or a generated
bundle you have edited.

## Usage

``` r
dmsa_reference_write(reference, dir)
```

## Arguments

- reference:

  A `dmsa_reference`.

- dir:

  Directory to create and write into.

## Value

`dir`, invisibly.

## Examples

``` r
f <- tempfile(fileext = ".csv")
dmsa_reference_template(f)
#> template written to /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpEgjJxj/file6048692f01e1.csv
#>   system + gene are required; module, w_g, anchor, role, confidence, source are optional.
#>   w_g: +1 raises system activation, -1 opposes it, 0 off-axis; continuous values in [-1, 1] are allowed.
#>   anchor: TRUE for the genes that DEFINE activation for that system.
ref <- dmsa_reference_csv(f, quiet = TRUE)

d <- tempfile()
dmsa_reference_write(ref, d)
list.files(d)
#> [1] "anchors.csv"  "manifest.dcf" "polarity.csv" "systems.csv" 
nrow(dmsa_reference_read(d)$systems)
#> [1] 8
```
