# Write a template csv for a user-supplied reference

Produces a small, filled-in example with every accepted column, so the
shape is unambiguous. Overwrite the rows with your own.

## Usage

``` r
dmsa_reference_template(path, with_modules = TRUE)
```

## Arguments

- path:

  Where to write it.

- with_modules:

  Include the module column.

## Value

`path`, invisibly.

## Examples

``` r
f <- tempfile(fileext = ".csv")
dmsa_reference_template(f)
#> template written to /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpV1rBC7/file5bfa7560a470.csv
#>   system + gene are required; module, w_g, anchor, role, confidence, source are optional.
#>   w_g: +1 raises system activation, -1 opposes it, 0 off-axis; continuous values in [-1, 1] are allowed.
#>   anchor: TRUE for the genes that DEFINE activation for that system.
head(utils::read.csv(f), 3)
#>     system       module  gene w_g anchor   role confidence   source
#> 1 HPA axis hypothalamic   CRH   1   TRUE driver       high PMID:...
#> 2 HPA axis hypothalamic CRHBP  -1  FALSE  brake       high PMID:...
#> 3 HPA axis    pituitary  POMC   1  FALSE driver       high PMID:...
## the filled-in rows show every accepted column; replace them with your own
## panel, then read the result back as a reference bundle
dmsa_reference_csv(f, quiet = TRUE)
#> dmsa reference: file5bfa7560a470.csv
#>   8 genes across 2 systems, 6 modules
#>   polarity: 8 gene-system pairs, continuous (6 activating, 2 braking, 0 off-axis)
#>   anchors: user (2 genes)
```
