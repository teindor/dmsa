# Polarity table for one system, in the form dmsa_align() wants

Polarity table for one system, in the form dmsa_align() wants

## Usage

``` r
dmsa_polarity_for(reference, system_id)
```

## Arguments

- reference:

  A `dmsa_reference`.

- system_id:

  The system to extract.

## Value

data.frame with `gene` and `w_g`, or NULL if the bundle carries no
polarity.

## Examples

``` r
f <- tempfile(fileext = ".csv")
dmsa_reference_template(f)
#> template written to /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpV1rBC7/file5bfa48f6003c.csv
#>   system + gene are required; module, w_g, anchor, role, confidence, source are optional.
#>   w_g: +1 raises system activation, -1 opposes it, 0 off-axis; continuous values in [-1, 1] are allowed.
#>   anchor: TRUE for the genes that DEFINE activation for that system.
ref <- dmsa_reference_csv(f, quiet = TRUE)
## w_g is the gene's sign toward its system's activation tone: +1 raises it,
## -1 opposes it. dmsa_align() multiplies w_g by each CpG's
## methylation-to-expression direction to get the probe's aligned sign.
dmsa_polarity_for(ref, "HPA axis")
#>    gene w_g
#> 1   CRH   1
#> 2 CRHBP  -1
#> 3  POMC   1
#> 4  MC2R   1
#> 5 NR3C1  -1
```
