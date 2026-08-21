# Alpha panel gene -\> system-activation polarity (bundled, DRAFT)

Curated w_g per gene: +1 gene product increases the system's activation
tone, -1 decreases it (e.g. NR3C1 in the HPA axis: more GR = stronger
negative feedback = LOWER axis tone), 0 = not on the activation axis
(readouts, specificity controls, ambiguous).

## Usage

``` r
alpha_polarity()
```

## Value

data.frame with system_id, system, gene, w_g_draft, role, confidence

## Details

STATUS: draft pending PI approval and per-gene citations. Every use
prints a reminder until the table is finalised.

## Examples

``` r
pol <- alpha_polarity()
head(pol[, c("system", "gene", "w_g", "role")], 4)
#>     system  gene w_g   role
#> 1 HPA axis   CRH   1 driver
#> 2 HPA axis   UCN   1 driver
#> 3 HPA axis  POMC   1 driver
#> 4 HPA axis PCSK1   1 driver
# CRH drives the HPA axis (w_g = +1); NR3C1 is its brake, so more GR means
# stronger negative feedback and LOWER axis tone (w_g = -1)
pol[pol$gene %in% c("CRH", "NR3C1") & pol$system_id == 2,
    c("gene", "w_g", "role")]
#>     gene w_g   role
#> 1    CRH   1 driver
#> 10 NR3C1  -1  brake
```
