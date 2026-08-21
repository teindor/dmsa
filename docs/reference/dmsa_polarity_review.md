# Rows of a polarity table that need a human decision

Three kinds of row should never be used without someone looking at them:
a sign the curator and a database disagree about, a sign no evidence
could settle, and a non-zero sign held at low confidence. This returns
them.

## Usage

``` r
dmsa_polarity_review(
  x = "alpha",
  which = c("all", "disagreement", "unresolved", "low_confidence")
)
```

## Arguments

- x:

  A `dmsa_polarity`, or anything
  [`dmsa_polarity()`](https://teindor.github.io/dmsa/reference/dmsa_polarity.md)
  accepts - including a `dmsa_selection`, which scopes the review to
  that selection's genes.

- which:

  `"all"`, `"disagreement"` (a database contradiction or a framing
  choice only the analyst can settle), `"unresolved"` or
  `"low_confidence"`.

## Value

data.frame of class `dmsa_polarity_review`.

## Examples

``` r
## The three kinds of row nobody should use unseen: a database contradiction,
## a sign no evidence settles, and a non-zero sign held at low confidence.
flagged <- dmsa_polarity_review(dmsa_select(systems = "hpa"))
table(flagged$review_flag)
#> 
#>     low_confidence_signed unresolved_needs_evidence 
#>                         7                         1 
## a disagreement is a framing call the analyst has to make, not a bug
dmsa_polarity_review(which = "disagreement")$gene
#> [1] "FOXL2"   "ADORA2A" "ADORA2B" "REST"    "PHOX2A"  "PHOX2B"  "EP300"  
#> [8] "FOXP2"   "PAX8"   
```
