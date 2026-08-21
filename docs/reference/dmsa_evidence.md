# Module-level evidence and citations for a cascade, selection or frame

Every module label in the bundled cascade was checked against the
literature. This returns what that check concluded, so a module-level
finding can be reported with the strength of its own definition
attached: whether the label was retained or revised, whether the locked
membership is homogeneous, the evidence tier, and the citation keys
behind it.

## Usage

``` r
dmsa_evidence(x, which = c("all", "moderate", "flagged"))
```

## Arguments

- x:

  A `dmsa_sets`, `dmsa_selection`, or `dmsa_frame`.

- which:

  `"all"`, `"moderate"` (Moderate evidence only), or `"flagged"`
  (Moderate, heterogeneous, or measurement-defined).

## Value

data.frame of class `dmsa_evidence`, one row per module.

## Examples

``` r
# \donttest{
dmsa_evidence(dmsa_select(systems = "immune"), which = "flagged")
#> module evidence (3 modules)
#>  24.11   Moderate T-cell receptor loci and CCL2 chemotaxis   [Loz07]; [New12]
#>           status: renamed_for_mechanistic_precision 
#>  24.12   Moderate Leukocyte activation and inflammatory sign [Loz07]; [Pla18]
#>           status: renamed_from_association_label_to_process_label 
#>  24.13   Moderate CCL3 chemokine signaling                   [Loz07]; [New12]
#>           status: renamed_for_mechanistic_precision 
#>  evidence searches:
#>     https://app.undermind.ai/projects/4a4c1774-bd88-462d-889b-015f29a373b5?path=/Epigenetic%20immune%20ageing%20metabolic%20module%20evidence 
# }
```
