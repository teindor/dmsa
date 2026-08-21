# Which public resources DMSA can draft polarity from, and what each is good for

Which public resources DMSA can draft polarity from, and what each is
good for

## Usage

``` r
dmsa_polarity_sources()
```

## Value

data.frame of resource name, endpoint template and caveat.

## Examples

``` r
src <- dmsa_polarity_sources()
src[c("source", "resource")]
#>     source      resource
#> 1       go Gene Ontology
#> 2 omnipath      OmniPath
#> 3   trrust     TRRUST v2
#> 4   signor    SIGNOR 3.0
## the trap worth knowing about: two of OmniPath's inputs record an
## activating sign by default when the direction is actually unknown
src$caveat[src$source == "omnipath"]
#> [1] "signed pairwise edges; consensus_stimulation at curation_effort 1 is often a default, not a finding"
```
