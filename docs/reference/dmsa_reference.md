# Construct a DMSA reference bundle

Construct a DMSA reference bundle

## Usage

``` r
dmsa_reference(
  systems,
  polarity = NULL,
  anchors = NULL,
  anchor_method = c("curated", "graph_sink", "user", "none"),
  name = "unnamed",
  version = NA_character_,
  notes = NULL
)
```

## Arguments

- systems:

  data.frame with columns `gene`, `system_id`, `system`, and optionally
  `module_id` and `module`. One row per gene-within-system (a gene may
  appear in several systems).

- polarity:

  data.frame with columns `gene`, `system_id`, `w_g`. `w_g` may be
  `-1/0/+1` or continuous in `[-1, 1]`. Optional columns `evidence`,
  `role`, `confidence`, `source` are carried through.

- anchors:

  Optional data.frame with `system_id`, `gene` - the genes that define
  the direction of "more activation" for each system.

- anchor_method:

  Character, how the anchors were obtained. One of `"curated"`,
  `"graph_sink"`, `"user"`, or `"none"` (systems then carry no sign and
  only gene-level analysis is valid).

- name, version, notes:

  Provenance, printed and stored.

## Value

An object of class `dmsa_reference`.

## Examples

``` r
## One row per gene-within-system, a sign per gene, and the anchor genes that
## fix which way "more activation" points for each system.
sys <- data.frame(gene = c("CRH", "POMC", "NR3C1", "OXT", "OXTR"),
                  system_id = c("1", "1", "1", "2", "2"),
                  system = c(rep("HPA axis", 3), rep("Oxytocin", 2)))
pol <- data.frame(gene = sys$gene, system_id = sys$system_id,
                  w_g = c(1, 1, -1, 1, 0.6))
anc <- data.frame(system_id = c("1", "2"), gene = c("CRH", "OXT"))
dmsa_reference(sys, pol, anc, anchor_method = "curated", name = "toy",
               version = "0.1")
#> dmsa reference: toy (0.1)
#>   5 genes across 2 systems  (no module layer)
#>   polarity: 5 gene-system pairs, continuous (4 activating, 1 braking, 0 off-axis)
#>   anchors: curated (2 genes)
```
