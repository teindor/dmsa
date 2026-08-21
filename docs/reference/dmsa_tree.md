# Build a cascade tree for a set of probes from a reference bundle

Returns the `tree` argument
[`dmsa_cascade()`](https://teindor.github.io/dmsa/reference/dmsa_cascade.md)
expects: one row per probe, columns ordered outermost first.

## Usage

``` r
dmsa_tree(reference, genes, system_id = NULL, levels = NULL)
```

## Arguments

- reference:

  A `dmsa_reference`.

- genes:

  Character vector, the set-membership gene of each probe.

- system_id:

  Optional: restrict to one system (a gene can sit in several, which
  would otherwise duplicate probes).

- levels:

  Which levels to include, in order. Defaults to `c("system", "module")`
  when a module layer exists, else `"system"`. `"gene"` is always
  appended.

## Value

data.frame with one row per probe. Probes whose gene is not in the
reference get `NA` and should be routed to the flat unannotated arm.

## Examples

``` r
systems <- data.frame(
  gene = c("CRH", "NR3C1", "FKBP5", "OXT"),
  system_id = c("HPA", "HPA", "HPA", "OXT"),
  system = c("HPA axis", "HPA axis", "HPA axis", "Oxytocin"),
  module_id = c("HPA.a", "HPA.b", "HPA.b", "OXT.a"),
  module = c("Drive", "Receptors", "Receptors", "Ligand"))
ref <- dmsa_reference(systems, anchor_method = "none")

## a probe whose gene is not in the bundle gets NA and goes to the flat arm
dmsa_tree(ref, genes = c("CRH", "NR3C1", "FKBP5", "NOTINSET"),
          system_id = "HPA")
#>   system module     gene
#> 1    HPA  HPA.a      CRH
#> 2    HPA  HPA.b    NR3C1
#> 3    HPA  HPA.b    FKBP5
#> 4   <NA>   <NA> NOTINSET
```
