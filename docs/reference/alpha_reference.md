# The bundled Project Alpha tables as a reference bundle

The bundled Project Alpha tables as a reference bundle

## Usage

``` r
alpha_reference(modules = NULL)
```

## Arguments

- modules:

  Optional data.frame with `system_id`, `gene`, `module_id`, `module` to
  add a module layer (e.g. the HPA H / P / A / negative-feedback split).

## Value

A `dmsa_reference`.

## Examples

``` r
ref <- alpha_reference()
#> Warning: 10 polarity row(s) refer to gene-system pairs absent from `systems`; they will never be used
ref
#> dmsa reference: Project Alpha 2026 (draft polarity)
#>   549 genes across 30 systems  (no module layer)
#>   polarity: 116 gene-system pairs, signed +-1 (58 activating, 35 braking, 23 off-axis)
#>   anchors: curated (38 genes)
#>   note: polarity is a DRAFT pending per-gene citations
# the polarity dmsa_align() consumes for one system (2 = HPA axis)
head(dmsa_polarity_for(ref, "2"), 4)
#>    gene w_g
#> 1   CRH   1
#> 2   UCN   1
#> 3  POMC   1
#> 4 PCSK1   1
```
