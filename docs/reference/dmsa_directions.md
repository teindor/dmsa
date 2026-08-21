# The direction calls DMSA ships with

Returns per-probe methylation-to-expression direction calls from the map
bundled with dmsa, so an alignment can be built without obtaining a
direction resource separately.

## Usage

``` r
dmsa_directions(probes = NULL, genes = NULL, tissue = "blood")
```

## Arguments

- probes:

  Character vector of CpG identifiers. `NULL` (default) returns the
  whole map.

- genes:

  Optional character vector restricting the result. Either one gene per
  entry of `probes`, selecting those specific pairs, or a shorter vector
  naming genes to keep.

- tissue:

  Which bundled layer. Only `"blood"` ships at present.

## Value

A `data.frame` with `probe`, `gene`, `d` (`-1` or `+1`), `p_plus` (the
probability the direction is `+1`) and `tier` (evidence grade). The
map's provenance and version are attached as the `"dmsa_direction_map"`
attribute.

## Details

The bundled map is the blood layer of `cpgdirection`, restricted to
probe-gene pairs that carry an actual call. About seven in ten rows of
the full table abstain, and an abstention cannot be aligned to, so those
rows are not shipped. What remains is 488,204 pairs over 237,761 probes
and 5,376 genes. Other tissues, the brain bridge and the full table live
in `cpgdirection`
([doi:10.5281/zenodo.22024185](https://doi.org/10.5281/zenodo.22024185)
), which remains the complete resource.

A probe may map to more than one gene - a median of two, up to seven -
and every mapped pair is returned, because which gene a probe is read
against is a modelling decision rather than a lookup.

## Examples

``` r
# What the bundled map knows about one gene's probes.
avp <- dmsa_directions(genes = "AVP")
nrow(avp)
#> [1] 238
head(avp, 3)
#>        probe gene  d p_plus tier
#> 1 cg00052046  AVP -1 0.3130    A
#> 2 cg00176879  AVP -1 0.1472    A
#> 3 cg00308631  AVP -1 0.3286    A

# Coverage of a probe set you care about.
probes <- c("cg00000029", "cg00000165", "not_a_real_probe")
got <- dmsa_directions(probes)
got
#>        probe      gene d p_plus tier
#> 1 cg00000029 LOC643802 1 0.8368    A
#> 2 cg00000029      SCP2 1 0.7805    A
#> 3 cg00000165 LOC149351 1 0.4370    B

# Provenance travels with the map.
attr(dmsa_directions(), "dmsa_direction_map")[c("tissue", "version")]
#> $tissue
#> [1] "blood"
#> 
#> $version
#> [1] "2026.08"
#> 
```
