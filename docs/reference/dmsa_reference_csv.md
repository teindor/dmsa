# Read a user-supplied system / module / gene csv as a reference bundle

One flat file, one row per gene-within-module (or per gene-within-system
if there is no module layer). Column names are matched case- and
punctuation-insensitively, so `System`, `system_id` and `"System ID"`
are all accepted.

## Usage

``` r
dmsa_reference_csv(
  path,
  name = NULL,
  version = NA_character_,
  notes = NULL,
  min_genes = 3L,
  quiet = FALSE
)
```

## Arguments

- path:

  Path to the csv (anything
  [`utils::read.csv`](https://rdrr.io/r/utils/read.table.html) accepts).

- name, version, notes:

  Provenance for the bundle; `name` defaults to the file name.

- min_genes:

  Systems with fewer genes than this are dropped with a message - a set
  of two genes is not a set test. Set to 1 to keep everything.

- quiet:

  Suppress the summary message.

## Value

A `dmsa_reference`.

## Details

**Required:** a system column (`system`, `system_id`, `pathway`, `set`)
and a gene column (`gene`, `symbol`, `gene_symbol`).

**Optional:**

- module (`module`, `subsystem`, `submodule`, `module_id`) - adds the
  module level to the cascade.

- `w_g` (`wg`, `weight`, `polarity`, `sign`) - the gene's effect on
  system activation, `-1`/`0`/`+1` or continuous in `[-1, 1]`. Without
  it the bundle carries no polarity and only gene-level DMSA is valid.

- anchor (`anchor`, `is_anchor`, `driver`) - TRUE for the genes that
  define what "more activation" means. Without it, anchors are taken
  from `role == "driver"` if a role column exists, else inferred from
  the most positive `w_g` with a warning.

- role, confidence, source - carried through untouched.

## Examples

``` r
f <- tempfile(fileext = ".csv")
dmsa_reference_template(f)
#> template written to /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpEgjJxj/file604839f61a80.csv
#>   system + gene are required; module, w_g, anchor, role, confidence, source are optional.
#>   w_g: +1 raises system activation, -1 opposes it, 0 off-axis; continuous values in [-1, 1] are allowed.
#>   anchor: TRUE for the genes that DEFINE activation for that system.
ref <- dmsa_reference_csv(f)
#> dmsa reference: file604839f61a80.csv
#>   8 genes across 2 systems, 6 modules
#>   polarity: 8 gene-system pairs, continuous (6 activating, 2 braking, 0 off-axis)
#>   anchors: user (2 genes)
```
