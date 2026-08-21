# Split a SummarizedExperiment into the pieces DMSA tests

DMSA models probes as the multivariate response, so it wants a
`samples x probes` matrix and a `data.frame` of phenotypes with rows in
the same order. A `SummarizedExperiment` stores the assay the other way
round - `features x samples` - with the phenotypes in `colData`. This
function performs that transposition once and returns the two objects
DMSA takes, so the orientation is not something every analysis script
has to get right on its own.

## Usage

``` r
dmsa_se(se, assay = 1L, probes = NULL)
```

## Arguments

- se:

  A `SummarizedExperiment`, features (probes) in rows and samples in
  columns, as returned by resources such as `recountmethylation`.

- assay:

  Which assay to use: a name or a positive integer. Defaults to the
  first. On a methylation object this is normally the beta or M-value
  matrix; DMSA standardises probes before aggregating, so either scale
  works, but M-values are the better-behaved choice for linear models.

- probes:

  Optional character vector restricting the probes returned, in the
  order given. Probes absent from the assay are an error, since a
  silently shortened panel changes which family a unit is corrected in.

## Value

A list with two elements: `M`, a `samples x probes` numeric matrix, and
`data`, a `data.frame` of the object's `colData` with rows in the same
order as `M`. Pass them to
[`dmsa_triangulate()`](https://teindor.github.io/dmsa/reference/dmsa_triangulate.md)
or
[`dmsa_frame()`](https://teindor.github.io/dmsa/reference/dmsa_frame.md).

## Details

Probe names are taken from the assay's rownames and become the column
names of `M`, which is what
[`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
matches its direction calls against. An assay with no rownames is an
error rather than a warning: with no probe identifiers there is nothing
to align to, and a positional match would be a guess.

## Examples

``` r
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  set.seed(1)
  counts <- matrix(rnorm(20 * 30), nrow = 20,
                   dimnames = list(sprintf("cg%04d", 1:20), NULL))
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(mval = counts),
    colData = data.frame(age = rnorm(30), sex = rep(c("F", "M"), 15)))
  parts <- dmsa_se(se)
  dim(parts$M)        # samples x probes, the orientation DMSA wants
  head(parts$data, 3)
}
#>          age sex
#> 1 -0.3410670   F
#> 2  1.5024245   M
#> 3  0.5283077   F
```
