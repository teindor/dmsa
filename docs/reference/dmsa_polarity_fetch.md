# Draft polarity signs for a panel from public databases

Queries public signed-relationship resources for the genes of a cascade
and returns a DRAFT polarity table for a curator to check. It does not
decide anything on its own, and it is deliberately unwilling to guess.

## Usage

``` r
dmsa_polarity_fetch(
  sets,
  anchors = NULL,
  sources = c("go", "omnipath"),
  max_genes = 200L,
  quiet = FALSE
)
```

## Arguments

- sets:

  A `dmsa_sets` or `dmsa_selection` whose genes to draft signs for, or a
  character vector of gene symbols.

- anchors:

  Named list mapping `system_id` to the anchor genes that define that
  system's activation tone. Required for anything except
  `sources = "go"`: without an anchor there is no tone for a sign to
  point at.

- sources:

  Which resources to query, in order. See details.

- max_genes:

  Stop after this many genes. Drafting is one HTTP call per gene for
  `go` and `omnipath`, so raise it deliberately.

- quiet:

  Suppress progress messages.

## Value

A data.frame of drafted rows with `w_g`, `w_g_source`, `evidence`,
`citation` and `needs_review`, plus an `attr(, "unreachable")` listing
resources that did not respond.

## Details

What each resource can and cannot do is worth knowing before you rely on
it. `go` is the only one that gives a gene-to-process sign directly,
through `positive/negative regulation of X` terms, so it needs no path
chaining - but its coverage is uneven gene to gene. `signor` carries
explicit gene-to-phenotype rows, with a cell-biological vocabulary that
fits metabolic and immune systems better than behavioural ones.
`omnipath` has the broadest coverage and the weakest signs: CollecTRI
and DoRothEA, two of its inputs, both document assigning an activating
mode BY DEFAULT when they do not know, so a `consensus_stimulation` at
`curation_effort = 1` should be read as sign-unknown. `trrust` is the
most honest transcription-factor resource because it emits `Unknown`
rather than guessing.

Nothing returned here is a polarity. It is a proposal, and the
`needs_review` column is `TRUE` for every row.

## Examples

``` r
if (FALSE) { # \dontrun{
## Needs network access: each gene is one HTTP call to the chosen resources.
sel <- dmsa_select(systems = "hpa")
draft <- dmsa_polarity_fetch(sel, anchors = list("2" = c("CRH", "POMC")))
utils::write.csv(draft, file.path(tempdir(), "polarity_draft.csv"),
                 row.names = FALSE)
} # }
```
