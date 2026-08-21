# Locus figure: chromosome, gene with CpG lollipops, and per-CpG effects

Genomic coordinates are optional. Supply `pos` and the CpGs are placed
to scale and the chromosome strip is drawn; supply `gene_region` instead
and they are ordered along the gene (TSS1500 to 3'UTR) and spaced
evenly; supply neither and they keep the order you gave. The axis states
which of the three you got. Nothing else in the figure changes.

## Usage

``` r
dmsa_plot_locus(
  probes,
  gene = "",
  chrom = NA,
  gene_start = NA,
  gene_end = NA,
  chrom_length = NA,
  signal_p = 0.05,
  invert = c("+1", "-1", "none"),
  order_by = c("auto", "pos", "region", "given"),
  gene_model = NULL,
  transcripts = c("canonical", "all"),
  context = "",
  file = NULL,
  width = 1800,
  height = NULL,
  res = 190
)
```

## Arguments

- probes:

  data.frame for ONE gene. Required: `probe`, `b`, `se`, `d`
  (cpgdirection, +1/-1). Optional: `pos` (genomic coordinate), `chr`,
  `gene_region`, `p` (per-probe adjusted p), `signal` (logical,
  overrides `p`).

- gene:

  Gene symbol, for labelling.

- chrom:

  Chromosome, e.g. `"20"`. Taken from a `chr` column when present and
  not given here.

- gene_start, gene_end:

  Gene span in the same coordinate system as `pos`. Defaults to the
  probe range padded by 10%. Ignored when there are no coordinates.

- chrom_length:

  Length of the chromosome. Defaults to a built-in hg19 table; pass your
  annotation's value to be sure. Strip A is omitted when neither is
  available.

- signal_p:

  Probes with `p < signal_p` are drawn as carrying signal. If `probes$p`
  is absent, all probes count as signal.

- invert:

  Which `d` is reflected onto the common effect axis in strip C, exactly
  as in
  [`dmsa_plot_probes()`](https://teindor.github.io/dmsa/reference/dmsa_plot_probes.md).

- order_by:

  `"auto"` uses `pos` if present, else `gene_region`, else the supplied
  row order. Force one with `"pos"`, `"region"` or `"given"`.

- gene_model:

  Real exon, UTR and strand structure to draw under the probes. A
  [`dmsa_gene_model`](https://teindor.github.io/dmsa/reference/dmsa_gene_model.md)
  table, a path to a GTF/GFF3, `TRUE` to fetch the model from Ensembl,
  or `NULL` for none. Requires genomic positions for the probes; a model
  on a different chromosome, or one that could not be obtained, is
  dropped with a message rather than drawn.

- transcripts:

  `"canonical"` draws the annotation's own canonical transcript and
  names it under the axis; `"all"` stacks every transcript in the model,
  one lane each.

- context:

  One line of context drawn under the title - in a DMSA report this is
  the outcome the effects belong to. Without it the panel names a gene
  and a chromosome but never says what strip C's effects are effects OF,
  which is not recoverable from the figure. Empty string draws nothing.

- file:

  Optional png path.

- width, height, res:

  png dimensions. `height = NULL` scales to the number of strips and
  probes.

## Value

Invisibly, the plotted data frame, with the columns the figure actually
used (`x`, `sig`, `shown`, `flipped`).

## See also

[`dmsa_gene_model`](https://teindor.github.io/dmsa/reference/dmsa_gene_model.md)
to build the model,
[`dmsa_plot_locus_gviz`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus_gviz.md)
for the Gviz engine,
[`dmsa_probe_coords`](https://teindor.github.io/dmsa/reference/dmsa_probe_coords.md)
for probe positions.

## Examples

``` r
probes <- data.frame(probe = paste0("cg", 1:5), chr = "chr20",
                     pos = c(3082600, 3082750, 3083050, 3084600, 3084700),
                     b = c(-.04, .02, .03, -.05, -.06), se = rep(.015, 5),
                     d = c(-1, 1, 1, -1, -1), p = c(.01, .3, .04, .002, .001))
## invert = "+1" (the default) reflects the d = +1 probes, so a positive
## effect always means "moved methylation the way that LOWERS expression"
out <- dmsa_plot_locus(probes, gene = "AVP", context = "attachment anxiety",
                       file = tempfile(fileext = ".png"))
#> no gene model drawn for AVP: probes are on true coordinates but the exons are not shown. Pass gene_models = TRUE to dmsa_frame() (fetches from Ensembl), or gene_model = <table> from dmsa_gene_model() to stay offline.
out[, c("probe", "d", "shown", "flipped")]
#>   probe  d shown flipped
#> 1   cg1 -1 -0.04   FALSE
#> 2   cg2  1 -0.02    TRUE
#> 3   cg3  1 -0.03    TRUE
#> 4   cg4 -1 -0.05   FALSE
#> 5   cg5 -1 -0.06   FALSE
```
