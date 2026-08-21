# The locus panel drawn by Gviz

The same probes and the same gene model as
[`dmsa_plot_locus`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md),
handed to Gviz instead of drawn here. Use this when you want what Gviz
gives and this package does not: a real cytoband ideogram, transcript
collapsing, and the ability to stack any other `GdObject` you already
build - CpG islands, ATAC peaks, ChIP tracks, conservation - in the same
coordinate frame.

## Usage

``` r
dmsa_plot_locus_gviz(
  probes,
  gene = "",
  gene_model = NULL,
  chrom = NA,
  genome = "hg38",
  extra = list(),
  ideogram = FALSE,
  file = NULL,
  width = 1800,
  height = 1100,
  res = 190
)
```

## Arguments

- probes:

  Data frame with `probe`, `pos`, `b`, and optionally `d`, `p`, `chr` -
  the same frame
  [`dmsa_plot_locus()`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
  takes.

- gene:

  Gene symbol, used for the title and the model lookup.

- gene_model:

  A `dmsa_gene_model`, or `TRUE` to fetch one.

- chrom:

  Chromosome, e.g. `"chr5"`. Taken from the model or the probes when
  absent.

- genome:

  Genome label passed to Gviz, default `"hg38"`.

- extra:

  A list of further Gviz tracks to stack under the model.

- ideogram:

  Draw the cytoband ideogram. Needs network on first use, as Gviz
  fetches the band table from UCSC.

- file:

  Optional PNG path.

- width, height, res:

  PNG geometry.

## Value

Invisibly, the list of tracks plotted.

## Details

Gviz, GenomicRanges and IRanges are Suggests, not dependencies. This
function errors with the install line if they are absent; everything
else in the package keeps working without them.

## See also

[`dmsa_plot_locus`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
for the dependency-free panel.

## Examples

``` r
if (FALSE) { # \dontrun{
gm <- dmsa_gene_model("NR3C1")
dmsa_plot_locus_gviz(pr, gene = "NR3C1", gene_model = gm)
} # }
```
