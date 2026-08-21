# A gene's exon, intron and UTR structure on real coordinates

Returns one row per drawable feature - coding exon slice, 5' UTR, 3'
UTR, or an exon whose coding status the annotation did not resolve - for
every transcript of a gene, with strand and genomic start/end. This is
what
[`dmsa_plot_locus`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
draws the gene from, and it is deliberately a plain data frame so it can
be built, inspected, corrected and cached like any other reference layer
in this package.

## Usage

``` r
dmsa_gene_model(
  gene,
  source = c("auto", "ensembl", "gff", "txdb", "ensdb", "table"),
  file = NULL,
  db = NULL,
  orgdb = NULL,
  table = NULL,
  species = "homo_sapiens",
  genome = "hg38",
  cache = NULL,
  host = "https://rest.ensembl.org",
  quiet = FALSE
)
```

## Arguments

- gene:

  Gene symbol, e.g. `"NR3C1"`.

- source:

  One of `"auto"`, `"ensembl"`, `"gff"`, `"txdb"`, `"ensdb"`, `"table"`.
  `"auto"` takes a cached copy if there is one, then Ensembl.

- file:

  GTF/GFF3 path when `source = "gff"`.

- db:

  TxDb or EnsDb object when `source` is `"txdb"` / `"ensdb"`.

- orgdb:

  OrgDb for symbol lookup with `"txdb"` (e.g. `org.Hs.eg.db`).

- table:

  A data frame already in this shape, when `source = "table"`.

- species:

  Ensembl species name. Default `"homo_sapiens"`.

- genome:

  Genome label recorded on every row. Default `"hg38"` - the assembly
  the bundled Alpha coordinates use.

- cache:

  Directory for fetched models, or `NULL` to disable. Default is a
  `dmsa` folder under the session temp directory.

- host:

  Ensembl REST host. `NULL` disables the network entirely and returns an
  empty model, which is how you guarantee an offline build.

- quiet:

  Suppress progress messages.

## Value

A data frame with the columns in `.GM_COLS`, class `dmsa_gene_model`.
Zero rows if the model could not be obtained - never a fabricated one.

## See also

[`dmsa_gene_model_check`](https://teindor.github.io/dmsa/reference/dmsa_gene_model_check.md),
[`dmsa_plot_locus`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)

## Examples

``` r
if (FALSE) { # \dontrun{
## Needs network: one lookup plus one region query per gene.
gm <- dmsa_gene_model("NR3C1")
subset(gm, canonical)[, c("feature", "start", "end", "exon_rank")]

## Or from an annotation you already have, with no network at all:
gm <- dmsa_gene_model("NR3C1", source = "gff", file = "Homo_sapiens.gtf.gz")
} # }
```
