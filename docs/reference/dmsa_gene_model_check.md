# Check a gene model before drawing from it

Check a gene model before drawing from it

## Usage

``` r
dmsa_gene_model_check(gm, verbose = TRUE)
```

## Arguments

- gm:

  A `dmsa_gene_model`, or any data frame in that shape.

- verbose:

  Print the report.

## Value

Invisibly, a list with `ok` and one entry per check.

## Examples

``` r
## AVP sits on the minus strand, so its 5' UTR is the HIGH coordinate
gm <- data.frame(
  gene = "AVP", gene_id = "ENSG00000101200", chr = "chr20",
  gene_start = 3082556, gene_end = 3084724, strand = "-",
  transcript = "ENST00000380293", transcript_biotype = "protein_coding",
  canonical = TRUE, feature = c("utr5", "cds", "cds", "utr3"),
  start = c(3084665, 3084555, 3082700, 3082556),
  end   = c(3084724, 3084664, 3082802, 3082699),
  exon_rank = c(1, 1, 3, 3), source = "gff", genome = "hg38")
dmsa_gene_model_check(gm)$ok
#> gene model: AVP | chr20 | 1 transcript(s) | 4 feature(s) | gff
#>   [ok] required columns                         
#>   [ok] has at least one feature                 4 rows
#>   [ok] feature vocabulary is closed             
#>   [ok] start <= end on every row                0 reversed
#>   [ok] one strand per transcript                
#>   [ok] one chromosome                           chr20
#>   [ok] every feature inside the gene span       0 feature(s) outside
#>   [ok] a canonical transcript is marked         1 canonical transcript(s)
#>   [ok] no overlapping features within a transcript 0 overlap(s)
#>   -> usable
#> [1] TRUE
## a feature outside the gene span means a neighbour's exons survived
bad <- gm; bad$start[1] <- 3000000; bad$end[1] <- 3000100
dmsa_gene_model_check(bad, verbose = FALSE)$ok
#> [1] FALSE
```
