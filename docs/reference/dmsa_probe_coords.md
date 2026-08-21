# Probe genomic coordinates, without Bioconductor

Called with probe ids alone, reads them from
[`cpgdirection::cpgd_cpg_positions()`](https://rdrr.io/pkg/cpgdirection/man/cpgd_cpg_positions.html) -
the same package the `cpgdirection` calls `d` come from, so nothing new
is installed. Given `file` or `url` instead, reads any delimited
manifest with base R. No annotation packages either way.

## Usage

``` r
dmsa_probe_coords(
  probes,
  file = NULL,
  url = NULL,
  cache = NULL,
  probe_col = c("probeID", "Probe_ID", "probe", "IlmnID", "Name", "cpg", "cpg_id"),
  chr_col = c("CpG_chrm", "chrm", "chr", "CHR", "seqnames", "chromosome"),
  pos_col = c("CpG_beg", "start", "pos", "MAPINFO", "position", "CpG_end")
)
```

## Arguments

- probes:

  Character vector of probe ids to look up.

- file:

  Path to a local manifest (plain text or `.gz`). Takes precedence over
  `url`.

- url:

  URL of a plain-text manifest to download once. Public Infinium
  manifests in tsv form are published by the Zhou lab
  (<https://zwdzwd.github.io/InfiniumAnnotation>); any file with a probe
  column and chromosome and position columns will do.

- cache:

  Directory to keep the downloaded manifest in, so the download happens
  once. Defaults to `tools::R_user_dir("dmsa", "cache")`.

- probe_col, chr_col, pos_col:

  Column names in the manifest. Defaults cover the common spellings; if
  they are not found the function lists the columns it did find rather
  than guessing.

## Value

data.frame with `probe`, `chr`, `pos`, restricted to `probes` and in
that order. Probes not found come back with `NA`, which
[`dmsa_plot_locus()`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
handles.

## See also

[`dmsa_plot_locus`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md),
which needs none of this.

## Examples

``` r
if (FALSE) { # \dontrun{
co <- dmsa_probe_coords(my_probes)                    # via cpgdirection
co <- dmsa_probe_coords(my_probes, file = "EPIC.hg38.manifest.tsv.gz")
dmsa_write_coords(co)                                 # now a project asset
} # }
```
