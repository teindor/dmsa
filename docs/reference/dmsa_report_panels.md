# Build the whole report: panel 1 once, panels 2 and 3 per firing system

Build the whole report: panel 1 once, panels 2 and 3 per firing system

## Usage

``` r
dmsa_report_panels(
  systems,
  genes,
  probes,
  nulls = NULL,
  outdir = ".",
  alpha = 0.05,
  invert = "+1",
  genes_shown = c("selected", "top", "all")
)
```

## Arguments

- systems, genes, probes:

  Data frames as described in the panel functions. `genes` and `probes`
  must carry `system_id` and `outcome` so they can be split.

- nulls:

  Optional named list of permutation `max|z|` null vectors, named
  `paste(system_id, outcome, sep = "|")`.

- outdir:

  Directory for the png files.

- alpha:

  Gate level. The gate is `min(1, 2*min(p_dense, p_sparse))`, so the
  boundary drawn is at `alpha/2` on each arm.

- invert:

  Which `d` is reflected onto the common axis: `"+1"` (default) reflects
  the probes whose methylation predicts HIGHER expression, `"-1"`
  reflects the others, `"none"` shows raw slopes.

  With `"+1"` the axis reads: **positive = the exposure moves
  methylation in the direction that LOWERS this gene's expression**.
  Every d = -1 probe is shown as fitted (more methylation, less
  expression) and every d = +1 probe is reflected so that it means the
  same thing.

- genes_shown:

  Which genes go into panel 3: `"selected"` (default) or `"top"` (the
  single strongest) or `"all"`.

## Value

Invisibly, a data frame of the files written.

## Examples

``` r
## one system fires, so it gets a gene panel and a probe panel of its own
systems <- data.frame(system_id = 1:2, system = c("HPA axis", "Oxytocin"),
                      outcome = "anx", p_dense = c(0.004, 0.40),
                      p_sparse = c(0.02, 0.55), top_gene = c("NR3C1", "OXTR"))
genes <- data.frame(system_id = 1L, outcome = "anx",
                    gene = c("NR3C1", "FKBP5", "CRH"), z = c(3.4, -1.1, 0.6),
                    n_probes = c(6L, 4L, 3L), p_adj = c(0.008, 0.51, 0.88),
                    selected = c(TRUE, FALSE, FALSE))
probes <- data.frame(system_id = 1L, outcome = "anx", gene = "NR3C1",
                     probe = paste0("cg", 1:6), se = rep(.08, 6),
                     b = c(.21, .18, .09, -.14, -.20, .16),
                     d = c(-1, -1, -1, 1, 1, -1))

out <- dmsa_report_panels(systems, genes, probes, outdir = tempfile())
#> wrote 3 figure(s) to /private/var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T/RtmpEgjJxj/file604855516c2b
out[, c("panel", "system_id", "outcome")]
#>   panel system_id outcome
#> 1     1        NA    <NA>
#> 2     2         1     anx
#> 3     3         1     anx
```
