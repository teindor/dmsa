# Run a declared frame and deliver figures, tables and a summary

Call 2 of the two-call interface. Runs each requested level x outcome
with the locked statistics, applies correction ONLY within level-local
families, runs the always-on brain-bridge post-hoc on surviving probes,
and writes everything under `frame$outdir`.

## Usage

``` r
dmsa_report(frame, progress = NULL, beep = NULL)
```

## Arguments

- frame:

  a
  [`dmsa_frame`](https://teindor.github.io/dmsa/reference/dmsa_frame.md).

- progress:

  Show a 0-100 value declared on the frame, which itself defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html).

- beep:

  Completion sound. `NULL` (default) takes the frame's value. `TRUE` is
  beepr sound 8; a number picks another; `FALSE` is silent. Without
  beepr installed this is a no-op.

## Value

object of class `dmsa_report`: file paths and result tables, invisibly
printable.

## Examples

``` r
set.seed(1)
map <- data.frame(gene = rep(c("NR3C1", "FKBP5"), each = 2), system_id = 1L,
                  system = "HPA axis", probe = paste0("cg0", 1:4),
                  column = paste0("cg0", 1:4), best_direction = c(-1, 1, -1, 1),
                  p_plus = c(.1, .9, .1, .9))
dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2))
## methylation tracks anxiety in each probe's own expression direction
for (i in 1:4)
  dat[[map$column[i]]] <- plogis(rnorm(40) + .6 * dat$anx * map$best_direction[i])

## everything lands under outdir; figures and tables are switched off here
## only to keep the example quick
fr <- dmsa_frame(dat, map = map, outcome = "anx", covariates = "cov1",
                 B = 49, seed = 1, outdir = tempfile(),
                 plots = FALSE, tables = FALSE)
r <- dmsa_report(fr)
#> no curated polarity entries match the chosen systems - system-level scores weight every gene +1
#> DMSA report written to:
#>   /private/var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T/RtmpV1rBC7/file5bfa3c4b3390
#>   0 figure(s), 0 table(s), summary.md
r$results[, c("level", "unit", "n_probes", "p_omnibus")]
#>               level     unit n_probes  p_omnibus
#> anx system   system HPA axis        4 0.02000000
#> anx gene 1.2   gene    NR3C1        2 0.02000000
#> anx gene 1.1   gene    FKBP5        2 0.02400506
```
