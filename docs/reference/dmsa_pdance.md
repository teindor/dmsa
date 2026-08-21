# The p-dance test: does a finding's p-value survive a change of analysis set?

Perturbs the analysis set around a fixed finding and records what
happens to each engine's p-value. The focal unit's own probes are never
touched, so any movement is attributable to the analysis set rather than
to the evidence.

## Usage

``` r
dmsa_pdance(
  engines,
  set,
  focal,
  pool = NULL,
  form = c("dropout", "addition"),
  grid = NULL,
  R = 100L,
  alpha = 0.05,
  seed = NULL
)
```

## Arguments

- engines:

  Named list of functions. Each takes one argument - the members of a
  perturbed analysis set, in whatever form `set` is given - and returns
  a single p-value. Anything expressible this way can be tested, which
  is the point: DMSA and its competitors are put through an identical
  procedure.

- set:

  The full analysis set: a character vector of probe identifiers or an
  integer vector of column indices.

- focal:

  The focal unit's own members, a subset of `set`. Never perturbed.
  These carry the evidence the finding is about.

- pool:

  For `form = "addition"`, the reservoir of probes to add, disjoint from
  `set`. Ignored for dropout.

- form:

  `"dropout"` (default) or `"addition"`.

- grid:

  Perturbation sizes. For dropout, fractions of the non-focal probes to
  remove; defaults to `c(.01, .02, .03, .05, .10, .20, .30, .50)`, which
  spans the observed range of cross-cohort probe-set divergence and then
  some. For addition, counts of probes to add; defaults to
  `c(1, 2, 5, 10, 25, 50)`.

- R:

  Replicates per grid point. Default 100. The reported quantity is a
  spread, so this is the resolution of the answer.

- alpha:

  Significance threshold used only to report whether a replicate's
  verdict differs from the unperturbed one. Default 0.05.

- seed:

  Optional integer for reproducibility. The caller's RNG state is
  restored on exit.

## Value

A `data.frame` with one row per engine, grid point and replicate:
`engine`, `form`, `size` (the fraction dropped or the count added),
`rep`, `p`, and `flip` - whether that replicate's verdict at `alpha`
differs from the unperturbed verdict. The unperturbed run is included as
`size = 0`, `rep = 0`. The baseline p-value per engine is attached as
the `"baseline"` attribute.

## Details

Two forms, both reported in the same way:

- `"dropout"`:

  Removes a random fraction of the non-focal probes. The omitted probes
  carry no focal evidence, so a finding that depends on them was never
  about its own unit.

- `"addition"`:

  Adds probes drawn from `pool`, which sits outside the analysis set.
  This is the direction a competitive engine is most sensitive to,
  because added probes enter its background.

The two forms move competitive and whole-set engines in opposite
directions, which is why running only one of them understates the
problem.

## Examples

``` r
set.seed(1)
n <- 120L; K <- 40L
f <- rnorm(n)
M <- sapply(seq_len(K), function(j) if (j <= 6) 0.7 * f + rnorm(n) else rnorm(n))
colnames(M) <- sprintf("cg%04d", seq_len(K))
y <- 0.6 * f + rnorm(n)
focal <- colnames(M)[1:6]

# Two engines on identical evidence. The competitive one scores the focal
# unit against whatever else happens to be loaded; the unit-only one uses
# the focal probes and nothing else.
competitive <- function(u) {
  rest <- setdiff(u, focal)
  if (length(rest) < 2) return(NA_real_)
  t.test(abs(cor(M[, focal], y)), abs(cor(M[, rest], y)))$p.value
}
unit_only <- function(u) {
  summary(lm(rowMeans(M[, focal]) ~ y))$coefficients["y", 4]
}

pd <- dmsa_pdance(list(competitive = competitive, unit_only = unit_only),
                  set = colnames(M), focal = focal,
                  grid = c(0.1, 0.3, 0.5), R = 20, seed = 1)

# Both start from the same evidence, and it is never touched.
round(attr(pd, "baseline"), 4)
#> competitive   unit_only 
#>      0.0419      0.0008 

# How far each engine's p travels when only the OTHER probes change.
round(tapply(pd$p, pd$engine, function(x) diff(range(x, na.rm = TRUE))), 4)
#> competitive   unit_only 
#>       0.053       0.000 

# And how often the verdict at .05 reverses. The unit-only engine cannot
# move, because nothing it reads has changed.
tapply(pd$flip, pd$engine, sum)
#> competitive   unit_only 
#>          11           0 
```
