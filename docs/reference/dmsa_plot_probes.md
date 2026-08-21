# Panel 3: probe level, with the CpG-to-expression direction applied

The panel that declares a direction. Each probe's slope is reflected
according to its `cpgdirection` call so that every probe is read on one
axis; reflected probes are drawn as open symbols and the convention is
printed under the title.

## Usage

``` r
dmsa_plot_probes(
  probes,
  invert = c("+1", "-1", "none"),
  gene_summary = NULL,
  title = "",
  axis_label = NULL,
  alpha = 0.05,
  file = NULL,
  width = 1700,
  height = 1400,
  res = 190
)
```

## Arguments

- probes:

  data.frame for ONE system and outcome, restricted to the genes worth
  showing: `gene`, `probe`, `b`, `se`, `d` (the cpgdirection call, +1 or
  -1); optionally `p` and `selected`.

- invert:

  Which `d` is reflected onto the common axis: `"+1"` (default) reflects
  the probes whose methylation predicts HIGHER expression, `"-1"`
  reflects the others, `"none"` shows raw slopes.

  With `"+1"` the axis reads: **positive = the exposure moves
  methylation in the direction that LOWERS this gene's expression**.
  Every d = -1 probe is shown as fitted (more methylation, less
  expression) and every d = +1 probe is reflected so that it means the
  same thing.

- gene_summary:

  Optional data.frame with `gene` and `z`, and optionally `p_adj`: the
  pooled gene-level result, drawn as a diamond on its own row. This is
  usually the point of the figure - individual probes are rarely
  significant on their own, and the evidence lives in their agreement.

- title:

  System label.

- axis_label:

  Overrides the x-axis label.

- alpha:

  Gate level. The gate is `min(1, 2*min(p_dense, p_sparse))`, so the
  boundary drawn is at `alpha/2` on each arm.

- file:

  Optional png path. If `NULL` draws to the current device.

- width, height, res:

  png dimensions.

## Value

Called for its side effect of drawing a plot on the active graphics
device. Returns `NULL` invisibly.

## Examples

``` r
## Every probe is read on one axis. A CpG whose methylation predicts HIGHER
## expression (d = +1) is reflected, so positive always means the exposure
## moved methylation the way that LOWERS this gene's expression.
pr <- data.frame(
  gene  = "NR3C1",
  probe = paste0("cg", 10001:10006),
  b     = c(-0.021, -0.014, 0.018, -0.009, 0.026, -0.017),
  se    = c(0.008, 0.007, 0.009, 0.006, 0.010, 0.008),
  d     = c(-1, -1, 1, -1, 1, -1))
dmsa_plot_probes(pr, gene_summary = data.frame(gene = "NR3C1", z = -2.9),
                 title = "HPA axis", file = tempfile(fileext = ".png"))
```
