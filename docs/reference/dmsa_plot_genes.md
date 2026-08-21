# Panel 2: gene level inside one system

Panel 2: gene level inside one system

## Usage

``` r
dmsa_plot_genes(
  genes,
  null_max = NULL,
  title = "",
  alpha = 0.05,
  file = NULL,
  width = 1500,
  height = 1500,
  res = 190
)
```

## Arguments

- genes:

  data.frame for ONE system and outcome: `gene`, `z`, `n_probes`,
  optionally `p_adj` (max-T adjusted) and `selected`.

- null_max:

  Optional numeric vector of the permutation null of `max|z|` over the
  same genes; its .95 and .99 quantiles are drawn as the
  multiplicity-adjusted thresholds.

- title:

  System label.

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
set.seed(1)
genes <- data.frame(gene = paste0("G", 1:8), n_probes = sample(3:9, 8, TRUE),
                    z = c(rnorm(7), 3.6))
genes$p_adj <- 2 * pnorm(-abs(genes$z))
## the permutation null of max|z| over this same family of genes
null_max <- apply(matrix(rnorm(200 * 8), 200, 8), 1, function(r) max(abs(r)))
dmsa_plot_genes(genes, null_max = null_max, title = "HPA axis",
                file = tempfile(fileext = ".png"))
```
