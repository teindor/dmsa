# Declare a DMSA analysis frame

Call 1 of the two-call interface. Declares the data, the levels (system
/ module / gene / probe), the outcomes, the covariates and blocks,
optional moderation, and the direction map - validates all of it, runs
the test drive (with `autofix = TRUE` correcting what is safely
correctable, LOUDLY), and returns a `dmsa_frame` for
[`dmsa_report`](https://teindor.github.io/dmsa/reference/dmsa_report.md).

## Usage

``` r
dmsa_frame(
  data,
  methylation = NULL,
  map = "alpha",
  system = TRUE,
  module = FALSE,
  gene = TRUE,
  probe = TRUE,
  systems = NULL,
  sets = "alpha",
  outcome,
  covariates = "contract",
  random_effects = "cID",
  chip = TRUE,
  chip_effect = c("random", "fixed", "none"),
  outcome_levels = NULL,
  moderation = FALSE,
  mod = "",
  mod2 = "",
  outcome_label = NULL,
  predictor_labels = NULL,
  mod_label = NULL,
  mod2_label = NULL,
  type = c("linear", "non-linear", "exponential"),
  outcome_type = c("gaussian", "logistic", "multinomial"),
  frame_role = c("predictor", "outcome"),
  cpg_map = c("confidence", "full"),
  B = 1999,
  alpha = 0.05,
  seed = 1,
  correction = c("maxT", "minP"),
  weighting = c("combined", "reliability", "flat"),
  w_floor = 1.5,
  palette = "viridis",
  plots = TRUE,
  tables = TRUE,
  summary = TRUE,
  gene_models = FALSE,
  progress = interactive(),
  beep = interactive(),
  plot_type = c("png", "pdf"),
  table_type = c("html", "docx", "rtf"),
  outdir = "dmsa_output",
  autofix = TRUE
)
```

## Arguments

- data:

  data.frame holding outcomes, covariates, ids and (by default) the
  Alpha methylation columns.

- methylation:

  NULL (auto-detect the bundled map's columns in `data`), a character
  vector of column names, or a numeric matrix whose colnames are
  probe/column ids.

- map:

  "alpha" for the bundled Project Alpha map, or a data.frame.

- system, module, gene, probe:

  logical: which levels to analyse.

- systems:

  Optional subset of systems to analyse. Accepts the selection cascade's
  short names - `systems = c("hpa", "oxytocin")` - and also system ids
  or full names. Matching is case-insensitive and a unique prefix of a
  short name is enough. `NULL` (default) analyses every covered system.
  Everything below the system - modules, genes, probes - is taken in
  full; use
  [`dmsa_select()`](https://teindor.github.io/dmsa/reference/dmsa_select.md)
  and pass the result as `sets` to narrow a lower level. Run
  [`dmsa_systems()`](https://teindor.github.io/dmsa/reference/dmsa_systems.md)
  to print the short names.

- sets:

  The selection cascade declaring system \> module \> gene \> probe.
  `"alpha"` (default) uses the bundled, module-audited Project Alpha
  2026c cascade; a path, `data.frame`, `dmsa_sets` or `dmsa_selection`
  may be supplied instead, and `NULL` disables it (module membership
  then falls back to the legacy bundled module map). See
  [`dmsa_sets_template()`](https://teindor.github.io/dmsa/reference/dmsa_sets_template.md)
  for the schema.

- outcome:

  character vector of outcome column(s). With several, each is tested
  with the others as covariates (mutual adjustment). When
  `frame_role = "outcome"` these are the exposure/predictor column(s).

- covariates:

  "contract" for the Alpha student contract, or a character vector. chip
  is added as a fixed factor automatically when present.

- random_effects:

  column(s) defining exchangeable blocks (permutation blocks); default
  "cID".

- chip:

  how to handle the array/chip batch factor. `TRUE` (default) enters it
  as a fixed effect, using `chip_f` if present and otherwise building it
  from a `chip`/`chip_T*` column. `FALSE` leaves it out entirely - on
  the Alpha parent build the factor has 61 levels on about 400 arrays,
  so it costs 60 degrees of freedom that the control-probe surrogate
  variables already partly absorb, and 7 of those 61 chips hold a single
  array - fitted exactly by their own dummy (leverage 1, residual 0), so
  they are counted in `n` while contributing nothing to any estimate.
  `"pool"` keeps the batch adjustment but merges every level holding
  fewer than 3 arrays into one `chip_small` level, which recovers the
  degrees of freedom and removes the leverage-1 rows. Any other
  character string names the column to use. The before/after contracts
  forbid chip and override this argument.

- chip_effect:

  how the chip factor enters the model, once `chip` has named it.
  `"random"` (default, and what the Alpha covariate contract specifies)
  takes chip out of the fixed design and fits it as a one-way random
  intercept, `(1 | chip)`, by REML quasi-demeaning: a median of about 16
  effective degrees of freedom on the parent build, and zero on the
  probes where the chip variance is estimated as nil. `"fixed"` restores
  the pre-1.18.0 behaviour and enters chip as a fixed factor, at one
  degree of freedom per chip beyond the first. `"none"` drops the batch
  adjustment altogether, whatever `chip` says. Validity does not rest on
  the variance component being right: the permutation null is computed
  under the same transform, so the random-intercept fit moves power, not
  type-I error.

- outcome_levels:

  Optional length-2 character vector naming the two levels of a
  two-level outcome, e.g. `c("T1", "T4")` or `c("female", "male")`. Used
  only in the written report, so a sentence can say which level is which
  instead of falling back to `outcome = 0` / `outcome = 1`. `NULL`
  (default) uses that fallback. Purely cosmetic: no level name reaches
  the model, the permutation or the family-wise correction.

- moderation, mod, mod2:

  moderation switch and moderator column(s); see the layout note above.

- outcome_label, predictor_labels, mod_label, mod2_label:

  Optional display labels used in figures and in the written report
  instead of the raw column names. Either one per column, in order, or a
  named vector keyed by column name. Purely cosmetic: no label reaches
  the model, the permutation or the family-wise correction.

- type:

  "linear", "non-linear" (linear/quadratic/threshold arms,
  ACAT-combined), or "exponential" (adds an exp arm).

- outcome_type:

  "gaussian", "logistic", or "multinomial".

- frame_role:

  "predictor" (outcome ~ frame + covariates) or "outcome" (named columns
  predict the methylation frame).

- cpg_map:

  "confidence" or "full"; see the note above.

- B, alpha, seed, correction:

  inference settings (correction is applied within level-local families
  only).

- weighting:

  Character. Probe weighting engine within a unit: `"combined"`
  (default) fuses the flat and reliability statistics on one shared
  permutation stream, `"flat"` weights every usable aligned probe
  equally, `"reliability"` weights each probe by its item-rest
  correlation with the rest of its unit. Weights are computed from
  methylation alone, so the permutation null is unaffected.

- w_floor:

  Numeric. Lower bound applied to reliability weights before
  normalisation, so a single poorly behaved probe cannot be driven to
  zero influence. Ignored when `weighting = "flat"`. records it; FALSE:
  strict errors instead.

- palette, plots, tables, summary, plot_type, table_type, outdir:

  output styling for `dmsa_report`.

- gene_models:

  Draw real exon/intron structure under the probes in each locus panel.
  `TRUE` fetches each gene's model from Ensembl, which makes network
  calls, so it is off by default; you may also pass a table built once
  with
  [`dmsa_gene_model`](https://teindor.github.io/dmsa/reference/dmsa_gene_model.md)
  and reused for every gene.

- progress:

  Show a 0-100 Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html): a bar
  writes carriage returns, which make a piped log unreadable, so it is
  off under `Rscript` and in tests.

- beep:

  Sound a completion signal when
  [`dmsa_report()`](https://teindor.github.io/dmsa/reference/dmsa_report.md)
  finishes. `TRUE` uses beepr sound 8; a number picks another sound;
  `FALSE` is silent. Defaults to
  [`interactive()`](https://rdrr.io/r/base/interactive.html), so a run
  never makes noise during `R CMD check` or a testthat run. beepr is a
  Suggests: without it this is silently a no-op.

- autofix:

  TRUE: the test drive corrects what is safely correctable and

## Value

object of class `dmsa_frame`.

## Examples

``` r
set.seed(42)
map <- data.frame(gene = rep(paste0("G", 1:3), each = 3), system_id = 1L,
                  system = "Sim system")
map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
map$column <- paste0(map$probe, "_", map$gene)
map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
map$p_plus <- ifelse(map$best_direction > 0, 0.9, 0.1)
d <- data.frame(out1 = rnorm(60), cov1 = rnorm(60),
                cID = rep(1:30, each = 2))
sig <- 0.5 * outer(d$out1, map$best_direction * (map$gene == "G1"))
d[map$column] <- plogis(matrix(rnorm(60 * nrow(map)), 60) + sig)
dmsa_frame(d, map = map, outcome = "out1", covariates = "cov1",
           random_effects = "cID", B = 99, outdir = tempfile("dmsa_ex"))
#> DMSA frame - 60 rows, 9 mapped probes, 1 system(s)
#>   levels:    system > gene > probe 
#>   families:  Sim system (3 genes, 9 probes) 
#>   outcome(s): out1  
#>   model:     gaussian linear | outcome enters as a PREDICTOR of methylation  
#>               
#>   sets:     Project Alpha 2026c (module-audited)
#>   module evidence: 9 High, 0 Moderate
#>   polarity: 874 signed, 360 off-axis (curated 112, database 191, literature 228, heuristic 676, none 27)
#>     190 row(s) need a decision - dmsa_polarity_review()
#>   map:      cpg_map = confidence (maps agree) 
#>   engine:    combined (flat + reliability fused per lens; joint-null maxT) 
#>   inference: B = 99 | maxT within level-local families | blocks: cID (~2^108 permutations)
#>   pilot:     B = 49 in 0.1s -> ETA ~0 min for the full report
#>   corrections (test drive):
#>    - methylation: beta values detected -> converted to M-values (log2 b/(1-b), eps 1e-4)
```
