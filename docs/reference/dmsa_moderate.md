# Moderated DMSA (mDMSA)

Moderated DMSA (mDMSA)

## Usage

``` r
dmsa_moderate(
  data,
  probes,
  alignment,
  outcome,
  moderator,
  design,
  method = c("expected", "fixed"),
  center = TRUE,
  scale_outcome = TRUE,
  B = 999L,
  engine = c("auto", "lm", "lmer"),
  beta_input = TRUE,
  check = TRUE,
  seed = NULL
)
```

## Arguments

- data:

  data.frame.

- probes, alignment:

  as in
  [`dmsa_fit()`](https://teindor.github.io/dmsa/reference/dmsa_fit.md).

- outcome:

  Character. The outcome column (a psychological measure).

- moderator:

  Character. The moderator column.

- design:

  A
  [`dmsa_design()`](https://teindor.github.io/dmsa/reference/dmsa_design.md).
  Its `focal` is ignored and replaced by `score * moderator`; its
  `fixed`, `random`, `exchangeable` and `forbid` are used as declared.

- method:

  `"expected"` or `"fixed"` scoring.

- center:

  Center the moderator (default TRUE). Leaving a moderator uncentered
  makes the product term collinear with its parents and the simple
  slopes uninterpretable; the returned object records what was done.

- scale_outcome:

  Standardise the outcome (default TRUE).

- B:

  Permutations (Freedman-Lane: residuals of the no-interaction model are
  permuted within exchangeable blocks, preserving the dependence).

- engine:

  as in
  [`dmsa_fit()`](https://teindor.github.io/dmsa/reference/dmsa_fit.md).

- beta_input:

  TRUE if `probes` hold beta values.

- check:

  Run the design check (default TRUE).

- seed:

  Optional integer.

## Value

An object of class `dmsa_moderate`.

## Examples

``` r
set.seed(1)
n <- 100; probes <- paste0("cg", 1:6)
dat <- data.frame(sex_c = rep(c(-.5, .5), length.out = n),
                  chip = rep(1:10, each = 10),
                  cID = rep(seq_len(n / 2), each = 2), mod = rnorm(n))
dat[probes] <- plogis(matrix(rnorm(n * 6), n, 6))
al <- dmsa_align(data.frame(cpg = probes, d = rep(c(1, -1), 3),
                            p_plus = rep(c(.9, .1), 3)),
                 genes = rep("OXTR", 6), level = "gene")
dat$out <- 0.6 * dmsa_score(dat, probes, al) * dat$mod + rnorm(n)
des <- dmsa_design("placeholder", "sex_c", random = "chip", exchangeable = "cID")
dmsa_moderate(dat, probes, al, outcome = "out", moderator = "mod",
              design = des, B = 99, engine = "lm", seed = 3)
#> mDMSA: out ~ tone x mod
#>   n = 100   moderator centered (mean 0.109, SD 0.898)
#>   interaction b = 0.4547  se = 0.0939  t = 4.84  p = 4.955e-06
#>   Freedman-Lane permutation p = 0.01 (B = 99, block = cID)
#>   product-term VIF = 1.06
#>   engine lm
#>   simple slopes:
#>  at_moderator moderator_value  slope    se     t        p
#>         -1 SD          -0.789 -0.392 0.141 -2.79 0.006372
#>         +1 SD           1.007  0.424 0.113  3.77 0.000284
#>   Johnson-Neyman: significant outside [-0.44, 0.49]; 63% of the sample
#>     (moderator observed range -2.21 to 2.40)
```
