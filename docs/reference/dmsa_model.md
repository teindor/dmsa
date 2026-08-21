# One model, block-permutation calibrated, for any term

The upper level of DMSA 2.0. Fits `formula` once, then permutes the
reduced-model residuals within permutation blocks of equal size
(Freedman-Lane) to obtain an exact p-value for one named coefficient -
including an interaction.
[`dmsa_moderate`](https://teindor.github.io/dmsa/reference/dmsa_moderate.md)
is the special case of this for `score * moderator` with a declared
[`dmsa_design()`](https://teindor.github.io/dmsa/reference/dmsa_design.md);
use that when you want the design contract enforced and simple slopes
reported, and this when the question is an arbitrary term.

## Usage

``` r
dmsa_model(
  formula,
  data,
  term,
  block = NULL,
  B = 1999,
  seed = NULL,
  nulls = FALSE
)
```

## Arguments

- formula:

  Model formula.

- data:

  data.frame.

- term:

  Name of the coefficient to test, exactly as `model.matrix` names it
  (e.g. `"S:ACE"`).

- block:

  Block labels, one per row of `data` - the family or cluster
  identifier. Rows are permuted only within blocks of EQUAL SIZE,
  because concatenating a random order of unequal-size families is not a
  valid relabelling.

- B:

  Permutations.

- seed:

  Optional integer.

- nulls:

  Optional list of pre-computed permutation null matrices, as returned
  by an earlier call on the same frame. Supplying it reuses the shared
  permutation stream instead of drawing a new one, which is what keeps
  lens, engine and level results jointly calibrated.

## Value

An object of class `dmsa_model`.

## Examples

``` r
set.seed(1)
n <- 80
d <- data.frame(S = rnorm(n), ACE = rnorm(n),
                cID = rep(seq_len(n / 2), each = 2))   # couples
d$y <- 0.5 * d$S * d$ACE + rnorm(n)
## rows are relabelled only within equal-size families, so the couple
## dependence is preserved under the null
dmsa_model(y ~ S * ACE, d, term = "S:ACE", block = d$cID, B = 99, seed = 1)
#> DMSA one-model test of 'S:ACE'
#>   n = 80 in 40 blocks
#>   b = +0.4358 (SE 0.1200)   t = +3.63   partial r = +0.385
#>   p = 0.0100  (block permutation, B = 99)   parametric p = 0.0005075
```
