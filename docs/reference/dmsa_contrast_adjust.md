# Correct p-values across the contrasts of one categorical exposure

Correct p-values across the contrasts of one categorical exposure

## Usage

``` r
dmsa_contrast_adjust(x, p, direction = NULL)
```

## Arguments

- x:

  A `dmsa_contrasts` object.

- p:

  Named or ordered numeric vector of p-values, one per contrast, in the
  order of `x$contrasts`.

- direction:

  Optional character or numeric vector of the pooled directions, carried
  through so the table reports sign next to significance.

## Value

A data frame with raw and adjusted p-values.

## Examples

``` r
set.seed(1)
d <- data.frame(style = factor(sample(c("secure", "anxious", "avoidant"), 90,
                                      TRUE)),
                age = rnorm(90), cID = rep(1:45, each = 2))
des <- dmsa_design(focal = "style", fixed = "age", exchangeable = "cID")
cc <- dmsa_contrasts(d, des, reference = "secure")
# one p-value per directional contrast, corrected within the family
dmsa_contrast_adjust(cc, c(anxious = 0.01, avoidant = 0.04),
                     direction = c("+", "-"))
#>             contrast    p p_adj direction
#> 1  anxious vs secure 0.01  0.02         +
#> 2 avoidant vs secure 0.04  0.04         -
```
