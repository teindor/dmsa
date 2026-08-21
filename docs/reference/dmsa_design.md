# Declare a DMSA design: covariates, dependence, and permutation blocks

Declare a DMSA design: covariates, dependence, and permutation blocks

## Usage

``` r
dmsa_design(
  focal,
  fixed,
  random = NULL,
  exchangeable = NULL,
  forbid = character(),
  label = NULL
)
```

## Arguments

- focal:

  Character. The focal term(s) whose pooled aligned effect is tested.
  May be a main effect (`"IC_c"`), an interaction written with a colon
  (`"time:BSI_c"`), or several terms, in which case the first is the
  tested one and the rest are mutually-adjusted co-focal terms.

- fixed:

  Character vector of fixed-effect covariate columns. These are the
  build's contract, in full - no defaults are supplied, because a
  default is exactly what goes stale.

- random:

  Grouping factors for random intercepts: a character vector
  (`c("chip","cID")`, expanded to `(1|chip) + (1|cID)`) or a one-sided
  formula (`~ (1|chip) + (1|cID)`). `NULL` declares independent
  observations - which
  [`dmsa_check_design()`](https://teindor.github.io/dmsa/reference/dmsa_check_design.md)
  will challenge if it finds repeated block IDs.

- exchangeable:

  Character. The column defining the outermost independent unit for
  permutation: couples share a family, repeated measures share a person,
  and a person sits inside a family - so the exchangeable unit is the
  family in all three cases. Rows inside a block move together and keep
  their order, which is what makes the same mechanism valid
  cross-sectionally and longitudinally.

- forbid:

  Character vector of columns that must NOT enter the model (the build's
  never-list). Supplying these turns a documentation rule into an error.

- label:

  Optional name for printing.

## Value

An object of class `dmsa_design`.

## Examples

``` r
dmsa_design(focal = "IC_c",
            fixed = c("sex_c", "age", "Fib", "ctrlSV3", "ctrlSV5"),
            random = c("chip", "cID"), exchangeable = "cID",
            forbid = c("ctrlSV4", "plate"))
#> dmsa design
#>   focal      IC_c
#>   fixed       sex_c, age, Fib, ctrlSV3, ctrlSV5 
#>   random      (1|chip) + (1|cID) 
#>   exchange    cID 
#>   never       ctrlSV4, plate 
```
