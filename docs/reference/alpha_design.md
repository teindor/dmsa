# Project Alpha build contracts

The four 2026 Project Alpha builds, transcribed from the student pack's
`covariate_sets.csv` including its never-list. Supply `focal` for the
analysis at hand; everything else is fixed by the build.

## Usage

``` r
alpha_design(build, focal, drop = character())
```

## Arguments

- build:

  1 (parents T1), 2 (T1 to T4), 3 (October 7), or 4 (children).

- focal:

  Character focal term(s), passed to
  [`dmsa_design()`](https://teindor.github.io/dmsa/reference/dmsa_design.md).

- drop:

  Optional character vector of contract covariates to drop, for a
  declared deviation (e.g. dropping `Epi_T1` when the exposure is the
  immune-cell fraction, its compositional complement). Recorded in the
  returned object so it shows up in print and cannot be silent.

## Value

A `dmsa_design`.

## Details

Build 2 and 3 are longitudinal: their contract asks for the
`time x predictor` interaction, so pass e.g.
`focal = "time:BSI_Total_c"`.

## Examples

``` r
d <- alpha_design(1, focal = "IC_T1_c")
d
#> dmsa design: Alpha build 1 (parents T1)
#>   focal      IC_T1_c
#>   fixed       sex_c, age_at_array_T1, Epi_T1, Fib_T1, ctrlSV3_T1, ctrlSV5_T1 
#>   random      (1|chip_T1) + (1|cID) 
#>   exchange    cID 
#>   never       ctrlSV1_T1, ctrlSV2_T1, ctrlSV4_T1, plate_T1, submission_T1, submission_kit_T1, chip_row_T1 
# builds 2 and 3 are longitudinal: the contract asks for the time interaction
alpha_design(2, focal = "time:BSI_Total_c")
#> dmsa design: Alpha build 2 (T1 to T4)
#>   focal      time:BSI_Total_c
#>   fixed       sex_c, Epi_array, Fib_array, ctrlSV5_array 
#>   random      (1|ID) + (1|cID) 
#>   exchange    cID 
#>   never       submission_array, plate_array, chip_array, ctrlSV1_array, ctrlSV2_array, ctrlSV3_array, ctrlSV4_array, age_at_array_array, chip_row_array 
# a deviation from the contract is recorded on the object, never silent
alpha_design(1, focal = "IC_T1_c", drop = "Epi_T1")$dropped
#> [1] "Epi_T1"
```
