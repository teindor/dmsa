# Re-run the frame's test drive

Re-run the frame's test drive

## Usage

``` r
dmsa_test_drive(frame)
```

## Arguments

- frame:

  a `dmsa_frame`.

## Value

the frame's corrections table, invisibly (printed first).

## Examples

``` r
set.seed(1)
map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                  probe = c("cg01", "cg02"), column = c("cg01", "cg02"),
                  best_direction = c(-1, 1), p_plus = c(.1, .9))
dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2),
                  cg01 = plogis(rnorm(40)), cg02 = plogis(rnorm(40)))
fr <- dmsa_frame(dat, map = map, outcome = "anx", covariates = "cov1",
                 B = 49, seed = 1, outdir = tempfile())

## beta values were detected and converted; every such fix is on the record
dmsa_test_drive(fr)
#> Test drive - 1 recorded action(s)
#>        field                issue
#>  methylation beta values detected
#>                                          action
#>  converted to M-values (log2 b/(1-b), eps 1e-4)
#> design condition number: 1.8
#> block permutations: ~2^61
#> pilot (B = 49, system 1): 0.0s -> full run ETA ~0 min
```
