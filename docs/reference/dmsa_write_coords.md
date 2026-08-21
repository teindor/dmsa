# Save a coordinate table so the lookup never has to be repeated

Save a coordinate table so the lookup never has to be repeated

## Usage

``` r
dmsa_write_coords(coords, path = "probe_coords.csv")
```

## Arguments

- coords:

  Result of
  [`dmsa_probe_coords()`](https://teindor.github.io/dmsa/reference/dmsa_probe_coords.md).

- path:

  csv path.

## Value

`path`, invisibly.

## Examples

``` r
## run the lookup once, then keep the csv with the project
co <- data.frame(probe = c("cg00000029", "cg00000108", "cg00000109"),
                 chr = c("16", "3", "3"),
                 pos = c(53434200, 37417716, 171916037))

f <- dmsa_write_coords(co, tempfile(fileext = ".csv"))
#> wrote /var/folders/yd/5vhpvlls31x18db51sk7y_zh0000gn/T//RtmpEgjJxj/file6048681e2325.csv - keep it with the project; the lookup is no longer needed
utils::read.csv(f)
#>        probe chr       pos
#> 1 cg00000029  16  53434200
#> 2 cg00000108   3  37417716
#> 3 cg00000109   3 171916037
```
