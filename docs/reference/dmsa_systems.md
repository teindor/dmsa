# The systems available in a cascade, with the short names to select them by

The systems available in a cascade, with the short names to select them
by

## Usage

``` r
dmsa_systems(x = "alpha", pattern = NULL)
```

## Arguments

- x:

  A `dmsa_sets` or `dmsa_selection`, or anything
  [`dmsa_sets()`](https://teindor.github.io/dmsa/reference/dmsa_sets.md)
  accepts. A bare string that is neither `"alpha"` nor an existing file
  path is treated as `pattern`, so `dmsa_systems("hpa")` searches the
  bundled cascade.

- pattern:

  Optional case-insensitive regular expression; matched against the
  short name, the full system name and the system id.

## Value

A data.frame with `system_id`, `system_short`, `system` and the
module/gene/CpG counts. Printed as a table.

## Examples

``` r
# \donttest{
dmsa_systems()             # all 30
#>    system_id    system_short                                        system
#> 1          1        oxytocin                  Oxytocin, vasopressin & CICR
#> 2          2             hpa          HPA axis & glucocorticoid signalling
#> 3          3 steroidogenesis           Steroidogenesis (adrenal & gonadal)
#> 4          4    sex_steroids          Sex steroids: receptors & metabolism
#> 5          5             hpg                      HPG axis, puberty & KNDy
#> 6          6          libido                      Libido & sexual function
#> 7          7       pregnancy            Pregnancy, parturition & lactation
#> 8          8        dopamine                                      Dopamine
#> 9          9       serotonin                                     Serotonin
#> 10        10   noradrenaline    Noradrenaline & sympathetic nervous system
#> 11        11       monoamine  Monoamine synthesis, transport & degradation
#> 12        12          opioid                         Opioid & melanocortin
#> 13        13          orexin                     Orexin / hypocretin & MCH
#> 14        14 endocannabinoid                               Endocannabinoid
#> 15        15      kynurenine                            Kynurenine pathway
#> 16        16       glutamate              Glutamate & synaptic scaffolding
#> 17        17            gaba          GABA & excitation-inhibition balance
#> 18        18   neurotrophins Neurotrophins & activity-dependent plasticity
#> 19        19            gpcr           GPCR signalling & second messengers
#> 20        20      epigenetic                          Epigenetic machinery
#> 21        21       imprinted                        Imprinted genes & DMRs
#> 22        22      one_carbon                         One-carbon metabolism
#> 23        23       circadian                         Circadian & melatonin
#> 24        24          immune                    Immune, inflammation & HLA
#> 25        25 social_synaptic            Social-synaptic & psychiatric risk
#> 26        26         thyroid                    Thyroid & neurodevelopment
#> 27        27        appetite              Appetite & metabolic programming
#> 28        28          ageing          Ageing, epigenetic clocks & telomere
#> 29        29       oxidative               Oxidative stress & mitochondria
#> 30        30   pharmacogenes                                 Pharmacogenes
#>    n_modules n_genes n_cpgs
#> 1          9      90   1690
#> 2          8      49    538
#> 3          6      35    375
#> 4          8      51    698
#> 5          6      55    490
#> 6          6      33    304
#> 7          6      25    222
#> 8          6      34    476
#> 9          7      31    300
#> 10         7      31    240
#> 11         5      19    138
#> 12         5      17    201
#> 13         4      14    107
#> 14         4      17    204
#> 15         5      24    338
#> 16         5      61   1533
#> 17         7      43    642
#> 18         6      48    608
#> 19         7      40    827
#> 20         8      87   1250
#> 21         6      44    703
#> 22         8      50    438
#> 23         5      29    336
#> 24        13     127   1084
#> 25         5      34   1028
#> 26         6      21    328
#> 27         7      41    795
#> 28         4      42    602
#> 29         6      30    316
#> 30         3      12     74
dmsa_systems("stress|hpa") # search
#>    system_id system_short                               system n_modules
#> 2          2          hpa HPA axis & glucocorticoid signalling         8
#> 29        29    oxidative      Oxidative stress & mitochondria         6
#>    n_genes n_cpgs
#> 2       49    538
#> 29      30    316
# }
```
