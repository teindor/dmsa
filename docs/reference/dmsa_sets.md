# Load a selection cascade (system \> module \> gene \> probe)

The cascade declares the analysis structure DMSA corrects within.
Loading it is separate from loading a direction map: the cascade says
which probes form which unit, the direction map says which way each
probe points.

## Usage

``` r
dmsa_sets(x = "alpha", audit = NULL, polarity = NULL, name = NULL)
```

## Arguments

- x:

  Either `"alpha"` for the bundled, literature-audited Project Alpha
  2026c cascade, a path to a CSV (optionally `.gz`) in the same shape,
  or a `data.frame` already in that shape. See
  [`dmsa_sets_template()`](https://teindor.github.io/dmsa/reference/dmsa_sets_template.md)
  for the schema and
  [`dmsa_sets_check()`](https://teindor.github.io/dmsa/reference/dmsa_sets_check.md)
  to validate a file before use.

- audit:

  Optional module-level evidence table: a path or `data.frame` with one
  row per `module_id` and any of `evidence_strength`, `audit_status`,
  `citation_keys`, `evidence_note`, `deep_search_url`. Defaults to the
  bundled audit when `x = "alpha"`, and to `NULL` otherwise. When a user
  cascade already carries these columns they are used and no separate
  table is needed.

- polarity:

  Optional gene-to-system polarity table (path, `data.frame` or
  `dmsa_polarity`) supplying `w_g`. Defaults to the bundled Alpha
  polarity when `x = "alpha"`, to the cascade's own `w_g` column when it
  has one, and otherwise to none - in which case a system-level score
  weights every gene +1 and the print method says so. Use
  [`dmsa_polarity_fetch()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_fetch.md)
  to draft one for your own panel from public databases.

- name:

  Optional label, printed and carried into reports.

## Value

An object of class `dmsa_sets`.

## See also

[`dmsa_systems()`](https://teindor.github.io/dmsa/reference/dmsa_systems.md)
for the short names,
[`dmsa_select()`](https://teindor.github.io/dmsa/reference/dmsa_select.md)
to choose systems,
[`dmsa_evidence()`](https://teindor.github.io/dmsa/reference/dmsa_evidence.md)
for the module evidence table.

## Examples

``` r
# \donttest{
cas <- dmsa_sets()          # bundled Alpha 2026c
cas
#> dmsa selection cascade: Project Alpha 2026c (module-audited) 
#>   30 systems | 188 modules | 1234 genes | 16823 CpGs | 17193 rows
#>   data-column keys available: col_parent_T1, col_long, col_child_T4, col_child_maternal_T1, col_child_paternal_T1 
#> 
#>   use these short names in systems = c(...):
#>      1  oxytocin         Oxytocin, vasopressin & CICR                   9 mod    90 gene   1690 cpg
#>      2  hpa              HPA axis & glucocorticoid signalling           8 mod    49 gene    538 cpg
#>      3  steroidogenesis  Steroidogenesis (adrenal & gonadal)            6 mod    35 gene    375 cpg
#>      4  sex_steroids     Sex steroids: receptors & metabolism           8 mod    51 gene    698 cpg
#>      5  hpg              HPG axis, puberty & KNDy                       6 mod    55 gene    490 cpg
#>      6  libido           Libido & sexual function                       6 mod    33 gene    304 cpg
#>      7  pregnancy        Pregnancy, parturition & lactation             6 mod    25 gene    222 cpg
#>      8  dopamine         Dopamine                                       6 mod    34 gene    476 cpg
#>      9  serotonin        Serotonin                                      7 mod    31 gene    300 cpg
#>     10  noradrenaline    Noradrenaline & sympathetic nervous system     7 mod    31 gene    240 cpg
#>     11  monoamine        Monoamine synthesis, transport & degradation   5 mod    19 gene    138 cpg
#>     12  opioid           Opioid & melanocortin                          5 mod    17 gene    201 cpg
#>     13  orexin           Orexin / hypocretin & MCH                      4 mod    14 gene    107 cpg
#>     14  endocannabinoid  Endocannabinoid                                4 mod    17 gene    204 cpg
#>     15  kynurenine       Kynurenine pathway                             5 mod    24 gene    338 cpg
#>     16  glutamate        Glutamate & synaptic scaffolding               5 mod    61 gene   1533 cpg
#>     17  gaba             GABA & excitation-inhibition balance           7 mod    43 gene    642 cpg
#>     18  neurotrophins    Neurotrophins & activity-dependent plasticit   6 mod    48 gene    608 cpg
#>     19  gpcr             GPCR signalling & second messengers            7 mod    40 gene    827 cpg
#>     20  epigenetic       Epigenetic machinery                           8 mod    87 gene   1250 cpg
#>     21  imprinted        Imprinted genes & DMRs                         6 mod    44 gene    703 cpg
#>     22  one_carbon       One-carbon metabolism                          8 mod    50 gene    438 cpg
#>     23  circadian        Circadian & melatonin                          5 mod    29 gene    336 cpg
#>     24  immune           Immune, inflammation & HLA                    13 mod   127 gene   1084 cpg
#>     25  social_synaptic  Social-synaptic & psychiatric risk             5 mod    34 gene   1028 cpg
#>     26  thyroid          Thyroid & neurodevelopment                     6 mod    21 gene    328 cpg
#>     27  appetite         Appetite & metabolic programming               7 mod    41 gene    795 cpg
#>     28  ageing           Ageing, epigenetic clocks & telomere           4 mod    42 gene    602 cpg
#>     29  oxidative        Oxidative stress & mitochondria                6 mod    30 gene    316 cpg
#>     30  pharmacogenes    Pharmacogenes                                  3 mod    12 gene     74 cpg
#>   module evidence: 163 High, 25 Moderate, 7 flagged (heterogeneous or measurement-defined)
#>     dmsa_evidence() lists them with their citation keys
#>   polarity: 874 signed, 360 off-axis (curated 112, database 191, literature 228, heuristic 676, none 27)
#>     190 row(s) need a decision - dmsa_polarity_review()
#> 
#>   dmsa_select(systems = c("oxytocin")) -> modules, genes and probes default to "full"
dmsa_systems()                 # the short names to use in `systems =`
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
# }
```
