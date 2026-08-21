# Reproducing the manuscript

Everything below runs on simulated data and needs no access to
participant records. Real-data analyses are described at the end.

## 0. Environment

``` r

install.packages(c("data.table", "limma"))     # limma from Bioconductor
remotes::install_github("teindor/cpgdirection") # direction calls (optional)
install.packages("dmsa_1.7.2.tar.gz", repos = NULL, type = "source")
```

R ≥ 4.0. The benchmark was run under R 4.3.3 with `data.table` 1.14.10
and `limma` 3.58.1.

## 1. The locked benchmark (`inst/benchmark/`)

One generator, one fairness contract, every engine on identical data
with native inference. Statistics without native inference receive
permutation-calibrated critical values from the same null draws (1,000
null datasets per configuration; roast uses 200 rotations).

``` bash
Rscript inst/benchmark/panelA_ratio.R   300 1     # ratio sweep, chunk 1
Rscript inst/benchmark/panelA_ratio.R   300 2     # chunk 2
Rscript inst/benchmark/panelB_cascade.R 300       # level-local cascade
Rscript inst/benchmark/panelC_shape.R   300       # functional form
Rscript inst/benchmark/save_crosscheck.R          # 48 archived datasets
```

| panel | question | manuscript output |
|----|----|----|
| A | can a method read the biological direction, and how does power move with the promoter/body ratio? | Table 6 |
| B | does naming the responsible member cost error control? | Supplementary Table S6 |
| C | what does a monotone test fail to see? | Supplementary Table S5 |
| cross-check | 48 datasets + null reference, so any cell is re-derivable without re-simulating | Methods |

**The DMSA entry is called last, after every other engine has drawn from
the random stream.** The competitors therefore see exactly the numbers
they saw before DMSA existed in the harness, and their values must
reproduce exactly on re-run. If a competitor’s number moves, the harness
changed, not the method.

### A known trap, documented so it is not repeated

`camera()` must be called with `inter.gene.cor = NA`. limma’s default of
`0.01` against a true within-gene correlation of `.5–.6` makes camera
wildly anti-conservative — measured family-wise error under the null of
**.86 (12 units) and 1.00 (42 units)**, versus **.060 / .045** when the
correlation is estimated. `inst/benchmark/engines.R` does this
correctly; an earlier draft of the sample-size script did not, and every
camera number it produced was void.

Likewise, rank-based p-values have granularity `1/(B+1)`. Correcting
them with Holm across `m` units floors the adjusted p at `m/(B+1)` — at
`B = 499`, `m = 42` that is `.084`, so rejection is *arithmetically
impossible*. Use maxT on the engine’s own permutation stream, not Holm
on its rank p-values.

## 2. Sample size (`inst/scripts/`)

``` r

source("inst/scripts/dmsa_nsweep_v5_mac.R")               # DMSA engines
source("inst/scripts/dmsa_nsweep_v6_competitors_mac.R")   # competitors
```

`v5` calibrates on the running machine and prints a projected total
before starting; at `NSIM = 100`, `B = 499` expect roughly 3 h on 16
cores. `v6` reuses `v5`’s DMSA rows rather than recomputing them and
takes minutes. Both save after every grid row and resume on re-source.

Both report each engine’s **realised** family-wise error beside its
power. Read every power number against the null rate in the same table:
DMSA’s designed any-lens rule runs at .04–.12 (higher with larger
families), Holm-corrected competitors typically below .05.

Output feeds Table 5. Note that N80 tracks the Šidák-corrected analytic
requirement within about 6% in the coherent designs — the level-local
toll costs what theory says it should. The exception is the block design
(20% discordant probes), where the composite correlation understates the
evidence and DMSA names the gene at 614 against an analytic 1,146.

## 3. p-dancing

``` r

dmsa_pdance(...)     # see ?dmsa_pdance
```

Two forms. **Addition**: hold a finding’s evidence fixed and add K
probes the analyst might also have selected; track every engine’s
p. **Dropout**: omit f = 1–50% of unrelated probes at random and ask
when each engine’s verdict survives. The diagnostic is engine-agnostic
and is the basis for the manuscript’s replication argument.

## 4. Real-data analyses

Project Alpha methylation and phenotype data are held under IRB approval
(Reichman University 5-2020; Israel Ministry of Health Helsinki
Committee) and are **not** distributed here. The analysis scripts,
pre-registrations, covariate contracts, and aggregate result tables are
deposited on OSF; see the manuscript’s data availability statement for
access terms.

The external validation uses public data — MESA purified monocytes,
GSE56046/GSE56047 — and can be reproduced without any restricted access.

## Session information for the reported runs

Benchmark and container runs: R 4.3.3, Linux x86_64. Sample-size sweeps:
macOS, 16 cores,
[`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html). Seeds
are fixed in every script;
[`set.seed()`](https://rdrr.io/r/base/Random.html) is called per
replicate from a deterministic offset, so any single cell can be re-run
in isolation.
