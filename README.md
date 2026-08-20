# dmsa — Directional Methylation Set Analysis

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

**DMSA asks a directional question of a set of CpG probes: does this gene's — or
this biological system's — methylation move in the direction that raises or
lowers its expression?** It then names the gene responsible while paying the
multiplicity toll, and it computes each unit's answer from that unit's own
probes alone.

That last property is not incidental. Most set-based methylation engines compute
a unit's p-value partly from the other probes the analyst happened to load, so
the same finding moves when the analysis set changes — a failure mode we name
and measure as **p-dancing** (`dmsa_pdance()`). DMSA's unit statistics are
invariant to it.

---

## Why alignment matters

A CpG in a promoter and a CpG in a gene body carry opposite expression
consequences for the same change in methylation. Averaging them discards the
signal. At a 50/50 promoter/body composition the methylation-scale mean cancels
exactly, and every direction-blind engine is null **by construction** rather
than by weakness.

DMSA aligns each probe to its predicted expression consequence before
aggregating:

```
m_j  =  d_j  ×  w_g  ×  (2 p₊ − 1)
```

where `d_j ∈ {−1, +1}` is the probe's methylation-to-expression direction,
`w_g` the gene's curated polarity with respect to its system's activation tone,
and `p₊` the confidence in the direction call.

In the locked benchmark, DMSA holds power .83–.86 flat across the empirically
reported promoter/body range while direction-blind competitors track the ratio
and bottom out at the nominal rate; DMSA reports the correct expression-scale
direction in 100% of its detections at every ratio.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("teindor/dmsa")
```

Or from a release tarball:

```r
install.packages("dmsa_1.7.2.tar.gz", repos = NULL, type = "source")
```

Direction calls come from the companion package
[`cpgdirection`](https://github.com/teindor/cpgdirection) (Suggests). Everything
else is base R plus `stats`; `gt`, `lme4`, `lmerTest`, `MASS` and `nnet` are
optional and used only where noted.

---

## The two-call interface

Declare the design once, then run it.

```r
library(dmsa)

frame <- dmsa_frame(
  data      = dat,                       # participants × (outcomes, covariates)
  M         = meth,                      # participants × probes
  systems   = c("hpa", "oxytocin"),      # short names; see dmsa_systems()
  module    = TRUE,                      # system > module > gene > probe
  outcome   = "attachment_anxiety",
  covariates = "contract",               # sex, age, cell fractions, ctrlSVs, chip
  B         = 1999
)

frame                                    # prints the design + every autofix applied
res <- dmsa_report(frame)                # figures, tables, summary.md under outdir
```

### The selection set

`systems =` takes short names from the bundled selection set — the curated
`system > module > gene > probe` cascade that declares which probes form which
unit. **Everything below the chosen system defaults to full**: `"hpa"` takes all
8 modules, 49 genes and 538 CpGs the set assigns to the HPA axis.

```r
dmsa_systems()                                  # 30 systems and their short names
dmsa_select(systems = c("hpa", "oxytocin"))     # inspect before running
dmsa_select(systems = "hpa", modules = "2.6")   # narrow a level if you must
```

Short names, ids, full names and unique prefixes all resolve; ambiguity errors
rather than guesses. Bringing your own panel is a CSV in the same shape —
`dmsa_sets_template()` writes the schema, `dmsa_sets_check()` validates it,
`dmsa_sets()` loads it. Full instructions in
[`SELECTION_SETS.md`](SELECTION_SETS.md).

Module labels in the bundled set are literature-audited, and the audit travels
with the result: `dmsa_evidence()` gives the evidence tier, the membership flags
and the citation keys behind each module, and the counts print with the frame and
in `summary.md`. A module-level finding is only as good as the module's
definition, so the definition's evidence is reported next to it.

```r
dmsa_evidence(frame, which = "flagged")
```

### The test drive

`dmsa_frame()` runs a full test drive before anything else: numeric coercion
with a per-column loss report, a QR rank audit that drops aliased columns **by
name**, a block audit with equal-size permutation strata, a hard error when
fewer than ~2¹⁰ distinct permutations exist, and a B = 49 pilot with a printed
ETA. Every action lands in `frame$corrections` and is repeated in the report.

### Three lenses, one permutation stream

Each unit is tested three ways on a single shared permutation stream:

| lens | what it rewards |
|---|---|
| **coherence** | cross-probe sign agreement — a gene whose probes all point one way |
| **composite** | a block shift detected by the averaging statistic |
| **diffuse** | a confidence-weighted quadratic, thin signal spread across many probes |

An ACAT omnibus combines them. The three are three views of one problem, and a
hit along one path is the architecture, not a loophole.

### Level-local families

Correction happens **inside** level-local families — systems, then the genes of
a named system, then the probes of a surviving gene — by maxT with step-down.
This is what lets DMSA name a member without the toll that a flat correction
across every gene in the panel would impose.

### Weighting engines

```r
dmsa_triangulate(..., weighting = "combined")   # default
dmsa_triangulate(..., weighting = "flat")       # exact reproducibility reference
dmsa_triangulate(..., weighting = "reliability")
```

`reliability` weights each aligned probe by its item-rest correlation with the
rest of its unit, computed from **methylation alone** so the permutation null
stays exact and weighting changes power rather than type-I error. `combined`
fuses flat and reliability on one shared permutation stream, standardising each
statistic by its own permutation null before combination. `flat` is the
original equal-weight engine and reproduces the published anchors exactly.

**Which to use.** `combined` is the default because the shape of a real signal
is not known in advance: on real cohort data neither flat nor reliability
dominated — reliability gained on concentrated single-gene signals and lost on
diffuse system-level ones — so hedging keeps every finding. That insurance is
not free. In simulation, where the signal shape *is* known, flat weighting needs
fewer participants than combined in every cell tested: 326 against 337 for a
twelve-gene family, 472 against 764 for a three-system family. If you have
reason to believe your signal is coherent, `weighting = "flat"` will find it
with a smaller sample. If you don't, pay the premium.

---

## Other entry points

| function | purpose |
|---|---|
| `dmsa_change()` | longitudinal two-wave mDMSA: time × S × E ≡ ΔS × E |
| `dmsa_pdance()` | set-selection sensitivity — engine-agnostic p-dancing diagnostic |
| `dmsa_scores()` | aligned unit scores without inference |
| `dmsa_test()` | a single unit, one call |
| `dmsa_tree()` | the level-local cascade on its own |
| `dmsa_reference_csv()` | build or check a direction/polarity reference |
| `dmsa_sets()` / `dmsa_select()` | the selection cascade and a selection from it |
| `dmsa_systems()` | the systems available, with their short names |
| `dmsa_evidence()` | module-level evidence tier, flags and citations |
| `dmsa_sets_template()` / `dmsa_sets_check()` | build and validate your own set |
| `dmsa_polarity()` | the gene-to-system signs, with a source on every row |
| `dmsa_polarity_review()` | the rows a human still has to decide, and only those |
| `dmsa_polarity_check()` | the sign invariants: anchors, one-sidedness, sub-process drift |
| `dmsa_polarity_fetch()` | draft signs for your own panel from GO, OmniPath, TRRUST, SIGNOR |

Run `help(package = "dmsa")` for the full index (59 exported functions).

---

## Reproducing the paper

See [`REPRODUCE.md`](REPRODUCE.md). In brief:

- `inst/benchmark/` — the locked simulation contest (panels A–C, cross-check),
  one harness, fairness contract enforced.
- `inst/scripts/` — the sample-size sweeps: `dmsa_nsweep_v5_mac.R` (DMSA
  engines) and `dmsa_nsweep_v6_competitors_mac.R` (competitors, calibrated).

Real-participant data are **not** in this repository. Project Alpha methylation
and phenotype data are held under IRB approval (Reichman University 5-2020,
Israel Ministry of Health Helsinki Committee) and are available under the terms
described in the manuscript's data availability statement.

---

## Citation

See [`CITATION.cff`](CITATION.cff), or:

> Ein-Dor, T. (2026). *dmsa: Directional Methylation Set Analysis* (version
> 1.7.2). R package.

---

## License

MIT © 2026 Tsachi Ein-Dor. See [`LICENSE`](LICENSE).
