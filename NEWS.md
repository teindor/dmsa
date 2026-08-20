# dmsa 1.20.1

Test-suite fixes only. No change to any result.

* `test-frame.R` asserted the literal word "binary" in a design note. Rewording
  that note in 1.19.0 - necessary, because "switched to logistic" wrongly
  implied a logistic model is fitted - therefore failed a test whose behaviour
  had not changed. The test now asserts the claim (a note exists for the
  outcome, it identifies a two-level outcome, and it describes a contrast)
  rather than the prose.
* Three legitimate warnings in the reference tests were surfacing as uncaptured
  warnings in the suite. They are now captured with `expect_warning()`, so a
  clean run reports no warnings.

# dmsa 1.20.0

Three review items, all verified on real data.

* **A binary outcome now gets its response curve.** The shape scan used to skip
  a non-gaussian outcome entirely. The reasoning was half right: putting a 0/1
  outcome on the left of a least-squares fit is a linear probability model, and
  the Sasabuchi and Fieller statistics are undefined on it. But a binary outcome
  IS non-linear - probability is bounded at 0 and 1 - and on the LOGIT scale the
  quadratic term is perfectly well defined. Curvature is now tested by
  likelihood-ratio against the linear logistic model and the fitted response
  curve is drawn on the probability scale, with observed proportions in
  equal-count bins carrying Wilson intervals, the data as a rug, and the same
  solid-inside-90% / dashed-beyond extrapolation rule as the other shape
  figures. On the Alpha parent build this surfaced curvature at **p = .0062**
  that the old skip reported as nothing at all.

* **A duplicate unit no longer survives on an uncorrected family of one.** A
  gene with a single direction-called probe creates a PROBE family of one, where
  maxT has nothing to correct against. The identical measurement therefore came
  out at the raw p (.028) while its gene-level twin, corrected across a 6-gene
  family, failed (.108) - and because the probe was flagged a duplicate it was
  also excluded from the figures, so the run reported a survivor it never drew.
  Duplicates now inherit the adjusted p of the unit they duplicate, which is the
  multiplicity actually incurred; each unit's own value is kept in a
  `*_adj_own` column so the inheritance is auditable.

* **Legend collision in the moderation panels.** The legend was drawn at
  "topleft" over a full-height y range and sat on the scatter and the slope
  lines. Headroom is now reserved and the legend has an opaque backing.

* The overview figure title used the raw outcome column name rather than the
  supplied label, and wrapped by character count rather than measured width.

# dmsa 1.19.0

Five review comments from the Set 1 rev-3 output, all implemented.

* **Categorical outcomes were described as if they were continuous.** One
  sentence template said "Higher X was associated with ...", which produced
  **"Higher sex was associated with ..."** and **"Higher time was associated
  with ..."** in shipped output. Three cases are now told apart, and by the
  DATA rather than the column name - a two-level outcome that varies *within*
  the finest permutation block of a before/after build is a within-person
  change, one that does not is a group contrast:
  - continuous: "Higher X was associated with ..." (unchanged)
  - binary: "Relative to female, male showed methylation consistent with ..."
  - wave: "Within person, from T1 to T4, the methylation of MEG3 shifted
    toward a LOWER expression tone."
  New `outcome_levels = c("T1", "T4")` names the two levels; without it they
  are named by their coded values, which is ugly but never ambiguous.
* The design note `binary outcome under gaussian -> outcome_type switched to
  'logistic'` was misleading - methylation is the response and the outcome is a
  predictor of it, so **no logistic model is fitted**. It now says what actually
  happens.
* **Labels reached exactly one place in the summary.** With
  `outcome_label` set, the Results paragraph honoured it and the report header,
  the survivor list, the system line and the Methods paragraph all still printed
  the raw column name. All of them now use the label. Tables gain a
  `unit_label` column and keep `unit` as the join key; design notes keep raw
  names, because they refer to columns in the data file.
* **Figure text was cut off.** A module-level three-way title measured **15.25
  in of text into 6.6 in of usable width - a 2.31x overflow**, centred, so it
  was cut at both ends. Three separate causes: the p-values were concatenated
  into the title (they are now a subtitle); wrapping used a fixed CHARACTER
  count rather than measured width; and `strwidth()` measures against the
  panel's `par("cex")` while `mtext(outer = TRUE)` draws at device scale, which
  under-measured every outer title by 1.2x in a 2x2 layout. Titles are now
  fitted before the layout is set, the outer margin is sized to the number of
  lines they need, and long labels shrink rather than clip. Also fixed a
  non-terminating loop in the truncation fallback.
* **A single CpG is now named `GENE's cgXXXXXXX`** in figure titles, figure
  file names, axis labels, the summary prose and the tables. A bare probe id
  does not say which gene's level-local family it was corrected inside, and the
  family IS the gene.
* **Moderated non-linear models gain a surface figure** - 3D perspective plus a
  filled contour with the observed data drawn on it. Only for a moderated
  non-linear model: a plain shape scan has one predictor and therefore no second
  axis. Grid cells with no observation within one bandwidth are drawn grey (3D)
  or blank (contour), and the figure prints what fraction of the plotted surface
  is extrapolation - a quadratic-by-moderator surface extrapolates violently at
  the corners, and a 3D surface makes that look like data.

# dmsa 1.18.0

**Chip now enters as a random intercept, as the Project Alpha covariate
contract specifies.**

`tutor/covariate_sets.csv` states `1 parent T1, chip_T1, random, (1 | chip_T1)`.
Through 1.17.0 dmsa entered chip as a FIXED factor - 61 levels, 60 df on the
parent build, with 7 singleton chips fitted perfectly. That was an undocumented
deviation from the pre-specified contract.

Why chip is needed at all (parent build, n = 396):

* chip is **not** associated with any outcome (F ~ 1, p = .25-.73 across five
  outcomes), so it is a precision variable, not a confounder - omitting it
  cannot bias an effect estimate;
* but it is a large technical-variance component in methylation: **61%** of
  probes show chip signal beyond the control-probe SVs at p < .05 (5% expected)
  and **30%** at p < .001 (0.1% expected);
* the control-probe SVs do not substitute for it (ctrlSV3 R2 = .32,
  ctrlSV5 R2 = .36 against chip);
* families **span** 2-3 chips, so chip is crossed with the cID permutation
  block and blocking on cID does not absorb it.

Implementation: for a one-way random intercept the GLS solution is exactly OLS
on quasi-demeaned data (Balestra-Nerlove), with the variance ratio estimated by
REML. No mixed-model package is required - the REML step is a 1-D optimisation
in `R/ranef.R`. Validated against `lme4::lmer(REML = TRUE)` on the Alpha parent
build: **coefficients and t-statistics agree to 5e-13**, variance components to
4e-06, including the sigma^2_u = 0 boundary.

Because the transform is a linear map applied once, before any permutation,
Freedman-Lane block permutation is unaffected. Validity does not depend on the
variance ratio being exactly right: the null is generated under the same
transform, so the variance ratio affects power, not type-I error.

What it costs and what it buys: the random intercept spends a median of **~16
effective df** instead of 60, and **0** on the ~20% of probes where the chip
variance is estimated as nil. On the two deposit genes it is appropriately more
conservative than the fixed factor, and both still survive:

| gene | random (contract) | fixed (1.17.0) | chip omitted |
|---|---|---|---|
| PLAGL1 | adj .0260 survives | adj .0100 survives | adj .1360 |
| AVP | adj .0140 survives | adj .0080 survives | adj .0140 survives |

Omitting chip loses PLAGL1 altogether.

* New argument `chip_effect = c("random", "fixed", "none")`, default `"random"`.
  `"fixed"` reproduces <= 1.17.0 behaviour and is ~5x faster, which is useful
  for structural test runs; the Design notes now record which was used and flag
  `"fixed"` as a deviation from the contract.
* Fixed a silent bug introduced with this change: under `"random"` the chip
  column was dropped from the retained data, which would have left the grouping
  vector NULL and the random intercept never applied.

# dmsa 1.17.0

Five defects found by auditing a 29-model battery run on all four Project Alpha
builds. All five were silent: every one produced a report that compiled cleanly
and said something false or incomplete.

* **Surviving modules were never named in the report.** `hits.csv` carried them,
  the summary and the "How to report this" Results block did not. In one run
  four modules survived, the strongest at adjusted p = 0.002 - better than two
  of the three genes that *were* reported - and a reader of `summary.md` alone
  would never have known. Module survivors are now listed in the summary and
  get a Results paragraph, with the caveat that a module pools several genes.
  `"nothing survived"` is also no longer printed when a module did.

* **The significance star at system level marked the wrong quantity.** The star
  came from `selected`, which is driven by the best lens's *adjusted* p, but it
  was printed beside the ACAT *omnibus* p. That produced `omnibus 0.0658 *` and
  `omnibus 0.0820 *` - a significance mark on a p above alpha. The line now
  reads `omnibus 0.0820 (best diffuse adj 0.0220 *)`, so each star sits on the
  number it refers to.

* **The main-effect Methods block never stated the covariates.** The moderated
  block always did. A reader of a linear report could not tell what had been
  adjusted for - including in runs where a mistyped covariate was silently
  dropped, or where every covariate passed was unknown and the model was
  effectively unadjusted. The Methods text now names the covariates, the block
  structure, and warns explicitly when requested covariates were dropped.

* **"Coverage: 0 gene(s) ... untestable" in three of the four builds.** The
  untestable-gene count re-read the direction map from disk, so it still held
  parent-T1 column names and matched nothing in the child and before-after
  builds. Those reports claimed complete coverage while the parent build, on
  the same two systems and the same 18 genes, correctly reported 36 untestable.
  The raw map is now translated to the build before the intersection.

* **The console under-reported the table count by one** whenever `probes.csv`
  was written; it was never added to the tally.

# dmsa 1.16.0 (2026-08-19)

Everything in this release came out of running the 25-model Set 1 battery on all
four Project Alpha 2026d builds. Five of the twenty-five models failed and two
more passed while testing something other than what they claimed to.

## Three of the four builds could not be analysed at all (critical)

* **The direction map names probes in the parent-T1 form only.** The same CpG is
  `cg02187522_TC21_T1_AVP` in the parent build, `cg02187522_TC21_AVP` in the two
  before/after builds, and `cg02187522_TC21_T4_AVP` (plus `_maternal_T1_` and
  `_paternal_T1_` sets) in the child build. `dmsa_frame()` intersected the map's
  column names with the data's and found **zero** matches for the child, T1-T4
  and Oct-7 builds, stopping with "no mapped methylation columns found". The
  selection cascade lists every variant explicitly (`col_parent_T1`, `col_long`,
  `col_child_T4`, `col_child_maternal_T1`, `col_child_paternal_T1`), so the frame
  now translates the map to whichever form the data uses, prefers the variant
  belonging to the detected build when several match (the child workbook carries
  three), and records the translation in the corrections table. Passing
  `methylation` as an explicit vector of column names overrides the choice - that
  is how to analyse the maternal or paternal array inside the child workbook.

## Typing the contract out by hand gave a different model (serious)

* **`sex_c` was derived only under `covariates = "contract"`.** A caller who
  passed the same six columns explicitly got `sex_c: covariate absent from data
  -> dropped`, so the model was silently **not adjusted for sex** and returned
  different p-values from the identical-looking contract call. `sex_c` is a
  derived column, not a column of any workbook, and is now built whenever it is
  asked for by either route.
* **A missing covariate now raises a `warning()`**, not just a line in the
  corrections table. A mistyped covariate name is indistinguishable from an
  intended one, and dropping it quietly produces an unadjusted model that still
  reads as adjusted in the write-up.

## Moderated curvature was reported in prose and never drawn

* `moderation$selected` was keyed on the linear product's adjusted p alone. A run
  whose only survivors were **moderated curvature** (S^2 x moderator) named nine
  of them in the summary and produced **no figure at all**, and any user
  filtering `moderation$selected` saw none of them. The table now carries
  `curv_selected` and `any_selected`, and the figure loop draws both arms,
  skipping units flagged as repeating the same data at another level.
* The figure title now leads with the arm that actually carried the unit. A panel
  drawn for a curvature survivor used to be headed "family-adjusted p = 1.0000".
* `selected` is now guarded with `is.finite()`; an NA adjusted p produced NA rows.

## `chip` is an argument

* `dmsa_frame(chip = FALSE)` leaves the array batch factor out; a character
  string names the column to use; `TRUE` (default) keeps the previous
  auto-detection. On the parent build `chip_T1` has 61 levels across ~400 arrays
  with 7 singleton chips, so entering it spends 60 degrees of freedom on a
  nuisance factor that the control-probe surrogate variables already partly
  absorb. That is a judgement to be made and reported, so it had to be sayable.
  The before/after contracts forbid chip and still override the argument.

## Shared probes: audited, and the hazard is structurally absent

* 229 probes in the Alpha panel sit inside two panel genes and appear as a column
  under each name (IGF2+INS 31, GABRA5+GABRB3 25, HLA-DPA1+HLA-DPB1 23,
  PEG3+ZIM2 21, and others). Two such genes in one level-local family are not two
  independent units. `.frame_shared_probes()` now reads the cascade's explicit
  `shared_with` column rather than trying to infer sharing from duplicated probe
  names - which cannot work, because an Alpha column name embeds its gene, so one
  probe under two genes is two differently named columns.
* Auditing all 494 documented (probe, gene, partner) rows against the direction
  map: **not one shared probe survives into the map under both names** (62 under
  exactly one, 432 under neither). No DMSA analysis of this panel can
  double-count a shared probe. The check is retained as a guard for future maps
  and reports into `frame$shared_probes`.

## Two failures that gave uninformative errors

* A selected system left with no direction-called probe reached the aligner with
  a zero-row table and died inside it on `replacement has 1 row, data has 0`,
  which says nothing about what went wrong. It now stops with the system name,
  how many probes remain overall, and which systems they cover.
* When **every** mapped probe column has a missing value on the complete-case
  rows the frame now says so, and names the usual cause: the methylation set and
  the covariates were measured on different people. That is exactly the case in
  the Project Alpha child workbook, which carries the mother's and father's T1
  arrays but no maternal or paternal cell-composition or control-probe
  covariates - so those arrays cannot be analysed with the child T4 contract,
  and the frame now says that instead of failing obscurely.

## Smaller

* Moderation figure captions are wrapped to their panel (at most two lines).
  The four-line model statement on the left collided with the legend; a
  single unwrapped line then ran into the neighbouring panel's caption in the
  three-way layout and the Johnson-Neyman note ran off the panel edge. The full
  model statement belongs in the summary.
* Captions name the moderator by its display label throughout; one line still
  used the raw column name.
* The shape-scan summary said figures are drawn "for units whose adjusted
  quadratic p clears alpha". They are drawn for any surviving departure arm -
  quadratic, threshold or exponential. The sentence now says what the code does.
* The "no figures" message for a moderation run names both arms when
  `type != "linear"`.
* `dmsa_frame(chip = "pool")`, see above.
* "1 direction-called probes" now agrees in number.

# dmsa 1.15.0 (2026-08-18)

## `covariates = "contract"` is build-aware (important)

* **The contract was hard-coded to the parent T1 build.** Project Alpha ships
  four builds, each documenting a different covariate set in its own
  `Covariates` sheet. `.frame_contract()` returned the parent T1 set for all of
  them, so on the other three every contract column except `sex` was absent,
  silently dropped as "covariate absent from data", and the model ran on **sex
  alone**. On the two before/after builds that omits **`time`** - the design
  itself. The contract is now selected from the columns present:
  parent T1, child T4, or before-after. When no documented contract matches, it
  **errors** instead of degrading.

* **Repeated-measures builds get the right blocks.** The before/after builds
  list ID and cID as random effects. Permuting within cID alone would let two
  timepoints swap between members of a twin pair, breaking the within-person
  pairing; when the caller leaves `random_effects` at its default, those builds
  now block on `c("ID", "cID")` and say so.

* **chip is not entered where the contract forbids it.** The before/after builds
  list chip as NEVER include; it is now omitted there rather than added as a
  fixed factor.

* **The build and its contract are printed** in the frame's corrections table,
  so the model's covariate set is auditable from the output.

## Overview figure fixes

* **"Bold unit = survives" was wrong in every overview figure.** `axis()` takes
  a single `font`, not a vector, so `ifelse(selected, 2, 1)` applied the first
  element to every label: all units bold when the top row survived, none
  otherwise. Drawn with `mtext()` now, which vectorises.
* **A two-line wrapped title was clipped** at the top of the device.

# dmsa 1.14.0 (2026-08-18)

## The summary is what people will read, so it now says what it means

* **A family of one is marked as uncorrected.** At probe level a gene with a
  single direction-called CpG forms a family of one, where maxT has nothing to
  correct against and the "family-adjusted" p equals the raw permutation p. That
  number was being reported alongside genuinely corrected ones - and because it
  is the smallest, it looked like the strongest result in the run. The summary,
  the "How to report this" sentence and `print()` now all say
  **FAMILY OF 1, so uncorrected**.

* **The same data at three level names is no longer counted three times.** A
  module whose only direction-called gene is X, that gene, and (if it has one
  probe) that probe produce identical statistics. They were listed as separate
  survivors, multiplying the apparent number of findings: "6 survive" was really
  4 distinct results, and "2 survive" was 1. Duplicates are now detected by
  identical outcome, b and t, labelled `= <coarser unit>`, and a distinct count
  is given.

* **The Methods paragraph covers the moderated non-linear model.** It described
  only the two-way product, so a reader copying it would not know S^2 had ever
  been fitted. It now states the second model and its tested term.

* **The curvature paragraph interprets its sign** - positive means the curve
  becomes more convex as the moderator rises, negative more concave - instead of
  listing coefficients with no reading.

* **Curvature panels no longer say "Johnson-Neyman: slope".** The no-boundary
  branch kept the slope wording over a curvature axis.

* **Labels reach the rest of the report**: the moderation and curvature
  paragraphs, panel subtitles and the J-N/curvature axis were still printing raw
  column names while the titles used the labels.

* Singular "1 direction-called CpG".

# dmsa 1.13.0 (2026-08-18)

## Outcome types, labels, and the non-linear moderation figure

* **A nominal outcome is now refused instead of silently scored (bug).** The
  main path regresses methylation ON the outcome and tests a single 1-df
  coefficient. With `outcome_type = "multinomial"` a k-level nominal variable
  was coerced to numeric and analysed as an ordered, equally-spaced score -
  "delivery mode 5 is five times delivery mode 1". `dmsa_frame()` now errors and
  points at `dmsa_outcome()` (which implements a rank-based scheme for
  multinomial and ordinal outcomes) or a binary recode. The gaussian-with-few-
  levels note now says the variable is being treated as ordered.

* **`outcome_type` no longer implies a link function that does not exist.** The
  frame printed "model: logistic linear" while fitting OLS. It now states that
  the outcome enters as a PREDICTOR of methylation, so a binary outcome is a
  group contrast on a continuous response and no link function is needed or
  used. That is correct for the main path - and the shape scan, which puts the
  outcome on the LEFT, now refuses a non-gaussian outcome rather than fitting a
  linear probability model under Sasabuchi/Fieller machinery.

* **A covariate derived from an outcome is excluded up front (bug).** With
  `outcome = "sex"` the contract built `sex_c` from `sex`, then found the design
  rank-deficient, then dropped it - three messages leaving the user to work out
  that the covariate WAS the outcome. It is now excluded before fitting, with
  that stated as the reason. (The final model was correct either way.)

* **The non-linear moderation figure draws the non-linear model.** With
  `type != "linear"` the panel drew straight simple slopes while the table
  reported moderated curvature. It now fits and draws the curved model, the
  right-hand panel shows CURVATURE as a function of the moderator (J-N is about
  a slope, and with curvature the slope depends on the score too), and the title
  carries the moderated-curvature adjusted p.

* **Three-way panels name their `mod2` level.** The row label lived only in the
  outer note, so each panel was ambiguous on its own.

* **Display labels.** `dmsa_frame()` gains `outcome_label`, `predictor_labels`,
  `mod_label` and `mod2_label` - one per column, in order, or a named vector.
  They appear in figures and in the written report in place of raw column names.
  Purely cosmetic: no label reaches the model, the permutation or the correction.

# dmsa 1.12.0 (2026-08-18)

## Moderated non-linear is now a real test, and the shape prose reaches the report

* **`type` was silently ignored under `moderation = TRUE` (bug).** The moderation
  branch of `dmsa_report()` returns before the shape scan runs, so asking for a
  moderated non-linear model produced a moderated LINEAR one with no warning.
  There is now an actual moderated non-linear test: whether the CURVATURE of the
  tone score depends on the moderator, tested as the moderator by squared-score
  product with S, S^2, the moderator and S x moderator all present for
  marginality, maxT-corrected within the same level-local family. On the Alpha
  HPA/SNS panel 6 of 106 units survive it.

* **"How to report this" now covers the shape scan.** It was main-effect only, so
  a run with `type = "non-linear"` handed the user a Methods paragraph
  describing a purely linear analysis - which is what they would then publish.
  There is now a Results-shape paragraph stating the incremental test and its
  adjusted p, a Methods-shape paragraph, and a reference list (Fieller 1943;
  Ganzach 1997; Lind & Mehlum 2010; Morris et al. 2023; Sasabuchi 1980;
  Simonsohn 2018). The linear Methods paragraph now says outright that the
  unit-level tests are linear in the aligned score.

* **The shape table ranked by the quadratic arm only (bug).** Units surviving on
  the threshold or exponential arm were pushed off the bottom of the truncated
  table while still getting a figure, so the table and the figures disagreed.
  Ranking is now by the best adjusted departure p across arms, and any surviving
  unit is always shown.

* **A family of one is marked as such.** At probe level a gene with a single
  direction-called probe forms a family of one, where maxT cannot correct
  anything and the "adjusted" p equals the raw p. The table now says so on the
  row instead of presenting an uncorrected number as corrected.

# dmsa 1.11.0 (2026-08-18)

## The scan now covers the whole declared hierarchy, and the exponential arm is real

* **Probe level added, so `system > module > gene > probe` is true.** `.rp_units()`
  had no probe branch, so the shape scan and the moderated path both stopped
  silently at the gene level while the declared design said otherwise. Only the
  main-effect path ever reached probes. There is now a probe branch, with
  family = the probes of one gene, matching the family the main-effect path
  corrects within. Both the shape scan and moderation use it. On the Alpha
  HPA/SNS panel the scan goes from 22 units to 53: 2 systems, 9 modules,
  11 genes, 31 probes.

* **Every departure arm is maxT-corrected, not just the quadratic.** An
  `type = "exponential"` run reported its own arm raw while the quadratic beside
  it in the same table was corrected. Each arm now gets Westfall-Young maxT
  inside its level-local family on its own permutation stream.

* **The figure gate no longer keys on the quadratic alone.** A run whose
  exponential arm survived drew no figure, because the gate asked only whether
  the quadratic had survived. It now fires on any surviving departure arm, and
  the exponential fit is drawn and carries its p-values in the legend.

* **The panel carries numbers; the summary carries the argument.** The figure
  header had grown to seven lines of prose that collided with the legend. It is
  now three short lines - unit, shape verdict with the Fieller interval, and the
  solid/dashed convention - with every arm's raw and adjusted p in a single
  legend. The reasoning about marginality, Sasabuchi and Fieller stays in
  `summary.md` where it can be read properly.

* **A caution on the exponential arm.** exp(score) is a monotone transform that
  correlates strongly with the linear score, so the incremental test has low
  power and a null there is weak evidence. It asks whether the association
  accelerates, not whether it is non-monotone - only the shape test speaks to
  that. Stated in the summary whenever the arm is present.

# dmsa 1.10.0 (2026-08-18)

## The shape scan runs at every level, and is corrected like everything else

* **Scanned at system, module and gene, not only system.** An aligned tone score
  exists at every level of the cascade, so there was never a reason to ask about
  functional form at one level only. `type = "non-linear"` now scans each
  declared level.

* **Westfall-Young maxT inside level-local families.** The scan previously
  reported raw p-values and was corrected nowhere - the one part of DMSA sitting
  outside the error control the rest of the package is built on. That is
  untenable once the scan covers dozens of units. The quadratic departure term
  is now maxT-corrected within its level-local family on a shared permutation
  stream, exactly like every other DMSA statistic. It changes conclusions: in the
  Alpha HPA/SNS panel, ADRA2A's quadratic term is p = .032 raw and p = .060
  adjusted, so uncorrected it reads as a finding and corrected it does not.

* **Fieller interval for the extremum.** Lind & Mehlum present containment of a
  confidence interval for the turning point within the observed range as the
  test equivalent to the slope-sign test, and warn that the delta method is
  severely biased in finite samples. The scan now reports a Fieller 95% interval
  and whether it is contained. An unbounded interval - which is common, and
  which the delta method would hide behind a tidy-looking standard error - is
  reported as unbounded, because it means the turning point is not determined.

* **Figures are drawn for units that survive their family**, not for every
  scanned unit, and the run says so when nothing survives.

* **Documented against the source literature.** The Sasabuchi hypotheses, the
  max-of-one-sided-p convention, the interior-interval suggestion and the
  Fieller recommendation were each checked against Lind & Mehlum rather than
  reconstructed from memory.

# dmsa 1.9.0 (2026-08-18)

## The shape scan now respects marginality, and will not call a curve a U

Reported shape numbers change in this release. Anything taken from a 1.8.9
shape table came from the wrong specification and should be re-run.

* **Departure terms are now fitted WITH the linear term (breaking).** The scan
  substituted the transformed score *for* the linear score, so the quadratic arm
  fitted `outcome ~ score^2 + covariates` with no linear term. Because the score
  is standardised, that imposes a symmetric curve whose turning point is pinned
  at the mean, and a significant coefficient there is not evidence of curvature
  (the marginality principle; Morris et al. 2023). Each departure term is now
  entered alongside the linear term, so the reported p is the 1-df incremental
  test of curvature over and above the linear component - equivalent to
  comparing the linear and quadratic models directly. The threshold arm is fixed
  the same way.

* **The omnibus no longer pools the linear arm.** ACAT over {linear, quadratic,
  threshold} answered "is there any association", which a significant linear
  term alone can drive and which therefore cannot support a claim about shape.
  The omnibus now combines the departure terms only and is reported as such,
  with the linear p in its own column.

* **A U-shape claim now requires a U-shape test.** A significant quadratic term
  is consistent with a relationship that rises throughout but with increasing
  steepness, so it does not establish a U or an inverted U (Lind & Mehlum 2010;
  Simonsohn 2018). The report adds the Sasabuchi slope-sign test - the fitted
  slope must take opposite signs at the two ends of the observed range, with the
  p-value taken as the LARGER of the two one-sided p-values - and declares a
  shape only when that test clears alpha *and* the turning point falls inside the
  observed range. Otherwise the summary says to report a curvilinear component
  and stop.

* **A trimmed shape test, because the extremes are the sparse points.** Lind &
  Mehlum evaluate the slope at the minimum and maximum of the observed range,
  which in a skewed score are exactly where the data thin out. The scan repeats
  the test on the central 90% and flags any shape that passes on the full range
  but not on the trimmed one, so a U resting on a handful of tail cases is
  visible as such.

* **The figure draws the models its p-values came from.** 1.8.9 drew a
  hierarchical quadratic curve while labelling it with the squared-term-only
  p-value. Every fit shown is now the same hierarchical model that produced the
  p beside it, the quadratic's turning point is marked when it lies inside the
  data, and fits stay solid across the central 90% and dashed beyond.

* **Two cautions written into the summary.** The Sasabuchi test uses parametric
  standard errors rather than the block permutation, so it is anticonservative
  under clustering; and curvature in a multivariable model can be induced by an
  omitted product term (Ganzach 1997), so a curvilinear result is worth
  re-checking with the relevant interaction present.

# dmsa 1.8.9 (2026-08-18)

## Bug fix: `type = "non-linear"` computed a result and threw it away

* **The shape scan was unreachable.** With `type = "non-linear"` or
  `"exponential"`, `dmsa_report()` ran an ACAT omnibus over linear, quadratic
  and threshold arms of each system's aligned score - and then discarded it. The
  rows were assembled into a local list that nothing ever read. Consequences: a
  non-linear run returned a results table **byte-identical** to the linear run
  (verified: max absolute difference in `p_omnibus` = 0), the returned object had
  no shape component, and `summary.md` never mentioned that a scan had happened.
  A user asking for a non-linear analysis got a linear one with no warning. The
  scan now reaches `summary.md`, `tables/shape.*`, `print()` and
  `report$shape`, with the per-arm p-values that `dmsa_model_omnibus()` was
  already computing and `.report_shape()` was dropping.

* **It says what it is, and what it is not.** Each arm is fitted as the ONLY
  score term, so the scan asks which functional form of the whole-system score
  tracks the outcome - it is **not** a test of whether curvature adds anything
  to a linear fit, and it does not re-specify any unit-level test. The summary
  states that the system, module, gene and probe p-values remain the linear
  ones, so the scan cannot be mistaken for a non-linear model of the units.

* **A figure per scanned system+outcome.** The aligned score against the
  covariate-adjusted outcome, with the linear, quadratic and threshold arms
  drawn over the points and each arm's own p-value in the legend. Fits are
  **solid across the central 90% of the score and dashed beyond it**, with the
  count of participants outside and a density strip along the axis: a quadratic's
  tails are set by whatever few observations sit out there, and curvature carried
  by a handful of cases should not be able to look like the finding.

# dmsa 1.8.8 (2026-08-18)

## A surviving moderation now gets a figure, and the model gets stated

* **Simple slopes with partial residuals, plus Johnson-Neyman.** 1.8.7 said a
  moderation run draws nothing. That was the wrong fix for the right
  observation: an interaction that survives correction is precisely the result
  a reader cannot picture from a coefficient. Each surviving unit now gets
  `moderation_<unit>_<outcome>.png` - the tone-score slope at three levels of
  the moderator drawn over the partial residuals it was actually fitted on, and
  a Johnson-Neyman panel showing over what range of the moderator the slope is
  distinguishable from zero. A three-way run draws the pair twice, with the
  second moderator held at -1 SD and +1 SD.

* **The J-N panel says how many participants are in the significant region.**
  A boundary near the top of a skewed moderator can be perfectly real and still
  describe almost nobody. The subtitle gives the count and the percentage, and a
  density strip shows where the participants actually sit.

* **The simple slopes are drawn at observed quantiles, not at +-1 SD.** The
  +-1 SD convention puts a line at an impossible moderator value whenever the
  moderator is bounded and skewed - an adversity count has a floor at zero, and
  mean minus one SD is below it. The panel now uses the 10th, 50th and 90th
  percentiles, falling back down a ladder of quantiles when ties in a
  zero-inflated moderator collapse two of them onto the same value.

* **The Methods paragraph now enumerates the model instead of gesturing at it.**
  "All lower-order terms present" is not auditable. It now names every term in
  the fitted model, lists the covariates, states that with several outcomes each
  is additionally adjusted for the others, and says plainly that the cluster
  variable enters as a permutation block rather than as a random intercept - no
  mixed model is fitted, so there is no variance component to report.

# dmsa 1.8.7 (2026-08-18)

## Moderation runs get the same reporting support as main-effect runs

* **`summary.md` now writes a "How to report this" section for mDMSA too.**
  1.8.6 added the section but returned early when there was no main-effect
  table, so the branch a user is *most* likely to need help wording - a
  slope-on-a-slope, reported from one collapsed score, under one lens - was the
  one branch that produced nothing. The moderated block states the three things
  the main-effect wording never has to: that the unit is collapsed into a single
  subject-level tone score BEFORE the interaction is fitted, that only the
  composite lens is reported because coherence and diffuse do not carry through
  a product, and what the coefficient means in words. With two moderators it
  says so correctly - the coefficient is the three-way term, so it reports how
  much the first moderator's moderation depends on the second, rather than
  mislabelling it as a change in the simple slope.

* **A moderation run with `plots = TRUE` says why it drew nothing.** There is no
  per-probe effect and no three-lens panel once the unit is a single score. The
  run now names that, and names the switch that gets the figures back.

# dmsa 1.8.6 (2026-08-18)

## Reporting the result, not just computing it

* **A dedicated module-level figure.** `modules_<outcome>.png` plots every
  module in the selected systems across all three lenses on one family-adjusted
  scale. Modules are **blocked and named by system, with a rule between blocks**,
  because that boundary is the panel's actual claim: a module is corrected
  against its own system's modules and no others. Each row carries the module's
  literature evidence tier and direction-called probe count in the right margin -
  Moderate rendered in red, so a result resting on a weaker module definition
  cannot be read without noticing. The module level was the one tier of the
  cascade with tables but no picture.

* **`summary.md` gains a "How to report this" section.** DMSA asks a reader to
  hold three unfamiliar things at once - which **lens** carried the finding,
  which **family** it was corrected inside, and a direction that is a claim
  about **expression tone** rather than about methylation level - and a user
  handed `adj .0120, coherence, d = -1` has no way to know how to write that
  up. The report now writes the sentences, with this run's own numbers in them:
  a Results paragraph per surviving unit, a Methods paragraph carrying the
  actual settings, and the direction stated in words. A run with no survivor
  gets the calibrated-null wording instead, which is the harder one to phrase.

## Fixes

* **Every table is also written as .docx**, and if `gt` is missing the run says
  so once, up front, instead of leaving you to notice the absent Word file when
  you go to paste it. csv to compute on, html to read, docx to paste into a
  manuscript - a results table retyped by hand is a table with typos in it.
* **The overview figure's title no longer runs off the right edge.** Outcome
  names are the user's own column names, and "Attachment_Anxiety_General_T1"
  already overran the device. The title wraps to the panel width.
* **A locus panel drawn without a gene model says so on the console**, naming
  the argument that fixes it. The caption always said it, but captions get
  skimmed, and a panel with no exons is half the point of the panel.

# dmsa 1.8.5 (2026-08-18)

* `dmsa_report()` now says where it wrote, every run, with the path resolved:
  `outdir` stays relative by default - the working directory is where every
  other R tool puts results, and changing that would be the more confusing
  choice - but "dmsa_output" is only findable if you already know what it is
  relative to. The message appears whether or not the result is printed, so
  assigning `rep <- dmsa_report(fr)` no longer ends a 35-minute run in silence.
  `print()` resolves the path too.

# dmsa 1.8.4 (2026-08-18)

* The progress bar advances through the probe level and the figure/table
  phases instead of sitting at 71% and jumping to 100% at the end. The two
  trailing phases were counted in the total but never ticked.

# dmsa 1.8.3 (2026-08-18)

## Two conveniences for long runs

A report is minutes of permutation with nothing on screen - the parent cohort
at B = 1999 across 19 systems takes 35 minutes - which leaves no way to tell
"working" from "hung".

* **A 0-100% progress bar.** `dmsa_frame(progress = )`, or override at run time
  with `dmsa_report(frame, progress = TRUE)`. It counts triangulation calls, the
  unit of work that actually costs time, so the bar is proportional rather than
  a decorative sweep.
* **A completion sound.** `dmsa_frame(beep = )`. `TRUE` uses beepr sound 8, a
  number picks another, `FALSE` is silent.

Both default to `interactive()`. A progress bar writes carriage returns, which
turn a piped log into soup, and a package that makes noise during `R CMD check`
or a testthat run is one people learn to mute - so both are live at the console
and off under `Rscript`, check and tests. beepr is a **Suggests**: without it
the sound is silently a no-op, never load-bearing.

# dmsa 1.8.2 (2026-08-18)

## Regression fix: maxT had lost its information advantage under a single engine

Caught by running the shipped package against the pre-registered proof of
concept and comparing every number to the deposited results.

A change on 15 August 2026 standardised each lens statistic by its own
permutation null before the family-wise step. That was introduced to fix a real
problem in the **combined** engine - a rank-p fusion floors the family-wise p
near (family size)/(B+1), which silently zeroes power for large families at
modest B - but it was applied to the **single-engine** path too, where it does
not belong.

Standardising divides out exactly what maxT's advantage is made of. A unit whose
probes all point the same way generates the **widest** permutation null in its
family and is therefore charged the least; rescale every unit to unit spread and
that width buys nothing, leaving minP-like behaviour under a maxT label.

The cost was the pre-registered anchor. AVP's family-adjusted coherence moved
from **.0120** - the reported value, and the one the p-dancing figures quote at
every perturbation - to **.1135**, which is the documented *combined*-engine
value of .125. PLAGL1, CAT and the exploratory sweep all shifted with it, in
both directions, and a gene that is not a finding (ADCY6) appeared to clear .05.

Standardisation now happens only where two engines are actually fused. With one
engine, maxT sees the raw statistic again.

**Every number in the pre-registered proof of concept now reproduces exactly**
from the shipped package at `weighting = "flat"`, `cpg_map = "full"`, B = 1999:
AVP coherence .0120/.0120, composite .0070/.1285, diffuse .0080/.0505; PLAGL1
coherence .0055/.0425, composite .0025/.0080, diffuse .0025/.1720. All twelve
raw and adjusted values match the deposit.

Three tests now pin the property: the concordant unit's toll stays small under
one engine, maxT stays strictly kinder than minP to it, and the two-engine
fusion keeps its standardisation.

# dmsa 1.8.1 (2026-08-18)

## The locus panel was mixing genome builds, silently

Found by running the whole pipeline on the real parent cohort rather than on
test fixtures - which is the only way this class of defect surfaces, because
neither layer is wrong on its own.

* **`dmsa_report()` took probe coordinates from `dmsa_probe_coords()`, which
  reads the companion package's EPIC manifest in hg19, and drew them against an
  hg38 gene model.** AVP moves 19,353 bp between the two builds, so its eight
  promoter CpGs were drawn 19 kb from their own gene, in the middle of the
  neighbouring intergenic space, with nothing reported. Coordinates now come
  from the bundled cascade, looked up by probe id: it is hg38, it is the same
  table the gene models are aligned to, and it is always present. The manifest
  lookup remains as the fallback only.
* **A genome-build guard.** `dmsa_plot_locus()` now measures the distance from
  the probes to the gene model and refuses to pass it off as a map when they
  cannot belong together: if every probe is further from the gene than the gene
  is long (floor 5 kb), it messages naming the build and prints CHECK GENOME
  BUILD on the figure itself. The rule is distance, not overlap, because real
  promoter probes legitimately sit outside the gene span - all eight of AVP's
  do - and a naive containment test would cry wolf on exactly the figure this
  feature exists to draw.

554 unit tests.

# dmsa 1.8.0 (2026-08-18)

## The locus panel draws the gene, not a bar shaped like one

Strip B used to draw a gene as a featureless grey rectangle, and - when a probe
carried no coordinate - to space probes evenly along it. Both told the reader
nothing while looking like they told them something. A CpG in the first exon
and a CpG 40 kb into an intron are different claims about regulation, and which
one it is decides the probe's direction call. The panel now shows it.

* `dmsa_gene_model()` returns real exon, UTR and strand structure on genomic
  coordinates, from any of four sources in one shape: **Ensembl REST** (no
  install, no annotation package, cached to disk), a **GTF/GFF3** you already
  have, a **TxDb** plus OrgDb, or an **EnsDb**. Every row records which.
* `dmsa_plot_locus(gene_model = ...)` draws it: coding exons as full-height
  boxes, UTRs half-height, introns as a line with chevrons running 5' to 3',
  and the TSS marked. These are the conventions every genome browser uses, so
  a reader already knows how to read the strip. `transcripts = "all"` stacks
  every isoform; the default draws the annotation's own canonical transcript
  and names it under the axis.
* `dmsa_plot_locus_gviz()` hands the same probes and the same model to **Gviz**
  for what Gviz gives and this does not: cytoband ideograms, transcript
  collapsing, and any other track stacked in the same coordinate frame.
* `dmsa_frame(gene_models = TRUE)` turns the models on for every locus panel in
  a report. Off by default, because a report should not make network calls
  without being asked; pass a prebuilt table instead to stay offline.
* `dmsa_gene_model_check()` validates a model before it is drawn, including the
  failure this is really guarding against: a region query keeps the NEIGHBOUR
  gene's exons and draws them as if they were this gene's.

Gviz, GenomicRanges, IRanges, GenomicFeatures, ensembldb, AnnotationDbi and
AnnotationFilter are **Suggests**. The Ensembl and GFF sources need none of
them - nothing beyond base R - so the figure is not reachable only by users who
can install a gigabyte of annotation.

## The evenly-spaced fallback is no longer a silent default

It looks like a map and is not one: two probes drawn adjacent may be 40 kb
apart. When the panel has no choice it now prints **NOT TO SCALE** in bold under
the axis, where the caveat survives being cropped into a slide, and messages
once on the console naming how to fix it. The bundled Alpha cascade carries an
hg38 position for every probe, and `dmsa_report()` now reads coordinates
straight from it, so the default path is always to scale.

## Fix

* The minus strand. On a gene that reads right to left the 5' UTR is the piece
  with the HIGHER coordinate, and the TSS is at the high end. Getting this
  backwards puts the promoter at the wrong end of half of all genes; it is
  pinned by tests on both strands.

# dmsa 1.7.2 (2026-08-17)

## `dmsa_sets_check("alpha")` failed the package's own data

The data sanity gate reported the bundled cascade as **not usable** on two
counts, and both were defects in the check rather than in the data. A validator
that cries wolf on the shipped set teaches users to ignore it, so:

* **Uniqueness is keyed on `probe_id`, not `cpg`.** EPIC v2 ships replicate
  designs for the same site - `cg09834343_BC21` and `_BC22` are two columns
  measuring one CpG - so 238 legitimate rows were being counted as duplicates.
  115 CpGs in the bundled set carry more than one design; they are now reported
  as a note, since replicate designs are correlated measurements of one site and
  that matters for a probe-level family.
* **Module evidence is read from the companion audit table.** It is shipped
  beside the cascade rather than inside it, so validating the raw CSV found no
  `evidence_strength` at all. `"alpha"` is now validated as it is actually
  loaded, and a cascade genuinely carrying no evidence gets a **WARN** rather
  than a FAIL - absent evidence changes how a result reads, it does not make the
  set unusable.

## Fixes found by `R CMD check --run-donttest`

* `dmsa_systems("stress|hpa")` - the form used in this package's own
  documentation - errored with "cascade file not found". A bare string that is
  neither `"alpha"` nor an existing file is now read as a search pattern, which
  is what the example always meant.
* `dmsa_polarity_fetch()`'s example moved to `\dontrun{}` and now writes to
  `tempdir()`. It made live HTTP calls during the check: 35 seconds of timeouts
  on a machine without network, and a stray CSV left in the check directory.
* `POLARITY_AUDIT.md`, `SELECTION_SETS.md` and `NEWS.md` added to
  `.Rbuildignore`; they are repository documents, not package files.

494 unit tests.

# dmsa 1.7.1 (2026-08-17)

* `dmsa_polarity()`, `dmsa_polarity_review()` and `dmsa_evidence()` accept a
  `dmsa_sets`, `dmsa_selection` or `dmsa_frame` object, not only a file path.
  The print banner recommends `dmsa_polarity_review()`, so the object already in
  the user's hand is the obvious argument; passing one used to fail with
  "cannot coerce class dmsa_selection to a data.frame".
* The module-label check in `.pol_toward_system()` uses word boundaries. Without
  them "sex de**termination**" matched the brake vocabulary and every gonadal
  module was flagged.

# dmsa 1.7.0 (2026-08-17)

## Gene-to-system polarity becomes auditable data

`w_g` is the multiplier that decides which way a system score points, so a
single wrong sign does not degrade a result - it inverts it. Until now those
signs lived inside a curated CSV with no provenance. They are now a first-class
table in which every row states where its sign came from and what would change
it.

* `dmsa_polarity()` loads the table: 1,234 genes across 30 systems, each row
  carrying `w_g`, a `role`, a `confidence`, a `w_g_source` and, for 451 rows, a
  citation. Sources: 112 curated, 121 OmniPath, 46 GtoPdb, 23 TRRUST, 1 UniProt,
  228 literature, 676 rule-based, 27 unresolved. Signs: 583 at +1, 291 at -1,
  360 held deliberately off-axis at 0.
* `dmsa_polarity_review()` surfaces the 190 rows a human still has to decide -
  low-confidence signs, unresolved genes, framing choices and the one database
  contradiction left standing - and nothing else. `which = "disagreement"` and
  `which = "unresolved"` narrow it.
* `dmsa_polarity_check()` runs the invariants: an anchor may never lower its own
  system's tone, no system may be entirely one-sided, no sign may point at its
  module's sub-process rather than the system, and no row may carry a sign
  without a source. It reports the four known, real exceptions in the Alpha
  panel rather than passing them silently.
* `dmsa_polarity_fetch()` drafts signs for a user's own cascade from GO,
  OmniPath, TRRUST and SIGNOR. When a resource is unreachable it returns
  `unresolved` at `w_g = 0`, never a guess, and flags every row for review.
* `dmsa_polarity_sources()` states each resource's caveat, including the one
  that matters most: CollecTRI and DoRothEA both assign an **activating** mode
  by default when the mode is unknown, so OmniPath's `consensus_stimulation` at
  `curation_effort = 1` is a default rather than a finding.

## The rule the audit turned on

A gene is signed against **the system's declared tone**, not against its
immediate molecular partner. FKBP5 restrains GR and GR mediates HPA negative
feedback, so FKBP5 is **+1** to the HPA system - a brake on a brake. The role
vocabulary therefore carries `brake-of-brake` and `feedback-enabler` as
first-class labels, and separates roles that commit to a direction (`driver`,
`brake`, ...) from roles whose sign is inherited from what they serve
(`accessory`, `transducer`, `capacity`, `readout`).

A database may raise the evidence grade of a sign it agrees with. It may never
flip a sign on its own. All 38 database contradictions were adjudicated one at a
time against the literature; four proved to be defects at source. `w_g` for the
Project Alpha panel is unchanged from the curated values in every one of the
112 rows the PI set by hand. `POLARITY_AUDIT.md` records the method, the
rejected resources and the escalations.

* The polarity table travels with `dmsa_sets()` and narrows with
  `dmsa_select()`, so a selection reports the coverage of its own genes.
* `print()` on a set, a selection or a frame names the unsigned and flagged
  counts instead of hiding them.

# dmsa 1.6.0

## The selection cascade becomes a first-class object

`system > module > gene > probe` now ships with the package as a curated
selection set, and it is the object `dmsa_frame()` resolves `systems =` against.
The point is not convenience: because families are declared in this file and
never harvested from whatever the analyst happened to load, the family a unit is
corrected within is a property of the map rather than of the session. That is the
mechanism behind DMSA's invariance to set selection, and it is now visible.

* `dmsa_sets()` loads a selection cascade - the bundled, module-audited Project
  Alpha 2026c set (30 systems, 188 modules, 1,234 genes, 16,823 CpGs) or any CSV
  in the same shape.
* `dmsa_systems()` prints the systems with the **short name** for each, which is
  what `systems =` now takes: `systems = c("hpa", "oxytocin")`. Ids, full names
  and unique prefixes resolve too; ambiguity errors with the candidates named
  rather than picking one.
* **Everything below the chosen system defaults to full.** `"hpa"` takes all 8
  modules, 49 genes and 538 CpGs the set assigns to it. `dmsa_select()` narrows a
  level when a narrower question was declared, and says so in its print output,
  because narrowing changes the family and therefore the toll.
* `dmsa_select(columns = "auto", data = dat)` keys the cascade to the caller's
  matrix through the `col_*` columns instead of requiring renamed columns.
* `dmsa_frame()` gains `sets =` (default `"alpha"`). The module level now takes
  its membership from the cascade rather than the legacy `modules_alpha.csv`, so
  the family a unit is corrected within is the family the user selected from.
  `sets = NULL` restores the legacy path.

## Module evidence travels with the result

Every module label in the bundled set was checked against the literature. That
check is now part of the output, on the principle that a module-level finding is
only as good as the module's definition.

* `dmsa_evidence()` returns the evidence tier (High / Moderate), the audit status
  (`supported_as_named`, `renamed_*`, `heterogeneous_locked_membership`,
  `measurement_defined_module`), the citation keys and the evidence note, for a
  set, a selection or a frame. `which = "flagged"` reduces it to the modules that
  need a caveat when reported.
* The counts print with `dmsa_sets()`, `dmsa_select()` and `dmsa_frame()`, and
  `summary.md` from `dmsa_report()` gains a **Module evidence** table covering the
  modules the run actually reported on.

## Bringing your own panel

* `dmsa_sets_template()` writes the CSV skeleton plus the full instructions;
  `dmsa_sets_check()` validates a candidate file against the five structural
  rules before it is trusted (one system per module, one module per gene, unique
  short names, constant module metadata, no duplicate rows) and reports CpGs
  annotated to more than one gene rather than rejecting them.
* `SELECTION_SETS.md` documents the schema, the resolution rules and the
  separation between the selection set, the direction map and the polarity table.

## Fixes

* `module_id` is read as text throughout. Parsed as a number, `24.10` collapses
  onto `24.1` and two different immune modules become one - a silent, total
  corruption of the family structure.
* `systems =` no longer errors when a system id is character rather than numeric.
* `data.table` moved from an undeclared `::` call to `Imports`.
* Documented the previously undocumented `weighting`, `w_floor`, `gate`,
  `sparse_reach`, `family_correction` and `nulls` arguments.

# dmsa 1.5.0 (2026-08-15)

* ENGINE DEFAULT CONFIRMED: `combined`. The evidence points two ways and both
  directions are on the record. On REAL data, pre-registered application to the
  PTSD and children cohorts showed neither flat nor reliability weighting
  dominates - reliability gains on concentrated single-gene signals (FKBP5,
  SGK1) and loses on diffuse system and module signals - so hedging the engine
  choice keeps every real finding. In SIMULATION, however, flat weighting needs
  FEWER participants than combined in all fourteen sample-size cells, by 3% at
  a twelve-gene family (326 vs 337) and by 62% on a three-system family (472 vs
  764). The default is `combined` because the shape of a real signal is not
  known in advance and the premium for that insurance is roughly 10% more
  participants in a typical cell; `weighting = "flat"` is a first-class option
  for anyone who knows their signal is coherent, and remains the exact
  reproducibility reference for the published anchors.

* Documentation rendered: `man/` now ships 48 Rd files, so `?dmsa_frame` and
  the rest resolve. The package previously carried roxygen comments with no
  generated documentation.

* `inst/scripts/artifact_probe_mask.R` - checks EPIC probes carrying a finding
  against five published artifact lists (cross-reactive, cross-hybridising,
  SNP at target, INDEL at target, 1000G variant in probe body), with the
  provenance of every list documented in the header.

* Two tests in `test-design.R` now qualify an internal call as `dmsa:::`, so
  the suite runs outside the package namespace. 343 tests, 0 failures.

* FIXED: the combined engine's rank-p fusion had a DISCRETENESS FLOOR. Rank
  p-values are granular at 1/(B+1); under maxT across K units, ties at the
  extreme floored the family-wise p at about K/(B+1). At K = 42, B = 399 that
  is .105, so power was exactly 0.000 regardless of effect size (found on an
  independent K = 42 sweep, 15 Aug 2026; the defect was invisible in our own
  tests, which never used a family larger than five). Engines are now fused on
  a CONTINUOUS scale: each lens/engine statistic is standardised by its own
  permutation null (robust centre/scale) before combination, so no ties and no
  floor. Verified at K = 42, B = 399.

* NEW: `p_unit` / `p_unit_adj` - a unit-level test on ONE joint null across
  lenses (and engines). The per-lens adjusted p-values correct across UNITS
  only, so the previously used rule - "name a unit if the minimum of the three
  lens-adjusted p's is below .05" - is not family-wise controlled across
  lenses and inflates the error rate roughly threefold (measured: .06 where
  .05 was intended, and up to .13 in some grids). `p_unit_adj` takes the
  maximum standardised lens statistic per unit and applies maxT over units:
  one rule, valid across lenses and units. Measured family-wise error .030
  (flat) and .040 (combined) at K = 12 against the old rule's .060. USE
  `p_unit_adj` FOR NAMING; the per-lens columns remain for reporting which
  lens carried a finding.

# dmsa 1.4.0 (2026-08-15)

* THREE ENGINES; the COMBINED engine is now the DEFAULT. Pre-registered
  testing on the PTSD and children cohorts showed neither flat nor
  reliability dominates on real data: reliability gains on concentrated
  single-gene signals (FKBP5, SGK1) but loses on diffuse system/module
  signals (the GR-module PTSD lead, the child HPA-system finding). The
  combined engine fuses the two per lens by an ACAT (Cauchy) combination of
  their rank-p's on ONE shared permutation stream, so the family-wise maxT is
  an exact JOINT NULL - not a post-hoc ACAT of two adjusted p's. It keeps
  every real signal (validated: GR module .0445 flat / .118 reliability /
  .064 combined; FKBP5 .681 / .189 / .368) for the small toll of hedging the
  engine choice, in the same spirit as the three-lens triangulation. Both
  per-permutation regressions are shared across engines (computed once), so
  combined costs ~2x lens work, not 2x total. `weighting = c("combined",
  "reliability", "flat")` throughout; flat remains the exact reproducibility
  reference. Calibration of the combined omnibus verified (outcome-free
  weights + pooled-rank fusion keep type I at alpha). Frame print names the
  active engine. +2 assertions (combined calibration + both-archetype
  recovery).

# dmsa 1.3.0 (2026-08-15)

* TWO ENGINES, reliability weighting now the DEFAULT. Each aligned probe is
  weighted by its item-rest (centrality) correlation with the rest of its
  unit - `dmsa_relweights()`, computed from methylation ALONE so the
  permutation null stays exact and weighting changes power, not type-I error.
  In simulation it gains power on units whose probe quality is heterogeneous
  and ties the flat engine when it is uniform; it therefore ships as the
  default. `weighting = "flat"` restores the original equal-weight engine and
  is the exact reproducibility reference. Threaded through `dmsa_frame()`,
  `dmsa_report()`, `dmsa_triangulate()`, `dmsa_scores()` and `dmsa_change()`
  via `weighting` and `w_floor` (leading-eigenvalue coherence floor, default
  1.5). Guards: units with < 3 usable probes, no shared axis, or all-zero
  clamped weights fall back to flat. `dmsa_test()` gains a `w` argument.
  Frame print reports the active engine. +7 assertions (test-weight.R):
  calibration, flat reproducibility, power gain on heterogeneous units,
  guard fallbacks.

# dmsa 1.2.0 (2026-08-14)

* `dmsa_change()` - longitudinal mDMSA on the 2-wave identity: time x S x E
  == Delta-S x E. Per-unit aligned scores at each wave from the SAME map,
  Delta-score x exposure through the composite machinery (Freedman-Lane,
  family blocks), maxT/minP across the level-local family on one shared
  permutation stream; baseline score as regression-to-the-mean guard;
  `mod2` adds the second moderator (Delta-S x E x mod2). Validated on the
  Oct-7 cohort (n = 172) before shipping, per the pre-set condition.
* The frame layout `moderation`/`mod`/`mod2` now carries the three-way end
  to end (planted-interaction test locked in the suite); in a longitudinal
  frame Delta IS time, so `mod`/`mod2` are the two remaining factors.
* Frame fixes from the real child-data run: compositional check no longer
  trips on zero-variance aliased columns; `chip_f` autofix accepts any
  `chip`/`chip_T*` column; factor covariates and chip are droplevels-ed
  after complete-case subsetting (empty levels made designs singular).

# dmsa 1.1.2 (2026-08-14)

* A completed moderation run with zero survivors now SAYS so: print and
  summary report units tested, families, B, and the strongest interaction
  (unit, b, t, raw and adjusted p) - a calibrated null no longer looks like
  a run that did nothing. Sub-minute runtimes print in seconds.

# dmsa 1.1.1 (2026-08-14)

* Fix: the untestable-gene count in `frame$untestable` and the report Design
  notes was computed across ALL systems in the data instead of the chosen
  `systems` subset (spotted by comparing the Mac and cloud anchor runs -
  the statistics were never affected, only the coverage sentence).

# dmsa 1.1.0 (2026-08-14)

The two-call user interface.

* `dmsa_frame()` - declare everything once: data, levels (system / module /
  gene / probe), outcomes, covariates ("contract" default), random-effect
  blocks, moderation, map, styling. Runs the TEST DRIVE automatically
  (autofix = TRUE): numeric coercion with per-column loss report, QR rank
  audit dropping aliased columns by name (compositional groups flagged),
  block audit with equal-size permutation strata, confounded-block flags,
  hard error only when fewer than ~2^10 distinct permutations exist, binary
  outcome under gaussian autofixed to logistic, B = 49 pilot with a printed
  ETA. Every action lands in `frame$corrections`, printed by `print(frame)`
  and repeated in the report's Design notes.
* `dmsa_report(frame)` - runs `dmsa_triangulate()` at each requested level
  with correction ONLY inside level-local families (systems named; modules of
  the system; genes of the system; probes of surviving genes), writes
  figures/ tables/ summary.md under `frame$outdir`, and returns the paths.
* Moderation layout: `moderation = FALSE` (default), `mod = ""`, `mod2 = ""`.
  `moderation = TRUE` + `mod` tests frame x mod through the composite lens
  (the one lens that carries to products), maxT across units on a shared
  permutation stream; `mod2` declares frame x mod x mod2 - the 3-way engine
  ships in 1.2.0 after tonight's longitudinal test drive, and 1.1.0 raises a
  clear error if `mod2` is set.
* `frame_role = "predictor" | "outcome"` - methylation as predictor (default)
  or as the outcome of the named columns; the composite lens flips the score
  to the left-hand side, narration switches to expression-tone language.
* `cpg_map = "confidence" | "full"` - both alignments are always built; any
  gene or system whose pooled direction flips between maps is analysed under
  the chosen map but warned about, stored in `frame$map_conflicts`, and
  repeated in the report.
* Brain-bridge post-hoc, always on: surviving units' probes go through
  `cpgdirection::cpg_brain_bridge()`; hit rule fixed as `bridge_usable |
  has_t2_direction | !is.na(bridge_grade)`; `ensemble_pctile` is context.
  Errors never kill a report.
* `dmsa_triangulate()` now returns NA (not a p = 1 placeholder) for units
  with zero usable direction-called probes; `dmsa_frame()` tracks them in
  `frame$untestable` and the report counts them per family.
* `dmsa_model(nulls = TRUE)` returns its permutation null (`null_t`) for
  family maxT across separately fitted models.
* The former figure helper `dmsa_report()` is now `dmsa_report_panels()`.
* gt moved into Suggests: tables use gt when installed, otherwise a
  self-contained HTML fallback (CSVs are always written).

Anchor (real Alpha data, systems 1 + 21, B = 1999, seed 1, cpg_map = "full"):
AVP coherence adj .0120, PLAGL1 composite adj .0080 - identical to the
locked triangulation run. New from the always-on bridge: two PLAGL1 probes
carry a formal brain bridge (cg08263357, ensemble percentile 94; cg11532302).
