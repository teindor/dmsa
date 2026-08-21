# Package index

## The two-call interface

Declare the design once, run it once.

- [`dmsa_frame()`](https://teindor.github.io/dmsa/reference/dmsa_frame.md)
  : Declare a DMSA analysis frame
- [`dmsa_report()`](https://teindor.github.io/dmsa/reference/dmsa_report.md)
  : Run a declared frame and deliver figures, tables and a summary
- [`dmsa_report_panels()`](https://teindor.github.io/dmsa/reference/dmsa_report_panels.md)
  : Build the whole report: panel 1 once, panels 2 and 3 per firing
  system

## Bioconductor interoperability

Working from a SummarizedExperiment.

- [`dmsa_se()`](https://teindor.github.io/dmsa/reference/dmsa_se.md) :
  Split a SummarizedExperiment into the pieces DMSA tests

## Alignment and testing

The directional chain and the three lenses.

- [`dmsa_align()`](https://teindor.github.io/dmsa/reference/dmsa_align.md)
  : Build the aligned sign s_j for a set of probes
- [`dmsa_directions()`](https://teindor.github.io/dmsa/reference/dmsa_directions.md)
  : The direction calls DMSA ships with
- [`dmsa_triangulate()`](https://teindor.github.io/dmsa/reference/dmsa_triangulate.md)
  : Three complementary lenses on the same family
- [`dmsa_test()`](https://teindor.github.io/dmsa/reference/dmsa_test.md)
  : Fixed-sign or expected-sign DMSA pooled test
- [`dmsa_pdance()`](https://teindor.github.io/dmsa/reference/dmsa_pdance.md)
  : The p-dance test: does a finding's p-value survive a change of
  analysis set?
- [`dmsa_score()`](https://teindor.github.io/dmsa/reference/dmsa_score.md)
  : Subject-level aligned tone score
- [`dmsa_scores()`](https://teindor.github.io/dmsa/reference/dmsa_scores.md)
  : Aligned set scores, and the alternatives worth comparing against
- [`dmsa_relweights()`](https://teindor.github.io/dmsa/reference/dmsa_relweights.md)
  : Outcome-free reliability weights for aligned probes
- [`dmsa_perm_pvalue()`](https://teindor.github.io/dmsa/reference/dmsa_perm_pvalue.md)
  : Permutation p-value for a DMSA statistic
- [`dmsa_omnibus()`](https://teindor.github.io/dmsa/reference/dmsa_omnibus.md)
  : Directionless omnibus test of a categorical focal exposure over a
  set

## Selection and hierarchy

- [`dmsa_sets()`](https://teindor.github.io/dmsa/reference/dmsa_sets.md)
  : Load a selection cascade (system \> module \> gene \> probe)
- [`dmsa_select()`](https://teindor.github.io/dmsa/reference/dmsa_select.md)
  : Select a set of units from a cascade
- [`dmsa_cascade()`](https://teindor.github.io/dmsa/reference/dmsa_cascade.md)
  : Gated hierarchical DMSA with level-specific inference
- [`dmsa_systems()`](https://teindor.github.io/dmsa/reference/dmsa_systems.md)
  : The systems available in a cascade, with the short names to select
  them by
- [`dmsa_levels()`](https://teindor.github.io/dmsa/reference/dmsa_levels.md)
  : Gated hierarchical DMSA on the pooled aligned statistic
- [`dmsa_tree()`](https://teindor.github.io/dmsa/reference/dmsa_tree.md)
  : Build a cascade tree for a set of probes from a reference bundle
- [`dmsa_sets_check()`](https://teindor.github.io/dmsa/reference/dmsa_sets_check.md)
  : Validate a candidate cascade before using it
- [`dmsa_sets_template()`](https://teindor.github.io/dmsa/reference/dmsa_sets_template.md)
  : Write a cascade template and the instructions for building one
- [`dmsa_cascade_null()`](https://teindor.github.io/dmsa/reference/dmsa_cascade_null.md)
  : Null statistics for cascade calibration

## Polarity and reference maps

- [`dmsa_polarity()`](https://teindor.github.io/dmsa/reference/dmsa_polarity.md)
  : Load a gene-to-system polarity table
- [`dmsa_polarity_check()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_check.md)
  : Validate a candidate polarity table
- [`dmsa_polarity_for()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_for.md)
  : Polarity table for one system, in the form dmsa_align() wants
- [`dmsa_polarity_fetch()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_fetch.md)
  : Draft polarity signs for a panel from public databases
- [`dmsa_polarity_review()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_review.md)
  : Rows of a polarity table that need a human decision
- [`dmsa_polarity_sources()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_sources.md)
  : Which public resources DMSA can draft polarity from, and what each
  is good for
- [`dmsa_reference()`](https://teindor.github.io/dmsa/reference/dmsa_reference.md)
  : Construct a DMSA reference bundle
- [`dmsa_reference_csv()`](https://teindor.github.io/dmsa/reference/dmsa_reference_csv.md)
  : Read a user-supplied system / module / gene csv as a reference
  bundle
- [`dmsa_reference_read()`](https://teindor.github.io/dmsa/reference/dmsa_reference_read.md)
  : Read a reference bundle from CSV files
- [`dmsa_reference_template()`](https://teindor.github.io/dmsa/reference/dmsa_reference_template.md)
  : Write a template csv for a user-supplied reference
- [`dmsa_reference_write()`](https://teindor.github.io/dmsa/reference/dmsa_reference_write.md)
  : Write a reference bundle to a directory
- [`dmsa_evidence()`](https://teindor.github.io/dmsa/reference/dmsa_evidence.md)
  : Module-level evidence and citations for a cascade, selection or
  frame
- [`alpha_polarity()`](https://teindor.github.io/dmsa/reference/alpha_polarity.md)
  : Alpha panel gene -\> system-activation polarity (bundled, DRAFT)
- [`alpha_gene_systems()`](https://teindor.github.io/dmsa/reference/alpha_gene_systems.md)
  : Alpha panel gene -\> system map (bundled)
- [`alpha_reference()`](https://teindor.github.io/dmsa/reference/alpha_reference.md)
  : The bundled Project Alpha tables as a reference bundle
- [`alpha_design()`](https://teindor.github.io/dmsa/reference/alpha_design.md)
  : Project Alpha build contracts

## Design, models and extensions

- [`dmsa_design()`](https://teindor.github.io/dmsa/reference/dmsa_design.md)
  : Declare a DMSA design: covariates, dependence, and permutation
  blocks
- [`dmsa_check_design()`](https://teindor.github.io/dmsa/reference/dmsa_check_design.md)
  : Validate a design against the data before fitting
- [`dmsa_model()`](https://teindor.github.io/dmsa/reference/dmsa_model.md)
  : One model, block-permutation calibrated, for any term
- [`dmsa_model_omnibus()`](https://teindor.github.io/dmsa/reference/dmsa_model_omnibus.md)
  : Omnibus across the score projections
- [`dmsa_fit()`](https://teindor.github.io/dmsa/reference/dmsa_fit.md) :
  Fit, pool and permute a DMSA set test under a declared design
- [`dmsa_moderate()`](https://teindor.github.io/dmsa/reference/dmsa_moderate.md)
  : Moderated DMSA (mDMSA)
- [`dmsa_change()`](https://teindor.github.io/dmsa/reference/dmsa_change.md)
  : Two-wave change-score DMSA: time x methylation x factor
- [`dmsa_outcome()`](https://teindor.github.io/dmsa/reference/dmsa_outcome.md)
  : Test an aligned tone score against a non-Gaussian outcome
- [`dmsa_contrasts()`](https://teindor.github.io/dmsa/reference/dmsa_contrasts.md)
  : Expand a categorical focal exposure into one directional design per
  contrast
- [`dmsa_contrast_adjust()`](https://teindor.github.io/dmsa/reference/dmsa_contrast_adjust.md)
  : Correct p-values across the contrasts of one categorical exposure
- [`dmsa_mvalues()`](https://teindor.github.io/dmsa/reference/dmsa_mvalues.md)
  : M-value transform with clamping
- [`dmsa_gate()`](https://teindor.github.io/dmsa/reference/dmsa_gate.md)
  : Gated hierarchical DMSA: one score, one model, level-local
  correction
- [`dmsa_balance()`](https://teindor.github.io/dmsa/reference/dmsa_balance.md)
  : Per-system coverage and polarity balance

## Multiplicity

- [`dmsa_lfdr()`](https://teindor.github.io/dmsa/reference/dmsa_lfdr.md)
  : Local false discovery rate by central matching
- [`dmsa_bfdr_select()`](https://teindor.github.io/dmsa/reference/dmsa_bfdr_select.md)
  : Select units by Bayesian FDR

## Figures

- [`dmsa_plot_genes()`](https://teindor.github.io/dmsa/reference/dmsa_plot_genes.md)
  : Panel 2: gene level inside one system
- [`dmsa_plot_locus()`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus.md)
  : Locus figure: chromosome, gene with CpG lollipops, and per-CpG
  effects
- [`dmsa_plot_locus_gviz()`](https://teindor.github.io/dmsa/reference/dmsa_plot_locus_gviz.md)
  : The locus panel drawn by Gviz
- [`dmsa_plot_probes()`](https://teindor.github.io/dmsa/reference/dmsa_plot_probes.md)
  : Panel 3: probe level, with the CpG-to-expression direction applied
- [`dmsa_plot_systems()`](https://teindor.github.io/dmsa/reference/dmsa_plot_systems.md)
  : Panel 1: dense arm against sparse arm, every system
- [`dmsa_gene_model()`](https://teindor.github.io/dmsa/reference/dmsa_gene_model.md)
  : A gene's exon, intron and UTR structure on real coordinates
- [`dmsa_gene_model_check()`](https://teindor.github.io/dmsa/reference/dmsa_gene_model_check.md)
  : Check a gene model before drawing from it

## Coordinates and annotation

- [`dmsa_probe_coords()`](https://teindor.github.io/dmsa/reference/dmsa_probe_coords.md)
  : Probe genomic coordinates, without Bioconductor
- [`dmsa_probe_annotation_template()`](https://teindor.github.io/dmsa/reference/dmsa_probe_annotation_template.md)
  : Write a template script that fetches probe coordinates from
  Bioconductor
- [`dmsa_write_coords()`](https://teindor.github.io/dmsa/reference/dmsa_write_coords.md)
  : Save a coordinate table so the lookup never has to be repeated
- [`dmsa_test_drive()`](https://teindor.github.io/dmsa/reference/dmsa_test_drive.md)
  : Re-run the frame's test drive
