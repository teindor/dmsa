## Spec 26/28: every statistic column is labelled for what it is - and the
## prose says in so many words that only the per-lens family-adjusted p
## names a finding, while the omnibus and the joint union are synthesis.

test_that("spec 26/28: summary.md carries the column glossary with the roles stated", {
  set.seed(21)
  pairs <- data.frame(cpg_id = sprintf("cg%07d", 1:3),
                      target_gene = c("NR3C1", "FKBP5", "CRH"),
                      best_direction = c(-1, 1, 1),
                      probability_plus1 = c(.05, .9, .9), usable = TRUE,
                      abstain_reason = NA_character_, best_evidence = "smr_high",
                      direction_tier = "S1", stringsAsFactors = FALSE)
  n <- 80
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n), cID = rep(1:40, each = 2))
  for (cl in pairs$cpg_id) d[[cl]] <- stats::plogis(stats::rnorm(n))
  od <- tempfile("dmsa_gloss")
  old <- options(dmsa.pair_table = pairs); on.exit(options(old))
  f <- suppressMessages(dmsa_frame(
    d, methylation = pairs$cpg_id, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = FALSE, B = 49, progress = FALSE, beep = FALSE, outdir = od))
  suppressWarnings(suppressMessages(dmsa_report(f)))
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  expect_true(any(grepl("^### The statistic columns, and which of them name a finding", sm)))
  ## every exported p column is in the table
  for (cn in c("p_coherence_adj", "p_union_exact", "p_omnibus", "p_omnibus_adj",
               "p_unit_adj", "fwer_realized", "selected"))
    expect_true(any(grepl(paste0("`", cn, "`"), sm, fixed = TRUE)), info = cn)
  ## the roles, in words
  expect_true(any(grepl("YES - the naming statistic", sm, fixed = TRUE)))
  expect_true(any(grepl("ANY of the three is below alpha", sm, fixed = TRUE)))
  expect_true(any(grepl("`selected` IS the verdict", sm, fixed = TRUE)))
  expect_true(any(grepl("never a verdict on its own", sm, fixed = TRUE)))
  expect_true(any(grepl("it does not gate naming", sm, fixed = TRUE)))
  expect_true(any(grepl("quote them, do not", sm, fixed = TRUE)))
  ## the glossary and the units table agree on what exists
  u <- utils::read.csv(file.path(od, "tables", "units.csv"))
  for (cn in c("p_coherence_adj", "p_union_exact", "p_omnibus", "p_omnibus_adj",
               "p_unit_adj", "fwer_realized", "selected", "best_lens"))
    expect_true(cn %in% names(u), info = cn)
})
