# The selection cascade: system > module > gene > probe.
# These tests protect three things that would fail silently if broken:
# module ids that look numeric, short-name resolution, and the default that
# everything below a chosen system is taken in full.

test_that("the bundled cascade loads with its declared shape", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  s <- dmsa_sets()
  expect_s3_class(s, "dmsa_sets")
  expect_equal(nrow(s$systems), 30L)
  expect_equal(nrow(s$modules), 188L)
  expect_equal(length(unique(s$cascade$gene)), 1234L)
  expect_equal(length(unique(s$cascade$cpg)), 16823L)
  ## every module belongs to exactly one system, every gene to one module
  expect_true(all(tapply(s$cascade$system_id, s$cascade$module_id,
                         function(z) length(unique(z))) == 1))
  expect_true(all(tapply(s$cascade$module_id, s$cascade$gene,
                         function(z) length(unique(z))) == 1))
})

test_that("module ids stay text: 24.10 is not 24.1", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  s <- dmsa_sets()
  expect_type(s$cascade$module_id, "character")
  ids <- unique(s$cascade$module_id)
  expect_true(all(c("24.1", "24.10", "24.11", "24.13") %in% ids))
  ## and they sort in human order, not lexical order
  im <- s$modules$module_id[s$modules$system_id == "24"]
  expect_equal(im[1:2], c("24.1", "24.2"))
  expect_equal(utils::tail(im, 1), "24.13")
})

test_that("systems resolve by short name, id, full name and unique prefix", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  f <- function(q) dmsa_select(systems = q)$systems$system_short
  expect_equal(f("hpa"), "hpa")
  expect_equal(f("HPA"), "hpa")                       # case-insensitive
  expect_equal(f(2), "hpa")                           # numeric id
  expect_equal(f("24"), "immune")                     # id as text
  expect_equal(f("Immune, inflammation & HLA"), "immune")
  expect_equal(f("oxy"), "oxytocin")                  # unique prefix
  expect_setequal(f(c("hpa", "oxytocin")), c("hpa", "oxytocin"))
})

test_that("ambiguous and unmatched system names error with the candidates named", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  expect_error(dmsa_select(systems = "ox"), "ambiguous")
  expect_error(dmsa_select(systems = "ox"), "oxidative")
  expect_error(dmsa_select(systems = "cortisol"), "did not match")
  expect_error(dmsa_select(systems = "cortisol"), "available short names")
})

test_that("everything below the chosen system defaults to full", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  sel <- dmsa_select(systems = "hpa")
  all_hpa <- dmsa_sets()$cascade
  all_hpa <- all_hpa[all_hpa$system_short == "hpa", ]
  expect_equal(nrow(sel$cascade), nrow(all_hpa))
  expect_true(all(sel$full))
  expect_equal(length(unique(sel$cascade$module_id)), 8L)
})

test_that("narrowing a level restricts and recounts honestly", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  sel <- dmsa_select(systems = "hpa", modules = c("2.6", "2.8"))
  expect_equal(sort(unique(sel$cascade$module_id)), c("2.6", "2.8"))
  expect_false(sel$full["modules"])
  expect_true(sel$full["genes"])
  expect_equal(sum(sel$modules$n_cpgs_selected),
               length(unique(sel$cascade$cpg)))
  g1 <- unique(sel$cascade$gene)[1]
  one <- dmsa_select(systems = "hpa", genes = g1)
  expect_equal(unique(one$cascade$gene), g1)
  expect_error(dmsa_select(systems = "hpa", modules = "99.9"), "did not match")
})

test_that("columns = 'auto' keys the cascade to the user's matrix", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  sel0 <- dmsa_select(systems = "pharmacogenes")
  d <- as.data.frame(matrix(0, 2, nrow(sel0$cascade)))
  names(d) <- sel0$cascade$col_child_T4
  sel <- dmsa_select(systems = "pharmacogenes", columns = "auto", data = d)
  expect_equal(sel$column_key, "col_child_T4")
  expect_error(dmsa_select(systems = "pharmacogenes", columns = "auto"),
               "needs `data`")
  expect_error(dmsa_select(systems = "pharmacogenes", columns = "col_nope"),
               "not a column")
})

test_that("module evidence is carried and can be filtered to the caveats", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  sel <- dmsa_select(systems = "immune")
  ev <- dmsa_evidence(sel)
  expect_s3_class(ev, "dmsa_evidence")
  expect_equal(nrow(ev), 13L)
  expect_true(all(c("evidence_strength", "citation_keys", "audit_status") %in%
                    names(ev)))
  fl <- dmsa_evidence(sel, which = "flagged")
  expect_true(nrow(fl) > 0 && nrow(fl) < nrow(ev))
  expect_true(all(fl$evidence_strength %in% "Moderate" |
                    grepl("heterogeneous|measurement_defined", fl$audit_status)))
  ## the citations are what a reader is meant to check - they must be present
  expect_true(all(nzchar(fl$citation_keys)))
  ## and the banner names the counts rather than hiding them
  expect_output(print(sel), "module evidence: ")
})

test_that("a user cascade round-trips through the template and validator", {
  tmp <- file.path(tempdir(), "sets_demo.csv")
  on.exit(unlink(c(tmp, sub("\\.csv$", ".md", tmp))), add = TRUE)
  expect_message(dmsa_sets_template(tmp), "wrote")
  expect_true(file.exists(sub("\\.csv$", ".md", tmp)))
  chk <- dmsa_sets_check(tmp, verbose = FALSE)
  expect_true(chk$ok)
  s <- dmsa_sets(tmp)
  expect_equal(nrow(s$systems), 1L)
  expect_equal(s$systems$system_short, "hpa")
  expect_equal(s$modules$evidence_strength, "High")
  sel <- dmsa_select(s, systems = "hpa")
  expect_equal(nrow(sel$cascade), 2L)
})

test_that("the validator rejects a gene that sits in two modules", {
  bad <- data.frame(
    system_id = "1", system = "S", module_id = c("1.1", "1.2"),
    module = c("A", "B"), gene = c("G1", "G1"),
    cpg = c("cg1", "cg2"), stringsAsFactors = FALSE)
  chk <- dmsa_sets_check(bad, verbose = FALSE)
  expect_false(chk$ok)
  expect_false(chk$checks[["each gene in one module"]]$ok)
  ## a missing required column stops before the structural checks
  chk2 <- dmsa_sets_check(bad[, c("system_id", "gene", "cpg")], verbose = FALSE)
  expect_false(chk2$ok)
  expect_false(chk2$checks[["required columns"]]$ok)
})

test_that("short names are derived when a user cascade omits them", {
  d <- data.frame(
    system_id = c("1", "2"),
    system = c("HPA axis & glucocorticoid signalling", "Dopamine"),
    module_id = c("1.1", "2.1"), module = c("A", "B"),
    gene = c("G1", "G2"), cpg = c("cg1", "cg2"), stringsAsFactors = FALSE)
  s <- dmsa_sets(d)
  expect_setequal(s$systems$system_short, c("hpa_axis", "dopamine"))
  expect_equal(dmsa_select(s, systems = "dopamine")$cascade$gene, "G2")
})

test_that("dmsa_frame accepts short names and carries the evidence", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  f <- system.file("extdata", "coverage_v4_full.csv", package = "dmsa")
  skip_if(!nzchar(f), "bundled direction map not installed")
  mp <- utils::read.csv(f, stringsAsFactors = FALSE)
  mp <- mp[is.finite(mp$best_direction) & mp$system_id %in% c(1, 2), ]
  skip_if(nrow(mp) < 10, "too few mapped probes")
  set.seed(11); n <- 60
  d <- data.frame(out1 = stats::rnorm(n), sex = rep(1:2, length.out = n),
                  cov1 = stats::rnorm(n), chip = rep(1:3, length.out = n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cc in mp$column) d[[cc]] <- stats::plogis(stats::rnorm(n))
  fr <- dmsa_frame(d, outcome = "out1", covariates = c("sex", "cov1"),
                   systems = c("hpa", "oxytocin"), module = TRUE, B = 49,
                   plots = FALSE, tables = FALSE, summary = FALSE, seed = 1)
  expect_equal(nrow(fr$systems), 2L)
  expect_true(!is.null(fr$module_evidence) && nrow(fr$module_evidence) > 0)
  expect_s3_class(dmsa_evidence(fr), "dmsa_evidence")
  expect_output(print(fr), "module evidence")
  ## sets = NULL keeps the legacy path available
  fr2 <- dmsa_frame(d, outcome = "out1", covariates = c("sex", "cov1"),
                    systems = 2, sets = NULL, B = 49, plots = FALSE,
                    tables = FALSE, summary = FALSE, seed = 1)
  expect_equal(nrow(fr2$systems), 1L)
  expect_null(fr2$module_evidence)
})

test_that("the bundled cascade passes its own validator", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  ## This is the data sanity gate a user is told to run before publishing, so a
  ## FAIL on the shipped set is a defect in the check, not in the data. It used
  ## to report two: EPIC v2 replicate probe designs read as duplicate rows, and
  ## module evidence read as absent because it lives in the companion audit file.
  chk <- dmsa_sets_check("alpha", verbose = FALSE)
  expect_true(chk$ok)
  expect_length(chk$warnings, 0L)
  expect_true(chk$checks[["no duplicate module+gene+probe_id rows"]]$ok)
  expect_true(chk$checks[["module evidence annotated"]]$ok)
  ## the replicate designs are reported rather than silently collapsed
  expect_equal(chk$replicate_cpgs, 115L)
  expect_equal(chk$shared_cpgs, 228L)
  ## and the same holds for an already-loaded object
  expect_true(dmsa_sets_check(dmsa_sets(), verbose = FALSE)$ok)
})

test_that("a cascade with no module evidence is usable, with a caveat", {
  ## Absent evidence changes how a result reads; it does not make the set
  ## unusable, and reporting it as a failure would train users to ignore
  ## failures.
  d <- data.frame(system_id = "1", system = "S", module_id = "1.1", module = "M",
                  gene = c("G1", "G2"), cpg = c("cg1", "cg2"),
                  stringsAsFactors = FALSE)
  chk <- dmsa_sets_check(d, verbose = FALSE)
  expect_true(chk$ok)
  expect_equal(chk$warnings, "module evidence annotated")
  expect_false(chk$checks[["module evidence annotated"]]$ok)
  expect_equal(chk$checks[["module evidence annotated"]]$level, "warn")
  ## with no probe_id column the key falls back to cpg
  expect_true("no duplicate module+gene+cpg rows" %in% names(chk$checks))
  ## a real duplicate is still a failure
  chk2 <- dmsa_sets_check(rbind(d, d), verbose = FALSE)
  expect_false(chk2$ok)
})

test_that("replicate probe designs for one CpG are not duplicates", {
  ## cg1 measured by two EPIC v2 designs is two columns of data, not a repeated
  ## row. Keying uniqueness on cpg rejected 238 legitimate rows in the bundled
  ## set.
  d <- data.frame(system_id = "1", system = "S", module_id = "1.1", module = "M",
                  gene = "G1", cpg = c("cg1", "cg1"),
                  probe_id = c("cg1_BC21", "cg1_BC22"), stringsAsFactors = FALSE)
  chk <- dmsa_sets_check(d, verbose = FALSE)
  expect_true(chk$checks[["no duplicate module+gene+probe_id rows"]]$ok)
  expect_equal(chk$replicate_cpgs, 1L)
  ## but the same probe twice is
  d2 <- d; d2$probe_id <- "cg1_BC21"
  expect_false(dmsa_sets_check(d2, verbose = FALSE)$ok)
})

test_that("dmsa_systems accepts a bare pattern as its first argument", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  ## The documented example dmsa_systems("stress|hpa") used to error with
  ## "cascade file not found: stress|hpa" - the string was taken as a file path.
  r <- dmsa_systems("stress|hpa")
  expect_setequal(r$system_short, c("hpa", "oxidative"))
  expect_equal(dmsa_systems("hpa")$system_short, "hpa")
  ## and the object forms still work
  expect_equal(nrow(dmsa_systems()), 30L)
  expect_equal(nrow(dmsa_systems(dmsa_sets())), 30L)
  expect_equal(nrow(dmsa_systems(dmsa_select(systems = "hpa"))), 1L)
  ## "alpha" is still the cascade, not a pattern
  expect_equal(nrow(dmsa_systems("alpha")), 30L)
  ## an unmatched pattern messages and returns no rows rather than erroring
  expect_message(out <- dmsa_systems("zzz"), "no system matches")
  expect_equal(nrow(out), 0L)
})
