# w_g is the multiplier that decides which way a system score points. These
# tests protect the invariants that make a wrong sign visible rather than silent.

test_that("the bundled polarity table loads and covers the cascade exactly", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  p <- dmsa_polarity()
  expect_s3_class(p, "dmsa_polarity")
  expect_equal(nrow(p$polarity), 1234L)
  expect_equal(length(unique(p$polarity$system_id)), 30L)
  expect_true(all(p$polarity$w_g %in% c(-1, 0, 1)))
  ## one row per gene-system pair, matching the cascade with no gaps
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  cas <- dmsa_sets()$cascade
  want <- unique(paste(cas$system_id, cas$gene))
  got <- paste(p$polarity$system_id, p$polarity$gene)
  expect_setequal(want, got)
  expect_false(anyDuplicated(got) > 0)
})

test_that("an anchor may never lower its own system's tone", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  p <- dmsa_polarity()$polarity
  expect_true(all(p$w_g[p$anchor] > 0))
  ## and the validator catches it if one ever does
  bad <- data.frame(system_id = "1", gene = c("A", "B"), w_g = c(-1, 1),
                    anchor = c(TRUE, FALSE), stringsAsFactors = FALSE)
  chk <- dmsa_polarity_check(bad, verbose = FALSE)
  expect_false(chk$checks[["every anchor has w_g > 0"]]$ok)
})

test_that("w_g outside [-1, 1] and duplicate rows are refused at load", {
  expect_error(dmsa_polarity(data.frame(system_id = "1", gene = "A", w_g = 2)),
               "must lie in")
  expect_error(dmsa_polarity(data.frame(system_id = "1", gene = "A", w_g = "x")),
               "non-numeric")
  expect_error(dmsa_polarity(data.frame(system_id = c("1", "1"),
                                        gene = c("A", "A"), w_g = c(1, -1))),
               "duplicate")
  expect_error(dmsa_polarity(data.frame(gene = "A", w_g = 1)), "missing required")
})

test_that("every sign carries a source, and the grade is derived from it", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  p <- dmsa_polarity()$polarity
  expect_false(any(p$w_g_source == "unstated"))
  expect_setequal(unique(p$grade[p$w_g_source == "curated"]), "curated")
  expect_setequal(unique(p$grade[p$w_g_source == "rule"]), "heuristic")
  expect_setequal(unique(p$grade[p$w_g_source == "omnipath"]), "database")
  expect_setequal(unique(p$grade[p$w_g_source == "unresolved"]), "none")
  ## an unresolved row must never carry a sign - that is the whole point of it
  expect_true(all(p$w_g[p$w_g_source == "unresolved"] == 0))
})

test_that("the PI's curated signs survived the audit unchanged", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  f <- system.file("extdata", "alpha_polarity.csv", package = "dmsa")
  skip_if(!nzchar(f), "legacy curated table not installed")
  old <- utils::read.csv(f, stringsAsFactors = FALSE,
                         colClasses = c(system_id = "character"))
  new <- dmsa_polarity()$polarity
  m <- merge(old, new, by = c("system_id", "gene"))
  expect_true(nrow(m) > 100)
  expect_equal(as.numeric(m$w_g_draft), as.numeric(m$w_g))
  expect_true(all(m$w_g_source == "curated"))
})

test_that("review surfaces the rows a human has to decide, and only those", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  r <- dmsa_polarity_review()
  expect_s3_class(r, "dmsa_polarity_review")
  expect_true(nrow(r) > 0)
  expect_true(all(nzchar(r$review_flag)))
  adj <- dmsa_polarity_review(which = "disagreement")
  ## a "disagreement" is either a database contradiction or a framing choice only
  ## the analyst can settle; both need a human, neither is a defect
  expect_true(all(grepl("adjudicate|framing", adj$review_flag)))
  expect_true(nrow(adj) < nrow(r))
  ## a disagreement must say what the alternative reading is, or it is not
  ## actionable
  expect_true(all(nzchar(adj$evidence)))
  un <- dmsa_polarity_review(which = "unresolved")
  expect_true(all(un$w_g == 0))
  expect_output(print(r), "needing a decision")
})

test_that("the validator reports one-sided and anchorless systems rather than passing them", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  chk <- dmsa_polarity_check(verbose = FALSE)
  ## these are known, real properties of the Alpha panel: pharmacogenes have no
  ## tone at all, and the imprinted panel is a readout rather than an axis
  expect_false(chk$checks[["every system has at least one anchor"]]$ok)
  expect_false(chk$checks[["no system is entirely one-sided"]]$ok)
  expect_match(chk$checks[["no system is entirely one-sided"]]$detail, "21")
  ## a two-sided system with anchors passes both
  ok <- data.frame(system_id = "1", gene = c("A", "B", "C"), w_g = c(1, -1, 0),
                   anchor = c(TRUE, FALSE, FALSE), w_g_source = "curated",
                   stringsAsFactors = FALSE)
  c2 <- dmsa_polarity_check(ok, verbose = FALSE)
  expect_true(c2$checks[["every system has at least one anchor"]]$ok)
  expect_true(c2$checks[["no system is entirely one-sided"]]$ok)
})

test_that("polarity travels with the sets object and narrows with a selection", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  s <- dmsa_sets()
  expect_s3_class(s$polarity, "dmsa_polarity")
  expect_equal(nrow(s$polarity$polarity), 1234L)
  expect_output(print(s), "polarity: ")
  sel <- dmsa_select(systems = "hpa")
  expect_equal(nrow(sel$polarity$polarity), 49L)
  expect_true(all(sel$polarity$polarity$system_id == "2"))
  ## coverage is computed against the selection, so unsigned genes are visible
  expect_equal(sum(sel$polarity$coverage$n_genes), 49L)
  expect_equal(sum(sel$polarity$coverage$n_missing), 0L)
})

test_that("a cascade with no polarity is usable and says so", {
  d <- data.frame(system_id = "1", system = "S", module_id = "1.1",
                  module = "M", gene = c("G1", "G2"), cpg = c("cg1", "cg2"),
                  stringsAsFactors = FALSE)
  s <- dmsa_sets(d)
  expect_null(s$polarity)
  expect_output(print(s), "polarity: none attached")
  ## and a cascade carrying its own w_g column picks it up
  d2 <- d; d2$w_g <- c(1, -1); d2$w_g_source <- "curated"
  s2 <- dmsa_sets(d2)
  expect_s3_class(s2$polarity, "dmsa_polarity")
  expect_equal(s2$polarity$polarity$w_g, c(1, -1))
})

test_that("the draft fetcher never invents a sign when a resource is unreachable", {
  ## No network in the test environment, so every resource must fail closed:
  ## unresolved rows with w_g 0, not guesses, and needs_review on every row.
  d <- data.frame(system_id = "1", system = "S", module_id = "1.1", module = "M",
                  gene = c("G1", "G2"), cpg = c("cg1", "cg2"),
                  stringsAsFactors = FALSE)
  suppressWarnings(
    draft <- dmsa_polarity_fetch(dmsa_sets(d), sources = "go", quiet = TRUE))
  expect_equal(nrow(draft), 2L)
  expect_true(all(draft$needs_review))
  expect_true(all(draft$w_g == 0 | draft$w_g_source != "unresolved"))
  expect_true(all(draft$w_g[draft$w_g_source == "unresolved"] == 0))
  expect_error(dmsa_polarity_fetch(dmsa_sets(d), sources = "nope"),
               "must name at least one")
})

test_that("the source catalogue states each resource's caveat", {
  s <- dmsa_polarity_sources()
  expect_true(all(c("go", "omnipath", "trrust", "signor") %in% s$source))
  expect_true(all(nzchar(s$caveat)))
  ## the default-to-activating problem must be recorded against OmniPath, since
  ## that is the trap a user would otherwise walk into
  expect_match(s$caveat[s$source == "omnipath"], "default")
  expect_match(s$caveat[s$source == "trrust"], "Unknown")
})

test_that("every sign points at the system, not at its module's sub-process", {
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  p <- dmsa_polarity()$polarity
  ts <- .pol_toward_system(p)
  ## A role that commits to a direction may never contradict its own sign. This
  ## is the mechanical signature of a gene scored against its module instead of
  ## its system: PRL was labelled "driver" at w_g -1 because prolactin drives an
  ## endocrine axis while lowering the libido tone the system actually names.
  expect_length(ts$role_conflict, 0L)
  ## and the vocabulary is closed, so a new role cannot slip past the check
  expect_length(ts$unknown_roles, 0L)
  ## brake-of-brake exists precisely so a +1 on an inhibitor is expressible
  expect_true(all(p$w_g[tolower(p$role) == "brake-of-brake"] > 0))
})

test_that("roles whose sign is inherited are allowed to disagree with their word", {
  ## GNAI2 is an accessory at -1 because the D2 receptor it serves is -1, and
  ## TERT is capacity at -1 because more telomerase means less ageing tone.
  ## Neither is an error, and the check must not fire on them.
  d <- data.frame(system_id = "1", gene = c("A", "B", "C"),
                  w_g = c(-1, -1, 1), role = c("accessory", "capacity", "driver"),
                  module_id = "1.1", module = "M", stringsAsFactors = FALSE)
  ts <- .pol_toward_system(d)
  expect_length(ts$role_conflict, 0L)
  ## but a driver at -1 is caught
  d$w_g[3] <- -1
  expect_length(.pol_toward_system(d)$role_conflict, 1L)
})

test_that("the module-label check uses word boundaries", {
  ## "sex deTERMINATION" must not match "termination", or the check cries wolf on
  ## every gonadal module.
  d <- data.frame(system_id = "5", gene = c("A", "B"), w_g = c(1, 1),
                  role = "driver", module_id = "5.6",
                  module = "Gonadal sex determination and differentiation",
                  stringsAsFactors = FALSE)
  expect_length(.pol_toward_system(d)$module_flip, 0L)
  ## a genuinely brake-labelled module with only positive signs does flip
  d$module <- "Dopamine reuptake and catabolism"
  expect_length(.pol_toward_system(d)$module_flip, 1L)
})

test_that("review accepts the object the user is already holding", {
  skip_if(is.null(.cas_builtin("cascade")), "bundled cascade not installed")
  skip_if(is.null(.pol_builtin()), "bundled polarity not installed")
  ## The print banner says "dmsa_polarity_review()", so the obvious argument is
  ## the sets or selection in hand. Passing one used to fail with an unhelpful
  ## "cannot coerce class dmsa_selection to a data.frame".
  sel <- dmsa_select(systems = c("hpa", "oxytocin"))
  r <- dmsa_polarity_review(sel)
  expect_s3_class(r, "dmsa_polarity_review")
  expect_setequal(unique(r$system_id), c("1", "2"))
  ## and it is scoped: fewer rows than the whole table
  expect_lt(nrow(r), nrow(dmsa_polarity_review("alpha")))
  expect_s3_class(dmsa_polarity_review(dmsa_sets()), "dmsa_polarity_review")
  expect_s3_class(dmsa_polarity(sel), "dmsa_polarity")
  ## a cascade with no polarity says so instead of coercing
  d <- data.frame(system_id = "1", system = "S", module_id = "1.1", module = "M",
                  gene = "G1", cpg = "cg1", stringsAsFactors = FALSE)
  expect_error(dmsa_polarity(dmsa_sets(d)), "no polarity table is attached")
})
