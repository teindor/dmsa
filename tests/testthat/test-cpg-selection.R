## PI request (2026-08-29): CpG-set selection by name pattern, for files that
## carry several biological groups (Alpha child file: "_maternal" /
## "_paternal" tag the parents' probes; the child's carry no tag).

.grp_pairs <- function() data.frame(
  cpg_id = c("cg0000001", "cg0000002", "cg0000003"),
  target_gene = c("OXTR", "NR3C1", "FKBP5"),
  best_direction = c(-1, 1, -1),
  probability_plus1 = c(.05, .9, .1),
  usable = TRUE, abstain_reason = NA,
  best_evidence = "measured", direction_tier = "M",
  stringsAsFactors = FALSE)

.grp_data <- function(n = 60) {
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  maternal_age = stats::rnorm(n),      # NOT a cpg column
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in c("cg0000001_BC21_T4_OXTR",
               "cg0000001_BC21_maternal_T1_OXTR",
               "cg0000001_BC21_paternal_T1_OXTR",
               "cg0000002_TC21_T4_NR3C1",
               "cg0000002_TC21_maternal_T1_NR3C1",
               "cg0000003_BC21_paternal_T1_FKBP5"))
    d[[cl]] <- stats::plogis(stats::rnorm(n))
  d
}

test_that("dmsa_cpg_columns: 'not the parents' is one line", {
  d <- .grp_data()
  child <- dmsa_cpg_columns(d, exclude = c("_maternal", "_paternal"))
  expect_setequal(child, c("cg0000001_BC21_T4_OXTR", "cg0000002_TC21_T4_NR3C1"))
  mom <- dmsa_cpg_columns(d, include = "_maternal")
  expect_setequal(mom, c("cg0000001_BC21_maternal_T1_OXTR",
                         "cg0000002_TC21_maternal_T1_NR3C1"))
  ## include then exclude compose
  par_no_dad <- dmsa_cpg_columns(d, include = c("_maternal", "_paternal"),
                                 exclude = "_paternal")
  expect_setequal(par_no_dad, mom)
  ## non-CpG columns are never candidates, even when a pattern matches them
  expect_false("maternal_age" %in%
                 dmsa_cpg_columns(d, include = "maternal"))
  ## works on a bare character vector and a matrix too
  expect_setequal(dmsa_cpg_columns(names(d),
                                   exclude = c("_maternal", "_paternal")),
                  child)
})

test_that("dmsa_frame(cpgs_exclude=) analyses only the child's probes", {
  set.seed(1)
  d <- .grp_data()
  old <- options(dmsa.pair_table = .grp_pairs()); on.exit(options(old))
  f <- dmsa_frame(d, direction_source = "cpgdirection",
                  cpgs_exclude = c("_maternal", "_paternal"),
                  outcome = "y", covariates = "cov1", blocks = "cID",
                  chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
                  summary = FALSE, progress = FALSE, beep = FALSE,
                  outdir = tempfile("dmsa_grp"))
  ## only the child's two CpGs are measurements; FKBP5 exists only as a
  ## paternal probe here and must be gone
  expect_setequal(unique(f$map$probe), c("cg0000001", "cg0000002"))
  expect_false("FKBP5" %in% f$map$gene)
  expect_true(any(grepl("cpg selection", f$corrections$field)))
  ## the covariate survives untouched
  expect_true("maternal_age" %in% names(f$data))
})

test_that("dmsa_frame(cpgs_include=) selects one parent's probes", {
  set.seed(1)
  d <- .grp_data()
  old <- options(dmsa.pair_table = .grp_pairs()); on.exit(options(old))
  f <- dmsa_frame(d, direction_source = "cpgdirection",
                  cpgs_include = "_maternal",
                  outcome = "y", covariates = "cov1", blocks = "cID",
                  chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
                  summary = FALSE, progress = FALSE, beep = FALSE,
                  outdir = tempfile("dmsa_grpM"))
  expect_setequal(unique(f$map$probe), c("cg0000001", "cg0000002"))
  expect_setequal(unique(f$map$gene), c("OXTR", "NR3C1"))
  ## and the measurement columns are the MATERNAL designs
  expect_true(all(grepl("_maternal",
                        unlist(f$measurements$replicates))))
})

test_that("patterns that remove everything are a hard error", {
  set.seed(1)
  d <- .grp_data()
  old <- options(dmsa.pair_table = .grp_pairs()); on.exit(options(old))
  expect_error(
    dmsa_frame(d, direction_source = "cpgdirection",
               cpgs_include = "_grandparental",
               outcome = "y", covariates = "cov1", blocks = "cID",
               chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
               summary = FALSE, progress = FALSE, beep = FALSE,
               outdir = tempfile("dmsa_grpE")),
    "removed every CpG-site column")
})

test_that("the bundled path honours the same filter", {
  set.seed(1)
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], 3)
  gene_tag <- "_FKBP5"                       # every FKBP5 column carries it
  n <- 60
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  f <- dmsa_frame(d, methylation = cols, direction_source = "bundled",
                  cpgs_include = gene_tag,
                  outcome = "y", covariates = "cov1", blocks = "cID",
                  chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
                  summary = FALSE, progress = FALSE, beep = FALSE,
                  outdir = tempfile("dmsa_grpB"))
  expect_equal(nrow(f$map), 3L)
  expect_error(
    dmsa_frame(d, methylation = cols, direction_source = "bundled",
               cpgs_exclude = gene_tag,
               outcome = "y", covariates = "cov1", blocks = "cID",
               chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
               summary = FALSE, progress = FALSE, beep = FALSE,
               outdir = tempfile("dmsa_grpB2")),
    "cpgs_include/cpgs_exclude")
})
