## Spec 43: requested exchangeability blocks can never silently disappear.
## Spec 45: categorical covariates stay categorical.

.saf_fixture <- function(n = 60, k = 3) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], k)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols, n = n)
}

.saf_run <- function(fx, covs = "cov1", re = "cID", data = fx$data)
  dmsa_frame(data, methylation = fx$cols, direction_source = "bundled",
             outcome = "y", covariates = covs,
             random_effects = re, chip = FALSE, B = 19, plots = FALSE,
             tables = FALSE, summary = FALSE, progress = FALSE, beep = FALSE,
             outdir = tempfile("dmsa_saf"))

test_that("a missing block column hard-errors even under autofix", {
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data; d$cID <- NULL
  e <- expect_error(.saf_run(fx, data = d))
  expect_match(conditionMessage(e), "not found in `data`")
  expect_match(conditionMessage(e), "unrestricted permutation")
  expect_match(conditionMessage(e), "No analysis was run")
})

test_that("all-singleton blocks hard-error rather than degenerate", {
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data; d$cID <- seq_len(fx$n)
  e <- expect_error(.saf_run(fx, data = d))
  expect_match(conditionMessage(e), "no permutable groups")
  expect_match(conditionMessage(e), "unrestricted permutation")
})

test_that("a valid block still runs, and no blocking is legitimate", {
  set.seed(1); fx <- .saf_fixture()
  expect_s3_class(.saf_run(fx), "dmsa_frame")
  expect_s3_class(.saf_run(fx, re = NULL), "dmsa_frame")
})

test_that("word-level categorical covariates are kept as factors", {
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data
  d$grp <- rep(c("control", "treated", "placebo"), length.out = fx$n)
  f <- .saf_run(fx, covs = c("cov1", "grp"), data = d)
  expect_s3_class(f$data$grp, "factor")
  expect_equal(nlevels(f$data$grp), 3L)
})

test_that("an explicit factor with numeric-looking levels stays categorical", {
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data
  d$site <- factor(rep(c("1", "2", "3"), length.out = fx$n))
  f <- .saf_run(fx, covs = c("cov1", "site"), data = d)
  expect_s3_class(f$data$site, "factor")   # NOT coerced to 1,2,3 continuous
  expect_equal(nlevels(f$data$site), 3L)
})

test_that("numbers stored as text are still coerced to numeric", {
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data
  d$numtxt <- as.character(round(stats::rnorm(fx$n), 3))
  f <- .saf_run(fx, covs = c("cov1", "numtxt"), data = d)
  expect_type(f$data$numtxt, "double")
})

test_that("a single-level categorical covariate is an explicit error", {
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data; d$const <- "onlylevel"
  expect_error(.saf_run(fx, covs = c("cov1", "const"), data = d),
               "cannot adjust anything")
})

## Spec 54: an existing report is never silently replaced.

test_that("dmsa_report refuses to overwrite an existing report", {
  od <- tempfile("dmsa_ovw")
  dir.create(file.path(od, "figures"), recursive = TRUE)
  dir.create(file.path(od, "tables"),  recursive = TRUE)
  cat("x\n", file = file.path(od, "figures", "overview_y.png"))
  f <- structure(list(outdir = od, progress = FALSE, beep = FALSE),
                 class = "dmsa_frame")
  e <- expect_error(dmsa_report(f))
  expect_match(conditionMessage(e), "already contains a DMSA report")
  expect_match(conditionMessage(e), "No report was written")
})

test_that("the overwrite guard does not fire on a fresh dir or when opted in", {
  od <- tempfile("dmsa_fresh")
  f  <- structure(list(outdir = od, progress = FALSE, beep = FALSE),
                  class = "dmsa_frame")
  ## the stub frame fails later for unrelated reasons; what matters is that it
  ## is NOT stopped by the overwrite guard
  e <- tryCatch(dmsa_report(f), error = function(e) conditionMessage(e))
  expect_false(grepl("already contains a DMSA report", e))

  dir.create(file.path(od, "figures"), recursive = TRUE,
             showWarnings = FALSE)
  cat("x\n", file = file.path(od, "figures", "overview_y.png"))
  e2 <- tryCatch(dmsa_report(f, overwrite = TRUE),
                 error = function(e) conditionMessage(e))
  expect_false(grepl("already contains a DMSA report", e2))
})

test_that("all-singleton blocks that routing cannot rescue name the culprit", {
  ## a per-row ID nested in couples: the nested pair keeps both columns as
  ## blocks (correctly), but ID makes every group a singleton - the error
  ## must say which column ALONE would permute
  set.seed(1); fx <- .saf_fixture()
  d <- fx$data
  d$rowID <- paste0("r", seq_len(fx$n))          # unique per row
  e <- expect_error(
    dmsa_frame(d, methylation = fx$cols, direction_source = "bundled",
               outcomes = "y", covariates = "cov1",
               random_effects = c("rowID", "cID"), chip = FALSE, B = 19,
               plots = FALSE, tables = FALSE, summary = FALSE,
               progress = FALSE, beep = FALSE, outdir = tempfile("saf_x")))
  msg <- conditionMessage(e)
  expect_match(msg, "no permutable groups")
  expect_match(msg, "`cID` alone gives")
})

## ---- random_effects routing (PI, 2026-08-29): lme4-style intent ------------

.rt_run <- function(d, cols, ...) dmsa_frame(
  d, methylation = cols, direction_source = "bundled", outcomes = "y",
  covariates = "cov1", B = 19, plots = FALSE, tables = FALSE,
  summary = FALSE, progress = FALSE, beep = FALSE,
  outdir = tempfile("dmsa_rt"), ...)

.rt_fix <- function(n = 60, chip_of_row) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], 3)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2),
                  chip_T1 = factor(chip_of_row))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(d = d, cols = cols)
}

test_that("CROSSED factors route: blocks on cID, (1|chip) intercept", {
  set.seed(1)
  fx <- .rt_fix(chip_of_row = rep(1:2, 30))         # partners on DIFFERENT chips
  f <- .rt_run(fx$d, fx$cols, random_effects = c("cID", "chip_T1"))
  expect_identical(f$block_cols, "cID")
  expect_true(nzchar(f$chip_random))                # the (1|chip) intercept
  i <- grep("random_effects", f$corrections$field)
  expect_match(f$corrections$issue[i], "crossed/independent")
  expect_match(f$corrections$issue[i], "\\(1 \\| cID\\) \\+ \\(1 \\| chip_T1\\)")
})

test_that("NESTED factors are diagnosed but behaviour is unchanged", {
  set.seed(1)
  fx <- .rt_fix(chip_of_row = rep(1:6, each = 10))  # couples share a chip
  f <- .rt_run(fx$d, fx$cols, random_effects = c("cID", "chip_T1"),
               chip = FALSE)
  ## both columns stay in the block (their interaction IS the cID grouping)
  expect_setequal(f$block_cols, c("cID", "chip_T1"))
  i <- grep("random_effects", f$corrections$field)
  expect_match(f$corrections$issue[i], "nested in")
  expect_match(f$corrections$issue[i], "\\(1 \\| chip_T1/cID\\)")
})

test_that("routing refuses a conflicting explicit chip, and 3+ factors", {
  set.seed(1)
  fx <- .rt_fix(chip_of_row = rep(1:2, 30))
  fx$d$plate <- factor(rep(1:4, each = 15))
  expect_error(.rt_run(fx$d, fx$cols,
                       random_effects = c("cID", "chip_T1"), chip = "plate"),
               "already claims that slot")
  expect_error(.rt_run(fx$d, fx$cols,
                       random_effects = c("cID", "chip_T1", "plate")),
               "supports two")
  expect_error(.rt_run(fx$d, fx$cols,
                       random_effects = c("cID", "chip_T1"),
                       chip_effect = "none"),
               "Pick one")
})
