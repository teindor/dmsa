## Spec 48: missing methylation is a declared policy, never a silent drop.
## Spec 49: the methylation scale is declared, not guessed.

.mth_fixture <- function(n = 60, k = 4) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], k)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols, n = n)
}

.mth_run <- function(fx, data = fx$data, ...)
  dmsa_frame(data, methylation = fx$cols, direction_source = "bundled",
             outcome = "y", covariates = "cov1",
             random_effects = "cID", chip = FALSE, B = 19, plots = FALSE,
             tables = FALSE, summary = FALSE, progress = FALSE, beep = FALSE,
             outdir = tempfile("dmsa_mth"), ...)

## ---- spec 48 --------------------------------------------------------------

test_that("clean methylation is unaffected by the policy argument", {
  set.seed(1); fx <- .mth_fixture()
  f <- .mth_run(fx)
  expect_s3_class(f, "dmsa_frame")
  expect_equal(nrow(f$data), fx$n)
})

test_that("missing methylation errors by default and names the policies", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data; d[[fx$cols[1]]][3] <- NA
  e <- expect_error(.mth_run(fx, data = d))
  msg <- conditionMessage(e)
  expect_match(msg, "missing or non-finite")
  expect_match(msg, "drop_probes")
  expect_match(msg, "common_complete_rows")
  expect_match(msg, "No analysis was run")
})

test_that("non-finite values that are not NA are caught too", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data; d[[fx$cols[2]]][5] <- Inf
  expect_error(.mth_run(fx, data = d), "missing or non-finite")
})

test_that("drop_probes drops the affected probe and records it", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data; d[[fx$cols[1]]][3] <- NA
  f <- .mth_run(fx, data = d, missing_methylation = "drop_probes")
  expect_s3_class(f, "dmsa_frame")
  expect_equal(nrow(f$data), fx$n)             # every sample retained
  expect_false(fx$cols[1] %in% f$map$column)   # the probe is gone
  expect_true(fx$cols[2] %in% f$map$column)
  expect_true(any(grepl("drop_probes", unlist(f$corrections))))
})

test_that("common_complete_rows drops the sample and keeps every probe", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data; d[[fx$cols[1]]][3] <- NA
  f <- .mth_run(fx, data = d, missing_methylation = "common_complete_rows")
  expect_s3_class(f, "dmsa_frame")
  expect_equal(nrow(f$data), fx$n - 1L)        # one sample dropped
  expect_true(all(fx$cols %in% f$map$column))  # every probe retained
  expect_true(any(grepl("common_complete_rows", unlist(f$corrections))))
})

test_that("common_complete_rows refuses to run on too few complete samples", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data
  d[[fx$cols[1]]][seq_len(fx$n - 10L)] <- NA
  e <- expect_error(.mth_run(fx, data = d,
                             missing_methylation = "common_complete_rows"))
  expect_match(conditionMessage(e), "too few to analyse")
  expect_match(conditionMessage(e), "drop_probes")
})

test_that("drop_probes still errors when nothing survives", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data
  for (cl in fx$cols) d[[cl]][7] <- NA
  e <- expect_error(.mth_run(fx, data = d, missing_methylation = "drop_probes"))
  expect_match(conditionMessage(e), "no probe is left to analyse")
})

test_that("an unknown policy is rejected", {
  set.seed(1); fx <- .mth_fixture()
  expect_error(.mth_run(fx, missing_methylation = "impute"))
})

## ---- spec 49 --------------------------------------------------------------

test_that("declared beta values outside [0, 1] are an error, not a clamp", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data; d[[fx$cols[1]]] <- d[[fx$cols[1]]] * 8 - 4   # M-value-like
  e <- expect_error(.mth_run(fx, data = d, methylation_scale = "beta"))
  expect_match(conditionMessage(e), "outside \\[0, 1\\]")
  expect_match(conditionMessage(e), "No analysis was run")
})

test_that("declared beta on real beta values runs and converts", {
  set.seed(1); fx <- .mth_fixture()
  f <- .mth_run(fx, methylation_scale = "beta")
  expect_s3_class(f, "dmsa_frame")
  expect_true(any(grepl("M-values", unlist(f$corrections))))
})

test_that("methylation_scale = 'M' converts nothing", {
  set.seed(1); fx <- .mth_fixture()
  d <- fx$data
  for (cl in fx$cols) d[[cl]] <- stats::rnorm(fx$n, 0, 2)    # honest M-values
  f <- .mth_run(fx, data = d, methylation_scale = "M")
  expect_s3_class(f, "dmsa_frame")
  expect_false(any(grepl("converted to M-values", unlist(f$corrections))))
})

test_that("'M' on beta-looking input is taken at its word", {
  ## auto would convert; declaring M must not. The declaration is the contract.
  set.seed(1); fx <- .mth_fixture()
  fa <- .mth_run(fx, methylation_scale = "auto")
  fm <- .mth_run(fx, methylation_scale = "M")
  expect_true(any(grepl("converted to M-values", unlist(fa$corrections))))
  expect_false(any(grepl("converted to M-values", unlist(fm$corrections))))
})

test_that("an unknown scale is rejected", {
  set.seed(1); fx <- .mth_fixture()
  expect_error(.mth_run(fx, methylation_scale = "beta_values"))
})
