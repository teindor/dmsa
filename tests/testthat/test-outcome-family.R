## Spec 50: the outcome family belongs to the OUTCOME, not to the frame.

.otf_fixture <- function(n = 60, k = 3) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], k)
  d <- data.frame(score = stats::rnorm(n),
                  depressed = rep(0:1, length.out = n),
                  cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols, n = n)
}

.otf_run <- function(fx, outcome, ...)
  dmsa_frame(fx$data, methylation = fx$cols, direction_source = "bundled",
             outcome = outcome,
             covariates = "cov1", blocks = "cID", chip = FALSE,
             B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
             progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_otf"), ...)

test_that("a single-outcome frame still stores the bare family string", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, "score")
  expect_identical(f$outcome_type, "gaussian")
  expect_null(names(f$outcome_type))
})

test_that("the binary autofix still fires for a single outcome", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, "depressed")
  expect_identical(f$outcome_type, "logistic")
})

test_that("a binary outcome no longer flips the family of a continuous one", {
  ## THE BUG: outcome_type was one scalar for the whole frame, and the autofix
  ## overwrote it. "depressed" turned "score" logistic.
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, c("score", "depressed"))
  expect_equal(length(f$outcome_type), 2L)
  expect_identical(f$outcome_type[["score"]], "gaussian")
  expect_identical(f$outcome_type[["depressed"]], "logistic")
})

test_that("order of the outcomes does not change either family", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, c("depressed", "score"))
  expect_identical(f$outcome_type[["score"]], "gaussian")
  expect_identical(f$outcome_type[["depressed"]], "logistic")
})

test_that("the correction note is filed against the outcome that moved", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, c("score", "depressed"))
  cr <- f$corrections
  expect_true(any(cr$field == "depressed"))
  expect_false(any(cr$field == "score" &
                   grepl("two-level", paste(cr[[2]], cr[[3]]))))
})

test_that("a named outcome_type is honoured per outcome", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, c("score", "depressed"),
                outcome_type = c(score = "gaussian", depressed = "logistic"))
  expect_identical(f$outcome_type[["score"]], "gaussian")
  expect_identical(f$outcome_type[["depressed"]], "logistic")
})

test_that("one unnamed value still applies to every outcome", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, c("score", "depressed"), outcome_type = "logistic")
  expect_true(all(f$outcome_type == "logistic"))
})

test_that("an unnamed multi-value outcome_type is refused", {
  set.seed(1); fx <- .otf_fixture()
  e <- expect_error(.otf_run(fx, c("score", "depressed"),
                             outcome_type = c("gaussian", "logistic")))
  expect_match(conditionMessage(e), "does not name them")
  expect_match(conditionMessage(e), "Positional matching is not accepted")
})

test_that("a named outcome_type must name every outcome", {
  set.seed(1); fx <- .otf_fixture()
  e <- expect_error(.otf_run(fx, c("score", "depressed"),
                             outcome_type = c(depressed = "logistic")))
  expect_match(conditionMessage(e), "does not declare a family")
  expect_match(conditionMessage(e), "score")
})

test_that("outcome_type may not name a non-outcome column", {
  set.seed(1); fx <- .otf_fixture()
  e <- expect_error(.otf_run(fx, "score",
                             outcome_type = c(score = "gaussian",
                                              cov1 = "logistic")))
  expect_match(conditionMessage(e), "not outcomes")
})

test_that("an unrecognised family is refused", {
  set.seed(1); fx <- .otf_fixture()
  e <- expect_error(.otf_run(fx, "score", outcome_type = "poisson"))
  expect_match(conditionMessage(e), "not recognised")
})

test_that("multinomial is still refused, named or not", {
  set.seed(1); fx <- .otf_fixture()
  expect_error(.otf_run(fx, "score", outcome_type = "multinomial"),
               "not supported")
  expect_error(.otf_run(fx, c("score", "depressed"),
                        outcome_type = c(score = "gaussian",
                                         depressed = "multinomial")),
               "not supported")
})

test_that("autofix = FALSE names the argument and value that fixes it", {
  set.seed(1); fx <- .otf_fixture()
  e <- expect_error(.otf_run(fx, c("score", "depressed"), autofix = FALSE))
  msg <- conditionMessage(e)
  expect_match(msg, "depressed")
  expect_match(msg, "logistic")
  expect_match(msg, "score = \"gaussian\"", fixed = TRUE)
})

test_that("the printer shows the family beside each outcome", {
  set.seed(1); fx <- .otf_fixture()
  f <- .otf_run(fx, c("score", "depressed"))
  out <- paste(utils::capture.output(print(f)), collapse = "\n")
  expect_match(out, "score \\[gaussian\\]")
  expect_match(out, "depressed \\[logistic\\]")
})

## ---- binary labels + storage types (PI, 2026-08-29) -----------------------
## "taking pills [1] was associated ... relative to not taking pills [0]":
## outcome_levels accepts a named list (per-outcome labels), binary detection
## works for factor/character/logical storage, and the direction sentence
## leads with the high level, carries the coded values in brackets, and is
## capitalized.

test_that("two-level outcomes are binary whatever the storage type", {
  base <- data.frame(num = c(0, 1, 0, 1), fac = factor(c("no", "yes", "no", "yes")),
                     chr = c("a", "b", "a", "b"), lgl = c(TRUE, FALSE, TRUE, FALSE),
                     con = c(1.2, 3.4, 2.2, 8.1))
  for (cc in c("num", "fac", "chr", "lgl"))
    expect_identical(.frame_outcome_kind(base, cc, list(), "parent")$kind,
                     "binary", info = cc)
  expect_identical(.frame_outcome_kind(base, "con", list(), "parent")$kind,
                   "continuous")
})

test_that("the binary direction sentence uses declared labels + coded values", {
  fr <- list(outcome = c("pills", "anx"),
             outcome_kind = list(pills = list(kind = "binary", lo = 0, hi = 1)),
             outcome_levels = list(pills = c("not taking pills",
                                             "taking pills")),
             labels = NULL)
  s <- .rp_dir_sentence(fr, "pills", +1, "GREB1")
  expect_match(s, "^Taking pills \\[1\\] was associated with methylation consistent with a HIGHER expression tone of GREB1, relative to not taking pills \\[0\\]")
  ## without labels: the coded-value fallback, no double brackets
  fr$outcome_levels <- NULL
  s2 <- .rp_dir_sentence(fr, "pills", -1, "GREB1")
  expect_match(s2, "Pills = 1 was associated")
  expect_false(grepl("\\[", s2))
})

test_that("outcome_levels as a list is validated against the outcomes", {
  d <- data.frame(y = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2),
                  cg01 = stats::plogis(rnorm(40)))
  map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                    probe = "cg01", column = "cg01",
                    best_direction = -1, p_plus = .1)
  expect_error(dmsa_frame(d, map = map, outcomes = "y", covariates = "cov1",
                          blocks = "cID", B = 19, plots = FALSE,
                          tables = FALSE, summary = FALSE, progress = FALSE,
                          beep = FALSE, outdir = tempfile(),
                          outcome_levels = list(WRONG = c("a", "b"))),
               "named by outcome")
})

test_that("gene_models accepts auto/TRUE/FALSE/table and refuses the rest", {
  d <- data.frame(y = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2),
                  cg01 = stats::plogis(rnorm(40)))
  map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                    probe = "cg01", column = "cg01",
                    best_direction = -1, p_plus = .1)
  expect_error(dmsa_frame(d, map = map, outcomes = "y", covariates = "cov1",
                          blocks = "cID", B = 19, plots = FALSE,
                          tables = FALSE, summary = FALSE, progress = FALSE,
                          beep = FALSE, outdir = tempfile(),
                          gene_models = "bogus"),
               "auto")
  f <- dmsa_frame(d, map = map, outcomes = "y", covariates = "cov1",
                  blocks = "cID", B = 19, plots = FALSE,
                  tables = FALSE, summary = FALSE, progress = FALSE,
                  beep = FALSE, outdir = tempfile())   # default = "auto"
  expect_identical(f$gene_models, "auto")
})

test_that("predictor_levels mirrors outcome_levels in the user's vocabulary", {
  set.seed(6)
  map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                    probe = "cg01", column = "cg01",
                    best_direction = -1, p_plus = .1)
  d <- data.frame(pills = rep(c(0, 1), 20), cov1 = rnorm(40),
                  cID = rep(1:20, each = 2),
                  cg01 = stats::plogis(rnorm(40)))
  ## declared with predictors= -> predictor_levels works and lands on the frame
  f <- dmsa_frame(d, map = map, predictors = "pills", covariates = "cov1",
                  blocks = "cID", B = 19, plots = FALSE,
                  tables = FALSE, summary = FALSE, progress = FALSE,
                  beep = FALSE, outdir = tempfile(),
                  predictor_levels = list(pills = c("no pills", "pills")))
  expect_identical(f$outcome_levels$pills, c("no pills", "pills"))
  ## both spellings at once -> refused in plain words
  expect_error(dmsa_frame(d, map = map, predictors = "pills",
                          covariates = "cov1", blocks = "cID",
                          B = 19, plots = FALSE, tables = FALSE,
                          summary = FALSE, progress = FALSE, beep = FALSE,
                          outdir = tempfile(),
                          predictor_levels = c("a", "b"),
                          outcome_levels = c("a", "b")),
               "ONCE")
  ## predictor_levels with outcomes= -> pointed to outcome_levels
  expect_error(dmsa_frame(d, map = map, outcomes = "pills",
                          covariates = "cov1", blocks = "cID",
                          B = 19, plots = FALSE, tables = FALSE,
                          summary = FALSE, progress = FALSE, beep = FALSE,
                          outdir = tempfile(),
                          predictor_levels = c("a", "b")),
               "outcome_levels")
})

## ---- label spellings + value-matched level labels (PI, 2026-08-29) --------
## predictor_labels names the `predictors = ` column (it used to silently
## label the covariates - the PI's own call errored with "1 label(s) for 6
## column(s)"); covariate_labels is the covariates' own argument; named
## level vectors are matched BY VALUE so a 1/2 coding cannot flip.

.lab_fix <- function() {
  map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                    probe = "cg01", column = "cg01",
                    best_direction = -1, p_plus = .1)
  d <- data.frame(pills = rep(c(0, 1), 20), sex = rep(c(1, 2), each = 20),
                  cov1 = rnorm(40), cID = rep(1:20, each = 2),
                  cg01 = stats::plogis(rnorm(40)))
  list(map = map, d = d)
}

test_that("predictor_labels names the predictors= column, one label suffices", {
  set.seed(8); fx <- .lab_fix()
  f <- dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                  covariates = "cov1", blocks = "cID", B = 19,
                  plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile(),
                  predictor_labels = "Used birth-control pills when met")
  expect_identical(unname(f$labels[["pills"]]),
                   "Used birth-control pills when met")
  ## with outcomes= the same argument errors, pointing to the right ones
  expect_error(dmsa_frame(fx$d, map = fx$map, outcomes = "pills",
                          covariates = "cov1", blocks = "cID",
                          B = 19, plots = FALSE, tables = FALSE,
                          summary = FALSE, progress = FALSE, beep = FALSE,
                          outdir = tempfile(),
                          predictor_labels = "x"),
               "covariate_labels")
})

test_that("covariate_labels names the covariates", {
  set.seed(8); fx <- .lab_fix()
  f <- dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                  covariates = "cov1", blocks = "cID", B = 19,
                  plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile(),
                  covariate_labels = c(cov1 = "Control 1"))
  expect_identical(unname(f$labels[["cov1"]]), "Control 1")
})

test_that("a length-2 mod_label is refused and pointed at mod_levels", {
  set.seed(8); fx <- .lab_fix()
  expect_error(dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                          covariates = "cov1", blocks = "cID",
                          moderation = TRUE, mod = "sex",
                          B = 19, plots = FALSE, tables = FALSE,
                          summary = FALSE, progress = FALSE, beep = FALSE,
                          outdir = tempfile(),
                          mod_label = c("Husband", "Wife")),
               "mod_levels")
})

test_that("named level labels are matched by VALUE and cannot flip", {
  set.seed(8); fx <- .lab_fix()
  ## deliberately out of order: "2" first - resolution is by value
  f <- dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                  covariates = "cov1", blocks = "cID",
                  moderation = TRUE, mod = "sex",
                  B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile(),
                  mod_levels = c("2" = "Wife", "1" = "Husband"),
                  predictor_levels = list(pills = c("1" = "taking pills",
                                                    "0" = "not taking pills")))
  expect_identical(f$mod_levels, c("Husband", "Wife"))    # low, high
  expect_identical(f$outcome_levels$pills,
                   c("not taking pills", "taking pills"))
  ## a name that is not one of the variable's levels errors, naming them
  expect_error(dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                          covariates = "cov1", blocks = "cID",
                          moderation = TRUE, mod = "sex",
                          B = 19, plots = FALSE, tables = FALSE,
                          summary = FALSE, progress = FALSE, beep = FALSE,
                          outdir = tempfile(),
                          mod_levels = c("0" = "Husband", "2" = "Wife")),
               "1 and 2")
})

test_that("a column given as covariate AND moderator is flagged and deduped", {
  set.seed(8); fx <- .lab_fix()
  ## moderation = TRUE: dropped from covariates, said out loud
  expect_message(
    f <- dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                    covariates = c("sex", "cov1"), blocks = "cID",
                    moderation = TRUE, mod = "sex",
                    B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                    progress = FALSE, beep = FALSE, outdir = tempfile()),
    "both as a covariate and as a moderator")
  expect_false("sex" %in% f$covariates)
  expect_true(any(grepl("covariate AND moderator", f$corrections$issue)))
  ## moderation = FALSE (the PI's actual call): live warning that NO
  ## moderation runs and sex stays a plain covariate
  expect_message(
    f2 <- dmsa_frame(fx$d, map = fx$map, predictors = "pills",
                     covariates = c("sex", "cov1"), blocks = "cID",
                     mod = "sex",
                     B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                     progress = FALSE, beep = FALSE, outdir = tempfile()),
    "moderation = FALSE")
  expect_true("sex" %in% f2$covariates)
})
