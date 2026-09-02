## Spec 51: the model orientation is declared, consistent, and printed.

.mor_fixture <- function(n = 60, k = 3) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], k)
  d <- data.frame(score = stats::rnorm(n), cov1 = stats::rnorm(n),
                  w = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols, n = n)
}

.mor_run <- function(fx, ...)
  dmsa_frame(fx$data, methylation = fx$cols, direction_source = "bundled",
             outcome = "score",
             covariates = "cov1", blocks = "cID", chip = FALSE,
             B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
             progress = FALSE, beep = FALSE, outdir = tempfile("dmsa_mor"), ...)

test_that("the main path is declared as methylation-on-outcome", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx)
  expect_match(f$model_orientation$main, "^methylation ~ score")
})

test_that("a linear frame with no moderation declares no secondary model", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx)
  expect_true(is.na(f$model_orientation$shape))
  expect_true(is.na(f$model_orientation$moderation))
})

test_that("type != linear declares the shape scan, outcome as response", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx, type = "non-linear")
  expect_match(f$model_orientation$shape, "^score ~ tone score")
})

test_that("moderation under frame_role = predictor faces the outcome", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx, moderation = TRUE, mod = "w")
  expect_match(f$model_orientation$moderation, "^score ~ tone score x moderator")
})

test_that("moderation under frame_role = outcome faces the tone score", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx, moderation = TRUE, mod = "w", frame_role = "outcome")
  expect_match(f$model_orientation$moderation, "^tone score ~ score")
})

test_that("frame_role = outcome with a shape scan is refused, not mixed", {
  ## THE BUG: moderation honoured frame_role, the shape scan never did, so one
  ## report carried a moderated model facing one way and a curvature test
  ## facing the other.
  set.seed(1); fx <- .mor_fixture()
  e <- expect_error(.mor_run(fx, type = "non-linear", frame_role = "outcome",
                             moderation = TRUE, mod = "w"))
  msg <- conditionMessage(e)
  expect_match(msg, "cannot be honoured together with type")
  expect_match(msg, "frame_role = \"predictor\"", fixed = TRUE)
  expect_match(msg, "type = \"linear\"", fixed = TRUE)
  expect_match(msg, "No analysis was run")
})

test_that("frame_role = outcome without moderation says it does nothing", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx, frame_role = "outcome")
  i <- which(f$corrections$field == "frame_role")
  expect_gte(length(i), 1L)
  expect_match(paste(unlist(f$corrections[i, ]), collapse = " "),
               "no effect on this frame")
})

test_that("frame_role = predictor without moderation files no such note", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx)
  expect_false(any(f$corrections$field == "frame_role"))
})

test_that("the printer names the model each active path fits", {
  set.seed(1); fx <- .mor_fixture()
  f <- .mor_run(fx, moderation = TRUE, mod = "w")
  out <- paste(utils::capture.output(print(f)), collapse = "\n")
  expect_match(out, "main\\s+methylation ~ score")
  expect_match(out, "moder\\.\\s+score ~ tone score x moderator")
  expect_match(out, "frame_role = predictor")
})

test_that("an unknown frame_role is refused", {
  set.seed(1); fx <- .mor_fixture()
  expect_error(.mor_run(fx, frame_role = "mediator"))
})

## ---- PI request (2026-08-29): predictors= / outcomes= spellings ------------

test_that("predictors= names the columns AND sets the orientation", {
  set.seed(1); fx <- .mor_fixture()
  f <- dmsa_frame(fx$data, methylation = fx$cols,
                  direction_source = "bundled",
                  predictors = "score",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile("mor_p"))
  expect_identical(f$frame_role, "outcome")
  expect_identical(f$outcome, "score")
  ## and it is EQUIVALENT to the old spelling
  set.seed(1)
  g <- .mor_run(fx, frame_role = "outcome")
  expect_identical(f$outcome, g$outcome)
  expect_identical(f$frame_role, g$frame_role)
})

test_that("outcomes= is the interface; legacy outcome= still works", {
  set.seed(1); fx <- .mor_fixture()
  f <- dmsa_frame(fx$data, methylation = fx$cols,
                  direction_source = "bundled",
                  outcomes = "score",
                  covariates = "cov1", blocks = "cID", chip = FALSE,
                  B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                  progress = FALSE, beep = FALSE, outdir = tempfile("mor_o"))
  expect_identical(f$outcome, "score")
  expect_identical(f$frame_role, "predictor")
  ## `outcome` is NOT in the signature any more...
  expect_false("outcome" %in% names(formals(dmsa_frame)))
  ## ...but every existing script that passes it keeps working, silently
  set.seed(1)
  g <- .mor_run(fx)                       # the helper still uses outcome=
  expect_identical(g$outcome, "score")
  ## and a genuinely unknown argument is still a loud error, not a typo sink
  expect_error(
    dmsa_frame(fx$data, methylation = fx$cols, outcomes = "score",
               covariate = "cov1"),
    "unknown argument")
})

test_that("contradictions and duplicates are refused with plain words", {
  set.seed(1); fx <- .mor_fixture()
  base <- function(...) dmsa_frame(fx$data, methylation = fx$cols,
    direction_source = "bundled", covariates = "cov1",
    blocks = "cID", chip = FALSE, B = 19, plots = FALSE,
    tables = FALSE, summary = FALSE, progress = FALSE, beep = FALSE,
    outdir = tempfile("mor_x"), ...)
  expect_error(base(predictors = "score", frame_role = "predictor"),
               "already sets the orientation")
  expect_error(base(outcome = "score", predictors = "score"),
               "Pick one")
  expect_error(base(outcome = "score", outcomes = "score"),
               "Pick one")
  expect_error(base(), "outcomes =|predictors =")
})
