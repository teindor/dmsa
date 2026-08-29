## Spec 46 + 65: a methylation matrix is never assumed to be row-aligned with
## the phenotype table when identifiers exist to check it against.

.al_fixture <- function(n = 60, k = 3) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], k)
  ids <- sprintf("S%02d", seq_len(n))
  d <- data.frame(sid = ids, y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2), stringsAsFactors = FALSE)
  M <- matrix(stats::plogis(stats::rnorm(n * length(cols))), n,
              dimnames = list(ids, cols))
  list(data = d, M = M, n = n)
}

.al_run <- function(fx, M, ...)
  dmsa_frame(fx$data, methylation = M, direction_source = "bundled",
             outcome = "y", covariates = "cov1",
             random_effects = "cID", chip = FALSE, B = 19, plots = FALSE,
             tables = FALSE, summary = FALSE, progress = FALSE, beep = FALSE,
             outdir = tempfile("dmsa_al"), ...)

test_that("shuffled rows with matching IDs are REORDERED, not analysed as-is", {
  set.seed(1); fx <- .al_fixture()
  Msh <- fx$M[sample(fx$n), , drop = FALSE]
  f <- .al_run(fx, Msh, sample_id = "sid")
  expect_s3_class(f, "dmsa_frame")
  expect_true(any(grepl("different order", f$corrections$issue)))
})

test_that("mismatched identifiers are a hard error", {
  set.seed(1); fx <- .al_fixture()
  Mbad <- fx$M; rownames(Mbad)[1:3] <- c("ZZ1", "ZZ2", "ZZ3")
  e <- expect_error(.al_run(fx, Mbad, sample_id = "sid"))
  expect_match(conditionMessage(e), "do not match")
  expect_match(conditionMessage(e), "No analysis was run")
})

test_that("sample_id without matrix rownames is a hard error", {
  set.seed(1); fx <- .al_fixture()
  M <- fx$M; rownames(M) <- NULL
  expect_error(.al_run(fx, M, sample_id = "sid"), "no rownames")
})

test_that("duplicated identifiers are refused on either side", {
  set.seed(1); fx <- .al_fixture()
  d2 <- fx$data; d2$sid[2] <- d2$sid[1]
  fx2 <- fx; fx2$data <- d2
  expect_error(.al_run(fx2, fx$M, sample_id = "sid"), "duplicated identifier")
})

test_that("a row-count mismatch is always an error", {
  set.seed(1); fx <- .al_fixture()
  expect_error(.al_run(fx, fx$M[-1, , drop = FALSE], sample_id = "sid"),
               "must describe the same samples")
})

test_that("no identifiers anywhere: positional, but RECORDED not silent", {
  set.seed(1); fx <- .al_fixture()
  M <- fx$M; rownames(M) <- NULL
  f <- .al_run(fx, M)
  a <- f$corrections[f$corrections$field == "sample alignment", ]
  expect_equal(nrow(a), 1L)
  expect_match(a$action, "POSITIONALLY")
})

test_that("rownames alone trigger announced auto-detection and reordering", {
  set.seed(1); fx <- .al_fixture()
  M2 <- fx$M[c(2:fx$n, 1), , drop = FALSE]
  expect_message(f <- .al_run(fx, M2), "auto-detected")
  expect_true(any(grepl("different order", f$corrections$issue)))
})
