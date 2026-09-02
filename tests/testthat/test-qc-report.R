## Spec 53: report failures are RECORDED, never swallowed.
##
## Thirteen places in the report used tryCatch(..., error = function(e) NULL):
## a component that failed simply vanished, and a reader could not tell
## "absent because not applicable" from "absent because it broke". Every one
## now goes through .rp_try(), which records (component, where, fallback,
## message) in a run-level ledger written to tables/qc_report.csv, listed in
## summary.md under "Report QC", counted once on the console, and carried on
## the report object as $qc. The analysis itself never stops for it.

.qc_frame <- function(od, seed = 21) {
  set.seed(seed)
  pairs <- data.frame(cpg_id = sprintf("cg%07d", 1:3),
                      target_gene = c("NR3C1", "FKBP5", "CRH"),
                      best_direction = c(-1, 1, 1),
                      probability_plus1 = c(.05, .9, .9), usable = TRUE,
                      abstain_reason = NA_character_, best_evidence = "smr_high",
                      direction_tier = "S1", stringsAsFactors = FALSE)
  n <- 80
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(1:40, each = 2))
  for (cl in pairs$cpg_id) d[[cl]] <- stats::plogis(stats::rnorm(n))
  d$cg0000001 <- stats::plogis(-1.4 * scale(d$y)[, 1] + stats::rnorm(n, 0, .4))
  old <- options(dmsa.pair_table = pairs); on.exit(options(old))
  suppressMessages(dmsa_frame(
    d, methylation = pairs$cpg_id, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = FALSE, B = 99, progress = FALSE, beep = FALSE, outdir = od))
}

test_that("spec 53: a clean run writes an EMPTY ledger and says so", {
  od <- tempfile("dmsa_qc_clean")
  f <- .qc_frame(od)
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  qf <- file.path(od, "tables", "qc_report.csv")
  expect_true(file.exists(qf))
  q <- utils::read.csv(qf, stringsAsFactors = FALSE)
  expect_identical(names(q), c("component", "severity", "outcome", "level",
                               "unit", "fallback", "message"))
  expect_equal(nrow(q), 0L)
  expect_s3_class(r$qc, "data.frame"); expect_equal(nrow(r$qc), 0L)
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  expect_true(any(grepl("^## Report QC$", sm)))
  expect_true(any(grepl("No report component failed or fell back", sm)))
  expect_true(any(grepl("tables/qc_report", r$tables)))
})

test_that("spec 53: a failing component is recorded, the report completes, the reader is told", {
  od <- tempfile("dmsa_qc_fail")
  f <- .qc_frame(od)
  ## make the manifest-coordinate lookup blow up: the locus panel then falls
  ## back to evenly spaced probes, and the failure must land in the ledger
  testthat::local_mocked_bindings(
    dmsa_probe_coords = function(...) stop("simulated manifest failure"),
    .package = "dmsa")
  msgs <- character()
  r <- withCallingHandlers(
    suppressWarnings(dmsa_report(f)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                            invokeRestart("muffleMessage") })
  hits <- r$results[r$results$level == "gene" & r$results$selected, ]
  skip_if(!nrow(hits), "no gene named at B = 99 on this seed")
  ## the analysis is intact and the panel was still drawn
  expect_true(nrow(r$results) > 0)
  expect_true(any(grepl("^locus_", basename(r$figures))))
  ## the ledger
  q <- utils::read.csv(file.path(od, "tables", "qc_report.csv"),
                       stringsAsFactors = FALSE)
  expect_true(nrow(q) >= 1L)
  i <- which(q$component == "locus: manifest coordinates")
  expect_true(length(i) >= 1L)
  expect_identical(q$severity[i[1]], "failed")
  expect_match(q$message[i[1]], "simulated manifest failure")
  expect_identical(q$unit[i[1]], hits$unit[1])
  expect_match(q$fallback[i[1]], "evenly spaced")
  ## the object carries it, summary.md lists it, the console counted it once
  expect_equal(nrow(r$qc), nrow(q))
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  expect_true(any(grepl("report component\\(s\\) failed or fell back", sm)))
  expect_true(any(grepl("locus: manifest coordinates", sm, fixed = TRUE)))
  expect_true(any(grepl("simulated manifest failure", sm, fixed = TRUE)))
  expect_equal(sum(grepl("report component\\(s\\) failed or fell back - see tables/qc_report.csv", msgs)), 1L)
})

test_that("spec 53: the ledger is reset between runs", {
  od <- tempfile("dmsa_qc_reset")
  f <- .qc_frame(od)
  dmsa:::.rp_qc_note("stale", "left over from a previous run")
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  expect_false(any(r$qc$component == "stale"))
})

test_that("spec 53: .rp_try records and returns NULL; success passes through", {
  dmsa:::.rp_qc_reset()
  expect_null(dmsa:::.rp_try("x", stop("nope"), unit = "U", fallback = "fb"))
  expect_identical(dmsa:::.rp_try("y", 42), 42)
  q <- dmsa:::.rp_qc_rows()
  expect_equal(nrow(q), 1L)
  expect_identical(q$component, "x"); expect_identical(q$unit, "U")
  expect_identical(q$fallback, "fb"); expect_match(q$message, "nope")
  expect_identical(q$severity, "failed")
})
