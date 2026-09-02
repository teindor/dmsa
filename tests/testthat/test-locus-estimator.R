## Spec 52: the locus figure draws the TESTED estimator.
##
## The panel used to draw a report-side OLS of each probe's raw M-value on
## the outcome, while the gene test pooled per-probe fits of the standardised
## (and, under a random intercept, RI-transformed) probe. Two estimators, one
## figure, no label. Now dmsa_triangulate() hands its own per-probe fits over
## (attr "probe_fits"), the panel draws those, and every number drawn is
## written to tables/locus_probes.csv.

## the engine's preprocessing of one probe: winsorise at 3 MADs, standardise
.le_prep <- function(y, winsor = 3) {
  md <- stats::median(y); s <- stats::mad(y)
  if (is.finite(s) && s > 0) y <- pmin(pmax(y, md - winsor * s), md + winsor * s)
  as.numeric(scale(y))
}
## ... and of the design: every non-intercept column standardised, so b is
## the slope per SD of the term on the standardised probe
.le_X <- function(formula, data) {
  X <- stats::model.matrix(formula, data)
  for (j in setdiff(colnames(X), "(Intercept)")) {
    s <- stats::sd(X[, j]); if (is.finite(s) && s > 0)
      X[, j] <- (X[, j] - mean(X[, j])) / s }
  X
}

.le_frame <- function(seed, od, ri = FALSE) {
  set.seed(seed); n <- 160; k <- 6
  pairs <- data.frame(cpg_id = sprintf("cg%07d", seq_len(k)),
                      target_gene = c(rep("NR3C1", 4), "FKBP5", "CRH"),
                      ## NR3C1's four probes all +1, so a chip shift common
                      ## to every probe survives into the ALIGNED score and
                      ## the random-intercept gamma is estimable from it
                      best_direction = c(rep(1, 4), -1, 1),
                      probability_plus1 = c(rep(.9, 4), .05, .9), usable = TRUE,
                      abstain_reason = NA_character_, best_evidence = "smr_high",
                      direction_tier = "S1", stringsAsFactors = FALSE)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2),
                  chip = factor(rep(sprintf("c%02d", 1:8), each = n / 8)))
  ## a chip effect, so a random intercept changes the per-probe fit
  ce <- stats::rnorm(8, sd = 0.8)[as.integer(d$chip)]
  for (j in seq_len(k))
    d[[pairs$cpg_id[j]]] <- stats::plogis(
      pairs$best_direction[j] * (if (j <= 4) 0.5 else 0) * scale(d$y)[, 1] +
        ce + stats::rnorm(n, 0, .7))
  old <- options(dmsa.pair_table = pairs); on.exit(options(old))
  suppressMessages(dmsa_frame(
    d, methylation = pairs$cpg_id, direction_source = "cpgdirection",
    outcome = "y", covariates = "cov1", blocks = "cID",
    chip = if (ri) TRUE else FALSE,
    B = 99, progress = FALSE, beep = FALSE, outdir = od))
}

test_that("spec 52: dmsa_triangulate hands over the per-probe fit it pooled", {
  set.seed(5); n <- 80; k <- 4
  M <- matrix(stats::rnorm(n * k), n, k, dimnames = list(NULL, paste0("cg", 1:k)))
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n))
  al <- dmsa_align(data.frame(cpg = colnames(M), d = rep(c(1, -1), 2),
                              p_plus = rep(c(.9, .1), 2)),
                   genes = rep("G", k), level = "gene")
  r <- dmsa_triangulate(M, d, rhs = c("y", "cov1"), term = "y",
                        units = rep("G", k), alignment = al, B = 49, seed = 1)
  pf <- attr(r, "probe_fits")
  expect_s3_class(pf, "data.frame")
  expect_identical(pf$column, colnames(M))
  ## the fit is of the WINSORISED, STANDARDISED probe on the term,
  ## covariate-adjusted - the engine's own preprocessing
  X <- .le_X(~ y + cov1, d)
  for (j in seq_len(k)) {
    fit <- stats::lm.fit(X, .le_prep(M[, j]))
    b <- unname(fit$coefficients["y"])
    se <- sqrt(sum(fit$residuals^2) / (n - ncol(X)) *
                 solve(crossprod(X))["y", "y"])
    expect_equal(pf$b[j], b, tolerance = 1e-10)
    expect_equal(pf$se[j], se, tolerance = 1e-10)
  }
  expect_equal(attr(pf, "dfr"), n - ncol(X))
  ## and the attribute survives the reordering of the output
  expect_true(!is.null(attr(r, "probe_fits")))
})

test_that("spec 52: the locus panel draws the engine's numbers and writes them", {
  od <- tempfile("dmsa_locus52")
  f <- .le_frame(seed = 3, od)
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  hits <- r$results[r$results$level == "gene" & r$results$selected, ]
  skip_if(!"NR3C1" %in% hits$unit, "NR3C1 not named at B = 99 on this seed")
  lp <- utils::read.csv(file.path(od, "tables", "locus_probes.csv"),
                        stringsAsFactors = FALSE)
  lp <- lp[lp$gene == "NR3C1", ]
  expect_equal(nrow(lp), 4L)
  expect_true(all(lp$estimator == "engine"))
  ## the drawn b/se equal the winsorised, standardised-probe fit the gene
  ## test pooled
  X <- .le_X(~ y + cov1, f$data)
  for (i in seq_len(nrow(lp))) {
    m <- .le_prep(f$M[, lp$column[i]])
    fit <- stats::lm.fit(X, m)
    expect_equal(lp$b[i], unname(fit$coefficients["y"]), tolerance = 1e-8,
                 info = lp$probe[i])
    expect_equal(lp$se[i], sqrt(sum(fit$residuals^2) / (nrow(X) - ncol(X)) *
                                solve(crossprod(X))["y", "y"]),
                 tolerance = 1e-8, info = lp$probe[i])
  }
  ## and NOT the old report-side OLS on raw M (different scale)
  old <- dmsa:::.rp_probe_fits(f, "y", lp$column)
  expect_false(isTRUE(all.equal(unname(old$b), lp$b, tolerance = 1e-6)))
  ## the nominal p in the table is the p of the drawn fit; the probe-level
  ## permutation p of the same probes agrees on which probes carry signal
  expect_equal(lp$p_nominal, 2 * stats::pt(-abs(lp$z), df = nrow(X) - ncol(X)),
               tolerance = 1e-10)
  pr <- utils::read.csv(file.path(od, "tables", "probes.csv"))
  m <- match(lp$probe, pr$unit)
  expect_false(anyNA(m))
  expect_equal(lp$p_nominal < .05, pr$p_coherence[m] < .05)
  ## summary.md points the reader at the table
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  expect_true(any(grepl("tables/locus_probes.csv", sm, fixed = TRUE)))
  expect_true(any(grepl("gene test's own per-probe estimates", sm, fixed = TRUE)))
})

test_that("spec 52: under a random intercept the panel follows the RI-transformed fit", {
  od <- tempfile("dmsa_locus52ri")
  f <- .le_frame(seed = 3, od, ri = TRUE)
  r <- suppressWarnings(suppressMessages(dmsa_report(f)))
  hits <- r$results[r$results$level == "gene" & r$results$selected, ]
  skip_if(!"NR3C1" %in% hits$unit, "NR3C1 not named at B = 99 on this seed")
  lp <- utils::read.csv(file.path(od, "tables", "locus_probes.csv"),
                        stringsAsFactors = FALSE)
  lp <- lp[lp$gene == "NR3C1", ]
  expect_true(all(lp$estimator == "engine"))
  ## the engine transformed X and Y by the RI gamma it estimated; a plain
  ## OLS of the preprocessed probe (no transform) must therefore DIFFER from
  ## what is drawn - the panel follows the test, not a re-fit
  X <- .le_X(stats::as.formula(paste("~", paste(
    dmsa:::.rp_rhs(f, "y"), collapse = "+"))), f$data)
  plain <- vapply(lp$column, function(cl)
    unname(stats::lm.fit(X, .le_prep(f$M[, cl]))$coefficients["y"]),
    numeric(1))
  expect_false(isTRUE(all.equal(unname(plain), lp$b, tolerance = 1e-6)))
  ## and it equals the fit the engine attached to the gene-level result
  al <- f$sets[[1]]$alignment
  tri <- dmsa_triangulate(f$M[, f$sets[[1]]$columns, drop = FALSE], f$data,
                          dmsa:::.rp_rhs(f, "y"), "y", f$sets[[1]]$map$gene,
                          al, block = f$block, B = 99, seed = f$seed,
                          ri_group = f$data[[f$chip_random]],
                          weighting = f$weighting %||% "combined",
                          w_floor = f$w_floor %||% 1.5,
                          correction = f$correction)
  pf <- attr(tri, "probe_fits")
  expect_equal(lp$b, pf$b[match(lp$column, pf$column)], tolerance = 1e-10)
  expect_equal(lp$se, pf$se[match(lp$column, pf$column)], tolerance = 1e-10)
})
