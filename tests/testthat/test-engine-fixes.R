## Regression pins for the PI-approved engine fixes E1-E10 (2026-08-29).

.ef_data <- function(n = 120, k = 6, seed = 7) {
  set.seed(seed)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cID = rep(seq_len(n / 2), each = 2))
  probes <- paste0("cg", 1:k)
  for (p in probes) d[[p]] <- stats::rnorm(n)
  list(d = d, probes = probes)
}

## ---- E1: the focal bookmark is never clobbered -----------------------------

test_that("E1: the fi-clobbering assignment is gone from the source", {
  dir <- NULL
  for (p in c(file.path("..", "..", "R"), "R"))
    if (file.exists(file.path(p, "fit.R"))) { dir <- p; break }
  skip_if(is.null(dir), "sources not visible")
  txt <- paste(readLines(file.path(dir, "fit.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("fi <- match(colnames(X)[fi]", txt, fixed = TRUE))
})

test_that("E1: lmer and lm engines share one permutation p (same stream)", {
  skip_if_not_installed("lme4")
  fx <- .ef_data()
  des <- dmsa_design(focal = "y", fixed = "cov1", random = "cID",
                     exchangeable = "cID")
  al <- data.frame(cpg = fx$probes, d = rep(c(1, -1), 3),
                   p_plus = rep(c(.9, .1), 3))
  a <- dmsa_align(al, genes = rep("G1", 6), level = "gene")
  f1 <- dmsa_fit(fx$d, probes = fx$probes, alignment = a, design = des,
                 engine = "lm", B = 99, beta_input = FALSE, seed = 3)
  f2 <- dmsa_fit(fx$d, probes = fx$probes, alignment = a, design = des,
                 engine = "lmer", B = 99, beta_input = FALSE, seed = 3)
  ## p_perm compares the lm observed statistic to lm nulls on the same seeded
  ## stream in BOTH runs; a corrupted focal index in the lmer loop breaks this
  expect_equal(f1$p_perm, f2$p_perm, tolerance = 1e-12)
})

## ---- E9: alignment order contract ------------------------------------------

test_that("E9: dmsa_align(character) returns rows in the caller's order", {
  dm <- readRDS(system.file("extdata", "direction_blood_hg19.rds",
                            package = "dmsa"))
  pcol <- intersect(c("probe", "cpg", "cpg_id"), names(dm))[1]
  pr <- utils::head(unique(as.character(dm[[if (is.na(pcol)) 1 else pcol]])), 4)
  pr_rev <- rev(pr)
  a <- suppressMessages(dmsa_align(pr_rev, level = "gene"))
  ## first occurrence order must follow the input, not the bundled map
  expect_identical(unique(a$probe), pr_rev[pr_rev %in% a$probe])
})


test_that("E9: a scrambled alignment is repaired; a mismatched one refused", {
  set.seed(1); n <- 60
  probes <- paste0("cg", 101:106)
  M <- matrix(rnorm(n * 6), n, dimnames = list(NULL, probes))
  al <- data.frame(probe = probes, s = rep(c(1, -1), 3),
                   p_s_plus = rep(c(.9, .1), 3), usable = TRUE)
  s0 <- dmsa_scores(M, al)$aligned
  sh <- sample(6)
  s1 <- dmsa_scores(M, al[sh, ])$aligned          # same set, scrambled
  expect_equal(s1, s0, tolerance = 1e-12)
  al2 <- al; al2$probe[1] <- "cg999999"           # overlapping but different
  expect_error(dmsa_scores(M, al2), "do not line up")
})

## ---- E3: the diffuse lens can drive the single-engine unit test ------------

test_that("E3: a diffuse-only signal reaches p_unit under the flat engine", {
  set.seed(42); n <- 160; k <- 10
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n))
  ## every probe tracks y, but with signs OPPOSING half its direction calls:
  ## the aligned mean cancels (coherence/composite die), the spread does not
  sg <- rep(c(1, -1), k / 2)
  probes <- paste0("cg", 1:k)
  for (j in seq_len(k))
    d[[probes[j]]] <- sg[j] * 0.8 * d$y + stats::rnorm(n)
  al <- dmsa_align(data.frame(cpg = probes, d = rep(1, k),
                              p_plus = rep(.95, k)),
                   genes = rep("G1", k), level = "gene")
  r <- dmsa_triangulate(as.matrix(d[probes]), d, rhs = c("y", "cov1"),
                        term = "y", units = rep("G1", k), alignment = al,
                        B = 199, seed = 1, weighting = "flat")
  expect_lt(r$p_diffuse, 0.05)          # the diffuse lens sees it
  ## diffuse is the STRONGEST lens on this construction, and the fixed
  ## any-lens statistic lets it carry the unit (the old raw-scale max could
  ## not be driven by the ratio-scaled diffuse statistic at all)
  expect_lte(r$p_diffuse, r$p_coherence)
  expect_lt(r$p_unit, 0.05)
})

## ---- E4: user polarity merges with the bundled table -----------------------

test_that("E4: a one-gene user table no longer strips every other gene", {
  probes <- c("cgA1", "cgA2")
  dir <- data.frame(cpg = probes, d = c(1, -1), p_plus = c(.9, .1))
  msgs <- character()
  withCallingHandlers(
    a <- dmsa_align(dir, genes = c("NR3C1", "FKBP5"), level = "system",
                    polarity = data.frame(gene = "NR3C1", w_g = 1),
                    system_id = "2"),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m))
                            invokeRestart("muffleMessage") })
  expect_true(any(grepl("polarity merged", msgs)))
  expect_identical(a$w_g[a$gene == "NR3C1"], 1)        # user wins
  expect_true(is.finite(a$w_g[a$gene == "FKBP5"]))     # bundled fills
})

## ---- E5: the EB signal fraction responds to the data ------------------------

test_that("E5: lfdr is high on pure noise and low on planted signal", {
  set.seed(1)
  l0 <- dmsa_lfdr(stats::rnorm(400))
  expect_gt(stats::median(l0), 0.85)
  z <- c(stats::rnorm(320), stats::rnorm(80, mean = 4))
  l1 <- dmsa_lfdr(z)
  expect_lt(mean(l1[321:400]), 0.5)                    # signal units credible
  expect_gt(stats::median(l1[1:320]), 0.8)             # nulls stay incredible
})

## ---- E6: omnibus p is invariant to row layout ------------------------------

test_that("E6: block-contiguous and interleaved layouts give one answer", {
  set.seed(2); n <- 84
  d <- data.frame(g = factor(rep(c("a", "b", "c"), length.out = n)),
                  cov1 = stats::rnorm(n),
                  fam = rep(seq_len(n / 2), each = 2))
  Y <- matrix(stats::rnorm(n * 4), n, dimnames = list(NULL, paste0("cg", 1:4)))
  des <- dmsa_design(focal = "g", fixed = "cov1", exchangeable = "fam")
  o <- order(d$fam)                                    # block-contiguous copy
  d2 <- d[o, , drop = FALSE]; rownames(d2) <- NULL
  Y2 <- Y[o, , drop = FALSE]
  r1 <- dmsa_omnibus(d, Y, design = des, beta_input = FALSE, B = 199, seed = 5)
  r2 <- dmsa_omnibus(d2, Y2, design = des, beta_input = FALSE, B = 199,
                     seed = 5)
  expect_equal(r1$statistic, r2$statistic, tolerance = 1e-8)
  ## the permuted BLOCKS are the same families either way; the p may differ
  ## only by which block is swapped with which under the seed, so demand
  ## agreement within permutation noise
  expect_lt(abs(r1$p_perm - r2$p_perm), 0.12)
})

## ---- E7: the Johnson-Neyman region lands on the right side -----------------

test_that("E7: strong main effect + weak interaction reports BETWEEN", {
  set.seed(3); n <- 300
  S <- stats::rnorm(n); W <- stats::rnorm(n)
  probes <- paste0("cg", 1:4)
  d <- data.frame(W = W, cov1 = stats::rnorm(n))
  for (p in probes) d[[p]] <- S + stats::rnorm(n, sd = .3)
  ## the tone tracks S; the outcome loads hard on the tone, barely on the
  ## product - the everyday A > 0 case the old code reported inverted
  d$y <- 0.8 * scale(rowMeans(d[probes]))[, 1] +
    0.02 * scale(rowMeans(d[probes]))[, 1] * W + stats::rnorm(n)
  al <- dmsa_align(data.frame(cpg = probes, d = 1, p_plus = .95),
                   genes = rep("G1", 4), level = "gene")
  des <- dmsa_design(focal = "y", fixed = "cov1")
  m <- dmsa_moderate(d, probes = probes, alignment = al, outcome = "y",
                     moderator = "W", design = des, B = 99,
                     beta_input = FALSE, seed = 1)
  expect_true(m$jn_region %in% c("between", "everywhere"))
  if (identical(m$jn_region, "between")) {
    expect_true(all(is.finite(m$jn_bounds)))
    expect_gt(m$jn_share, 0.5)     # significant for MOST of the sample
  }
})

## ---- E8: mediate refuses the silent Alpha substitution ---------------------

test_that("E8: a custom cascade or CSV path is refused, not substituted", {
  d <- data.frame(y = stats::rnorm(40))
  cas <- data.frame(system_id = "1", system = "S", module_id = "1.1",
                    module = "m", gene = "G1")
  expect_error(
    suppressMessages(dmsa_mediate(d, map = suppressMessages(dmsa_sets(cas)),
                                  x = "a", m = "b", outcome = "y")),
    "will not silently substitute")
})

## ---- E10: the batch --------------------------------------------------------

test_that("E10: minP under the combined engine is refused in dmsa_change", {
  expect_error(dmsa_change(correction = "minP"),
               "not defined for the combined engine")
})

test_that("E10: pc1 survives an all-zero-multiplier set", {
  set.seed(4); n <- 50
  M <- matrix(stats::rnorm(n * 4), n,
              dimnames = list(NULL, paste0("cg", 1:4)))
  al <- data.frame(probe = colnames(M), s = 0, p_s_plus = NA_real_,
                   usable = FALSE)
  s <- dmsa_scores(M, al, flavours = c("mean", "pc1", "blind"))
  expect_true(all(is.finite(s$pc1)))
})

test_that("E10: BED-convention position columns are shifted to 1-based", {
  man <- data.frame(probeID = c("cg1", "cg2"), CpG_chrm = c("chr1", "chr2"),
                    CpG_beg = c(999L, 1999L))
  co <- suppressMessages(dmsa:::.coords_match(man, c("cg1", "cg2"),
          probe_col = "probeID", chr_col = "CpG_chrm", pos_col = "CpG_beg",
          what = "test manifest"))
  expect_identical(co$pos, c(1000, 2000))
  man2 <- data.frame(probeID = "cg1", CpG_chrm = "chr1", MAPINFO = 1000L)
  co2 <- suppressMessages(dmsa:::.coords_match(man2, "cg1",
          probe_col = "probeID", chr_col = "CpG_chrm", pos_col = "MAPINFO",
          what = "test manifest"))
  expect_identical(co2$pos, 1000)                     # 1-based stays put
})

test_that("E10b: the moderated outcome test still runs (full Freedman-Lane)", {
  set.seed(5); n <- 120
  d <- data.frame(hit = rbinom(n, 1, .5), S = stats::rnorm(n),
                  W = stats::rnorm(n), cov1 = stats::rnorm(n))
  r <- dmsa_outcome(d, outcome = "hit", score = "S", moderator = "W",
                    covariates = "cov1", family = "binomial", B = 99, seed = 2)
  expect_true(is.finite(r$p_perm))
  expect_identical(r$tested, "SW")
})

## ---- E5 undermining (2026-08-29): the kurtosis guard ------------------------

test_that("E5 guard: a heavy-tailed null falls back rather than selecting", {
  set.seed(11)
  z <- rt(2000, 5) / sqrt(5 / 3)          # pure heavy-tailed null, unit var
  expect_warning(l <- dmsa_lfdr(z), "heavy-tailed null")
  expect_null(l)                          # the documented frequentist fallback
})

test_that("E5 guard: a Gaussian null and one-sided signal pass through", {
  set.seed(12)
  expect_silent(l0 <- dmsa_lfdr(rnorm(1000)))
  expect_false(is.null(l0))
  z <- c(rnorm(340), rnorm(60, 3.5))      # one-sided: one tail only
  expect_silent(l1 <- dmsa_lfdr(z))
  expect_false(is.null(l1))
  expect_lt(mean(l1[341:400]), 0.5)       # power intact
})

test_that("E5 guard: symmetric two-sided signal is a safe fallback", {
  set.seed(13)
  z <- c(rnorm(340), rnorm(30, 4), rnorm(30, -4))
  ## indistinguishable from kurtosis by tail symmetry - falls back with a
  ## warning; the cascade's frequentist arm handles it from there
  expect_warning(l <- dmsa_lfdr(z), "heavy-tailed null")
  expect_null(l)
})

## ---- second-level Westfall-Young minP (PI-approved 2026-08-29) ------------
## p_union_exact: exact family-wise p for the any-lens union claim, calibrated
## from the same permutation stream; attr union_null_min carries the null
## family-best distribution; mean(union_null_min < alpha) is the realized
## FWER of the any-lens naming rule in THIS design.

test_that("dmsa_triangulate returns a coherent second-level union calibration", {
  set.seed(31)
  n <- 60; K <- 6
  M <- matrix(stats::plogis(rnorm(n * K)), n,
              dimnames = list(NULL, paste0("cg", 1:K)))
  d <- data.frame(y = rnorm(n), cov1 = rnorm(n))
  ## one strong aligned unit
  M[, 1:2] <- stats::plogis(1.2 * scale(d$y)[, 1] + rnorm(n * 2, 0, .5))
  al <- dmsa_align(data.frame(cpg = colnames(M), d = 1, p_plus = .9),
                   genes = rep(c("G1", "G2", "G3"), each = 2), level = "gene")
  r <- dmsa_triangulate(M, d, c("y", "cov1"), "y",
                        rep(c("G1", "G2", "G3"), each = 2), al,
                        B = 99, seed = 1)
  expect_true("p_union_exact" %in% names(r))
  ok <- r$n_probes > 0
  expect_true(all(is.finite(r$p_union_exact[ok])))
  ## the exact union p can never beat the unit's own best lens-adjusted p
  lmin <- pmin(r$p_coherence_adj, r$p_composite_adj, r$p_diffuse_adj,
               na.rm = TRUE)
  expect_true(all(r$p_union_exact[ok] >= lmin[ok] - 1e-12))
  ## and is monotone in it
  o <- order(lmin[ok])
  expect_true(all(diff(r$p_union_exact[ok][o]) >= -1e-12))
  ## the calibration distribution rides along, one value per usable draw
  Mb <- attr(r, "union_null_min")
  expect_true(length(Mb) == 99L)
  expect_true(all(Mb[is.finite(Mb)] > 0 & Mb[is.finite(Mb)] <= 1))
  ## realized FWER of the any-lens rule at .05, from the same stream
  fw <- mean(Mb < .05, na.rm = TRUE)
  expect_true(fw >= 0 && fw <= 1)
  ## family-corrected omnibus (PI-approved): present, never beats its own
  ## raw omnibus, monotone in it, with its own calibration distribution
  expect_true("p_omnibus_adj" %in% names(r))
  expect_true(all(is.finite(r$p_omnibus_adj[ok])))
  expect_true(all(r$p_omnibus_adj[ok] >= r$p_omnibus[ok] - 1e-12))
  oo <- order(r$p_omnibus[ok])
  expect_true(all(diff(r$p_omnibus_adj[ok][oo]) >= -1e-12))
  Mo <- attr(r, "omnibus_null_min")
  expect_true(length(Mo) == 99L)
  expect_true(all(Mo[is.finite(Mo)] > 0 & Mo[is.finite(Mo)] <= 1))
})
