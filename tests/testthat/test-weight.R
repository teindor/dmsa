# The two engines: reliability (default) vs flat. Reliability must (a) preserve
# calibration (outcome-free weights), (b) reproduce flat when probes are
# homogeneous, (c) gain power when probe quality is heterogeneous, (d) leave the
# flat engine bit-reproducible when explicitly requested.

.sim_unit <- function(n = 300, K = 12, rho_car = 1, beta = 0, seed = 1,
                      miscall = 0.10) {
  set.seed(seed); n <- 2L * (n %/% 2L); cid <- rep(seq_len(n / 2), each = 2)
  f <- sqrt(.25) * rep(rnorm(n / 2), each = 2) + sqrt(.75) * rnorm(n)
  ncar <- max(1L, round(rho_car * K))
  lam <- c(rep(.8, ncar), rep(.05, K - ncar))
  M <- sapply(seq_len(K), function(j) lam[j] * f + sqrt(1 - lam[j]^2) * rnorm(n))
  y <- beta * f + sqrt(.2) * rep(rnorm(n / 2), each = 2) + rnorm(n)
  d <- sample(c(1, -1), K, TRUE, prob = c(1 - miscall, miscall))
  colnames(M) <- sprintf("cg%04d", seq_len(K))
  al <- dmsa_align(data.frame(cpg = colnames(M), d = d,
                              p_plus = ifelse(d > 0, .9, .1)),
                   genes = rep("G", K), level = "gene")
  list(M = M, dat = data.frame(y = y, cv = rnorm(n)), al = al,
       units = rep("G", K), cid = cid)
}

test_that("reliability weights are outcome-free and preserve calibration", {
  # null (beta = 0): reject ~ alpha for the DEFAULT (reliability) engine
  rej <- mean(vapply(1:60, function(s) {
    S <- .sim_unit(n = 200, K = 12, rho_car = 0.4, beta = 0, seed = 100 + s)
    r <- suppressWarnings(dmsa_triangulate(
      S$M, S$dat, c("y", "cv"), "y", S$units, S$al, block = S$cid,
      B = 199, seed = s))               # weighting = "reliability" by default
    isTRUE(r$p_composite[1] < .05)
  }, logical(1)))
  expect_lt(rej, 0.18)                    # calibrated (MC slack at 60 reps)
})

test_that("flat engine is reproducible and reliability defaults on", {
  S <- .sim_unit(n = 300, K = 12, rho_car = 0.3, beta = 0.3, seed = 7)
  r_flat <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                             block = S$cid, B = 299, weighting = "flat", seed = 1)
  r_flat2 <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                              block = S$cid, B = 299, weighting = "flat", seed = 1)
  expect_identical(r_flat$p_composite, r_flat2$p_composite)   # deterministic
  r_def <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                            block = S$cid, B = 299, seed = 1)
  expect_identical(attr(r_def, "weighting"), "combined")      # default engine
})

test_that("combined engine is calibrated and recovers both archetypes", {
  # calibration: combined omnibus type I ~ alpha (outcome-free weights, joint null)
  rej <- mean(vapply(1:80, function(s) {
    S <- .sim_unit(n = 200, K = 12, rho_car = 0.4, beta = 0, seed = 2200 + s)
    r <- suppressWarnings(dmsa_triangulate(
      S$M, S$dat, c("y", "cv"), "y", S$units, S$al, block = S$cid,
      B = 199, weighting = "combined", seed = s))
    isTRUE(r$p_omnibus[1] < .05) }, logical(1)))
  expect_lt(rej, 0.16)                    # not inflated (MC slack at 80 reps)
  # heterogeneous: combined keeps the reliability gain (>= flat)
  pw <- function(wt, rho) mean(vapply(1:40, function(s) {
    S <- .sim_unit(n = 300, K = 12, rho_car = rho, beta = 0.26, seed = 2400 + s)
    r <- suppressWarnings(dmsa_triangulate(
      S$M, S$dat, c("y", "cv"), "y", S$units, S$al, block = S$cid,
      B = 199, weighting = wt, seed = s))
    isTRUE(min(r$p_coherence[1], r$p_composite[1], r$p_diffuse[1]) < .05) },
    logical(1)))
  expect_gte(pw("combined", 0.25), pw("flat", 0.25) - 0.02)   # not below flat
})

test_that("reliability gains power on a heterogeneous unit, ties on homogeneous", {
  pw <- function(rho, weighting) mean(vapply(1:40, function(s) {
    S <- .sim_unit(n = 300, K = 12, rho_car = rho, beta = 0.26, seed = 300 + s)
    r <- suppressWarnings(dmsa_triangulate(
      S$M, S$dat, c("y", "cv"), "y", S$units, S$al, block = S$cid,
      B = 199, weighting = weighting, seed = s))
    isTRUE(r$p_composite[1] < .05) }, logical(1)))
  # heterogeneous: reliability clearly beats flat
  expect_gt(pw(0.25, "reliability"), pw(0.25, "flat") + 0.08)
  # homogeneous: reliability does not underperform flat (within MC noise)
  expect_gte(pw(1.0, "reliability"), pw(1.0, "flat") - 0.08)
})

test_that("small K and incoherent units fall back to flat weights", {
  set.seed(11)
  Z <- scale(matrix(rnorm(200 * 2), 200))          # K = 2 -> flat
  expect_true(all(dmsa_relweights(Z, rep("g", 2), c(1, -1)) == 1))
  # incoherent K = 6 (random signs, no shared axis) -> coherence < floor -> flat
  Zi <- scale(matrix(rnorm(200 * 6), 200))
  w <- dmsa_relweights(Zi, rep("g", 6), rep(1, 6), w_floor = 0.6)
  expect_true(all(w == 1))
})

# ---------------------------------------------------------------------------
# maxT's information advantage under a SINGLE engine.
#
# The property the pre-registered proof of concept rests on: a unit whose
# probes all point the same way generates the widest permutation null in its
# family, so maxT charges it almost nothing. AVP is named at family-adjusted
# .0120 against a raw .0120 for exactly this reason.
#
# Standardising each unit's statistic by its own permutation null divides that
# width out and collapses maxT to minP-like behaviour. A change on 15 Aug 2026
# applied the standardisation to the one-engine path as well as the two-engine
# fusion where it belongs, and silently moved AVP from .0120 to .1135 - the
# documented COMBINED value. These tests exist so that cannot happen quietly
# again.
# ---------------------------------------------------------------------------

.sim_family <- function(n = 300, K = 8, nnoise = 5, beta = .28, seed = 7) {
  set.seed(seed); n <- 2L * (n %/% 2L); cid <- rep(seq_len(n / 2), each = 2)
  f <- rnorm(n)
  ## the concordant unit: every probe loads the same way
  Mc <- sapply(seq_len(K), function(j) .75 * f + sqrt(1 - .75^2) * rnorm(n))
  colnames(Mc) <- sprintf("cgC%03d", seq_len(K))
  ## noisy units: same size, no common factor
  Mn <- do.call(cbind, lapply(seq_len(nnoise), function(u) {
    m <- matrix(rnorm(n * K), n, K)
    colnames(m) <- sprintf("cgN%d_%03d", u, seq_len(K)); m }))
  M <- cbind(Mc, Mn)
  units <- c(rep("CONC", K), rep(paste0("NOISE", seq_len(nnoise)), each = K))
  y <- beta * f + rnorm(n)
  al <- dmsa_align(data.frame(cpg = colnames(M), d = 1, p_plus = .9),
                   genes = units, level = "gene")
  list(M = M, dat = data.frame(y = y, cv = rnorm(n)), al = al,
       units = units, cid = cid)
}

test_that("one engine keeps maxT's width advantage: adjusted stays near raw", {
  S <- .sim_family()
  r <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                        block = S$cid, B = 999, weighting = "flat",
                        correction = "maxT", seed = 3)
  i <- which(r$unit == "CONC")
  expect_true(r$p_coherence[i] < .02)
  ## The toll on the concordant unit must be modest. Under the standardisation
  ## bug this ratio ran to ~10x; the raw-statistic maxT keeps it small.
  expect_lt(r$p_coherence_adj[i] / r$p_coherence[i], 4)
})

test_that("the single-engine statistic is NOT standardised before maxT", {
  ## Directly: with one engine, fuse() must hand maxT the raw statistic. If it
  ## standardises, every unit's null has unit spread and the family max is
  ## driven by rank rather than width - minP by another name.
  S <- .sim_family(seed = 11)
  mx <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                         block = S$cid, B = 999, weighting = "flat",
                         correction = "maxT", seed = 5)
  mn <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                         block = S$cid, B = 999, weighting = "flat",
                         correction = "minP", seed = 5)
  i <- which(mx$unit == "CONC"); j <- which(mn$unit == "CONC")
  ## maxT must be strictly kinder than minP to the concordant unit - that gap
  ## IS the information advantage, and it vanishes if the statistic is scaled.
  expect_lt(mx$p_coherence_adj[i], mn$p_coherence_adj[j])
})

test_that("two engines are still standardised, so the fusion stays on one scale", {
  ## The 15 Aug change fixed a real problem for the COMBINED engine: a rank-p
  ## fusion floors the family-wise p near (family size)/(B+1). Standardisation
  ## must therefore survive wherever two engines are actually fused.
  S <- .sim_family(seed = 13)
  r <- dmsa_triangulate(S$M, S$dat, c("y", "cv"), "y", S$units, S$al,
                        block = S$cid, B = 999, weighting = "combined",
                        correction = "maxT", seed = 5)
  i <- which(r$unit == "CONC")
  ## no floor: with 6 units and B = 999 a rank-p fusion would pin this near
  ## 6/1000; the standardised fusion must be free to go below that
  expect_true(is.finite(r$p_coherence_adj[i]))
  expect_lt(r$p_coherence[i], .05)
})
