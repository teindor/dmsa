## Fitting-engine tests: does dmsa_fit recover a planted aligned effect, does
## it keep its error rate under dependence, and does the moderated version
## behave under centering.

sim <- function(n = 240, K = 24, h = 0, seed = 11, fam = TRUE) {
  set.seed(seed)
  cid  <- if (fam) rep(seq_len(n / 2), each = 2) else seq_len(n)
  chip <- rep(seq_len(ceiling(n / 8)), each = 8)[seq_len(n)]
  x <- rnorm(n)
  if (fam) x <- x + rep(rnorm(n / 2, 0, .5), each = 2)   # family-level component
  d <- sample(c(-1, 1), K, TRUE)                          # true CpG -> expression sign
  M <- matrix(rnorm(n * K), n, K)
  for (k in seq_len(K)) M[, k] <- M[, k] + h * d[k] * x   # aligned effect
  beta <- 1 / (1 + 2^(-M))
  dd <- data.frame(cID = cid, chip = chip, x = x,
                   sex_c = rep(c(-.5, .5), length.out = n),
                   cov1 = rnorm(n))
  probes <- paste0("cg", seq_len(K))
  dd[probes] <- beta
  al <- dmsa_align(data.frame(cpg = probes, d = d, p_plus = ifelse(d > 0, .95, .05)),
                   genes = paste0("g", ceiling(seq_len(K) / 3)), level = "gene")
  list(data = dd, probes = probes, al = al)
}

des <- function() dmsa_design("x", c("sex_c", "cov1"),
                              random = c("chip", "cID"), exchangeable = "cID")

test_that("dmsa_fit recovers a planted aligned effect", {
  s <- sim(h = 0.35, seed = 3)
  f <- dmsa_fit(s$data, s$probes, s$al, des(), B = 199, engine = "lm", seed = 1)
  expect_s3_class(f, "dmsa_fit")
  expect_gt(f$z, 3)
  expect_lt(f$p_perm, .05)
  expect_equal(f$n, 240)
})

test_that("dmsa_fit is null when there is no effect", {
  s <- sim(h = 0, seed = 4)
  f <- dmsa_fit(s$data, s$probes, s$al, des(), B = 299, engine = "lm", seed = 2)
  expect_gt(f$p_perm, .05)
})

test_that("permutation respects the family block (calibrated under dependence)", {
  ps <- vapply(1:25, function(i) {
    s <- sim(h = 0, seed = 100 + i)
    dmsa_fit(s$data, s$probes, s$al, des(), B = 99, engine = "lm",
             seed = i, check = FALSE)$p_perm
  }, numeric(1))
  expect_lt(mean(ps < .05), 0.20)     # loose bound: no gross inflation
  expect_gt(mean(ps), 0.25)
})

test_that("the design check runs by default and blocks a bad design", {
  s <- sim(h = 0)
  s$data$const <- 1
  bad <- dmsa_design("x", c("sex_c", "const"), random = "chip", exchangeable = "cID")
  expect_error(dmsa_fit(s$data, s$probes, s$al, bad, B = 0), "constant")
})

test_that("a focal term absent from the model matrix is named, not silently dropped", {
  s <- sim(h = 0)
  d <- dmsa_design("not_a_column", "sex_c", exchangeable = "cID")
  expect_error(dmsa_fit(s$data, s$probes, s$al, d, B = 0), "absent from data")
})

test_that("engine='lmer' needs declared random effects", {
  s <- sim(h = 0)
  d <- dmsa_design("x", "sex_c", exchangeable = "cID")
  skip_if_not_installed("lme4")
  ## the design check fires first and is the better error; with it disabled the
  ## engine still refuses to invent a random structure that was not declared
  expect_error(dmsa_fit(s$data, s$probes, s$al, d, B = 0, engine = "lmer",
                        check = FALSE),
               "declares no random effects")
  expect_error(dmsa_fit(s$data, s$probes, s$al, d, B = 0, engine = "lmer"),
               "no random effects declared")
})

test_that("lmer and lm agree closely on the same data, and the null is lm-matched", {
  skip_if_not_installed("lme4")
  s <- sim(h = 0.3, seed = 7)
  f <- dmsa_fit(s$data, s$probes, s$al, des(), B = 99, engine = "auto", seed = 5)
  expect_match(f$engine, "lmer")
  expect_gt(f$agreement_lmer_lm, 0.9)
  expect_true(is.finite(f$z_lm))
  expect_length(f$null_z, 99)
})

test_that("B = 0 gives a descriptive fit with no permutation p", {
  s <- sim(h = 0.3)
  f <- dmsa_fit(s$data, s$probes, s$al, des(), B = 0, engine = "lm")
  expect_true(is.na(f$p_perm))
})

test_that("dmsa_score is standardised and orientation-aware", {
  s <- sim(h = 0.4, seed = 9)
  sc <- dmsa_score(s$data, s$probes, s$al)
  expect_equal(mean(sc), 0, tolerance = 1e-8)
  expect_equal(stats::sd(sc), 1, tolerance = 1e-8)
  expect_gt(abs(cor(sc, s$data$x)), 0.2)
})

test_that("mDMSA recovers a planted interaction and flags uncentered moderators", {
  set.seed(21)
  s <- sim(h = 0.5, seed = 12)
  sc <- dmsa_score(s$data, s$probes, s$al)
  s$data$mod <- rnorm(nrow(s$data), 10, 2)          # far from zero on purpose
  s$data$out <- 0.5 * sc * (s$data$mod - 10) + rnorm(nrow(s$data))
  d <- dmsa_design("placeholder", c("sex_c", "cov1"),
                   random = c("chip", "cID"), exchangeable = "cID")
  m <- dmsa_moderate(s$data, s$probes, s$al, outcome = "out", moderator = "mod",
                     design = d, B = 199, engine = "lm", seed = 3)
  expect_s3_class(m, "dmsa_moderate")
  expect_lt(m$p_perm, .05)
  expect_lt(m$vif_interaction, 5)                   # centered: product not collinear
  expect_true(m$centered)
  expect_equal(nrow(m$simple_slopes), 2)

  mu <- dmsa_moderate(s$data, s$probes, s$al, outcome = "out", moderator = "mod",
                      design = d, B = 0, engine = "lm", center = FALSE)
  expect_gt(mu$vif_interaction, m$vif_interaction)  # uncentering inflates it
  expect_equal(mu$b_int, m$b_int, tolerance = 1e-6) # but the interaction is unchanged
})

test_that("mDMSA refuses a never-list outcome or moderator", {
  s <- sim(h = 0.3)
  s$data$banned <- rnorm(nrow(s$data))
  d <- dmsa_design("placeholder", "sex_c", random = "chip",
                   exchangeable = "cID", forbid = "banned")
  expect_error(dmsa_moderate(s$data, s$probes, s$al, outcome = "banned",
                             moderator = "cov1", design = d, B = 0),
               "never-list")
})

test_that("Johnson-Neyman reports honestly when there is no region", {
  set.seed(31)
  s <- sim(h = 0, seed = 15)
  s$data$mod <- rnorm(nrow(s$data))
  s$data$out <- rnorm(nrow(s$data))
  d <- dmsa_design("placeholder", c("sex_c", "cov1"),
                   random = c("chip", "cID"), exchangeable = "cID")
  m <- dmsa_moderate(s$data, s$probes, s$al, outcome = "out", moderator = "mod",
                     design = d, B = 0, engine = "lm")
  expect_true(is.na(m$jn_share) || m$jn_share < 0.25)
})
