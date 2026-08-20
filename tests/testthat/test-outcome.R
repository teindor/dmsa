mk <- function(n = 240, family = "binomial", beta = 0.8, seed = 1) {
  set.seed(seed)
  S <- rnorm(n); cv <- rnorm(n); cid <- rep(seq_len(n / 2), each = 2)
  eta <- beta * S
  y <- switch(family,
    gaussian    = eta + rnorm(n),
    binomial    = rbinom(n, 1, plogis(eta - 0.2)),
    ordinal     = cut(eta + rlogis(n), c(-Inf, -1, 0.5, Inf),
                      labels = c("lo", "mid", "hi"), ordered_result = TRUE),
    multinomial = factor(apply(cbind(rlogis(n), 0.9 * eta + rlogis(n),
                                     -0.9 * eta + rlogis(n)), 1, which.max),
                         levels = 1:3, labels = c("A", "B", "C")),
    poisson     = rpois(n, exp(0.5 + eta)))
  data.frame(y = y, S = S, cv = cv, cid = cid)
}

test_that("binomial returns log-odds with a direction and an odds ratio", {
  o <- dmsa_outcome(mk(family = "binomial"), "y", "S", family = "binomial",
                    covariates = "cv", block = "cid", B = 99, seed = 2)
  expect_s3_class(o, "dmsa_outcome")
  expect_true(o$directional)
  expect_gt(o$estimate, 0)
  expect_lt(o$p_perm, .05)
  out <- capture.output(print(o))
  expect_true(any(grepl("log-odds", out)))
  expect_true(any(grepl("exp\\(coef\\)", out)))
})

test_that("ordinal keeps one coefficient and therefore one direction", {
  o <- dmsa_outcome(mk(family = "ordinal"), "y", "S", family = "ordinal",
                    covariates = "cv", block = "cid", B = 99, seed = 3)
  expect_true(o$directional)
  expect_equal(o$n_levels, 3L)
  expect_gt(o$estimate, 0)
  expect_lt(o$p_perm, .05)
})

test_that("multinomial refuses to claim a single direction and reports per category", {
  o <- dmsa_outcome(mk(family = "multinomial"), "y", "S", family = "multinomial",
                    covariates = "cv", block = "cid", B = 99, seed = 4)
  expect_false(o$directional)
  expect_true(is.na(o$estimate))
  expect_s3_class(o$per_category, "data.frame")
  expect_equal(nrow(o$per_category), 2L)      # K - 1 contrasts
  out <- capture.output(print(o))
  expect_true(any(grepl("NO single direction", out)))
  expect_true(any(grepl("family = 'ordinal'", out)))
})

test_that("poisson works and reports a rate ratio", {
  o <- dmsa_outcome(mk(family = "poisson", beta = 0.4), "y", "S",
                    family = "poisson", covariates = "cv", B = 99, seed = 5)
  expect_true(o$directional)
  out <- capture.output(print(o))
  expect_true(any(grepl("log-rate", out)))
})

test_that("a moderated non-Gaussian test tests the product term", {
  d <- mk(family = "binomial", beta = 0)
  set.seed(9); d$W <- rnorm(nrow(d))
  d$y <- rbinom(nrow(d), 1, plogis(0.9 * d$S * d$W - 0.2))
  o <- dmsa_outcome(d, "y", "S", family = "binomial", moderator = "W",
                    covariates = "cv", block = "cid", B = 199, seed = 6)
  expect_equal(o$tested, "SW")
  expect_lt(o$p_perm, .05)
  expect_true(o$centered)
})

test_that("outcome coding is validated with an actionable message", {
  d <- mk(family = "gaussian")
  expect_error(dmsa_outcome(d, "y", "S", family = "binomial", B = 0),
               "must be 0/1, logical, or a 2-level factor")
  d2 <- mk(family = "binomial"); d2$y <- factor(d2$y)
  expect_error(dmsa_outcome(d2, "y", "S", family = "ordinal", B = 0),
               ">= 3 levels")
  d3 <- mk(family = "gaussian"); d3$y <- d3$y - min(d3$y) + 0.5
  expect_error(dmsa_outcome(d3, "y", "S", family = "poisson", B = 0),
               "non-negative integers")
  expect_error(dmsa_outcome(d, "y", "NOPE", family = "gaussian", B = 0),
               "absent from data")
})

test_that("a rare outcome class warns rather than failing silently", {
  d <- mk(family = "gaussian"); set.seed(7)
  d$y <- c(rep(1L, 5), rep(0L, nrow(d) - 5))
  expect_warning(dmsa_outcome(d, "y", "S", family = "binomial", B = 0),
                 "rarer outcome class")
})

test_that("a two-level factor and 0/1 give the same answer", {
  d <- mk(family = "binomial")
  a <- dmsa_outcome(d, "y", "S", family = "binomial", covariates = "cv",
                    B = 0)
  d$y <- factor(d$y, labels = c("no", "yes"))
  b <- dmsa_outcome(d, "y", "S", family = "binomial", covariates = "cv",
                    B = 0)
  expect_equal(a$estimate, b$estimate, tolerance = 1e-8)
})

test_that("B = 0 gives model inference only", {
  o <- dmsa_outcome(mk(), "y", "S", family = "binomial", B = 0)
  expect_true(is.na(o$p_perm))
  expect_true(is.finite(o$p_model))
  out <- capture.output(print(o))
  expect_true(any(grepl("permutation not run", out)))
})

test_that("gaussian via dmsa_outcome agrees with lm on the same model", {
  d <- mk(family = "gaussian", beta = 0.5)
  o <- dmsa_outcome(d, "y", "S", family = "gaussian", covariates = "cv", B = 0)
  l <- coef(summary(lm(y ~ S + cv, d)))["S", 1]
  expect_equal(o$estimate, l, tolerance = 1e-8)
})
