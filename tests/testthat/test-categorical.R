mkcat <- function(n = 120, k = 3, seed = 1) {
  set.seed(seed)
  lv <- c("secure", "anxious", "avoidant", "disorganized")[seq_len(k)]
  data.frame(style = factor(sample(lv, n, TRUE), levels = lv),
             age = rnorm(n), cID = rep(seq_len(n / 2), each = 2))
}

des3 <- function(d) dmsa_design(focal = "style", fixed = "age",
                                exchangeable = "cID")

test_that("a 3-level exposure becomes two directional contrasts", {
  d <- mkcat()
  cc <- dmsa_contrasts(d, des3(d))
  expect_s3_class(cc, "dmsa_contrasts")
  expect_equal(length(cc$designs), 2L)
  expect_equal(cc$reference, "secure")
  expect_setequal(cc$contrasts, c("anxious", "avoidant"))
  ## each expanded design has exactly ONE focal column, which is what makes it
  ## directional
  for (dd in cc$designs) expect_equal(length(dd$focal_vars), 1L)
})

test_that("the contrast column excludes the third group instead of pooling it", {
  d <- mkcat()
  cc <- dmsa_contrasts(d, des3(d))
  cn <- cc$designs[["anxious"]]$focal_vars
  col <- cc$data[[cn]]
  expect_true(all(is.na(col[d$style == "avoidant"])))
  expect_equal(levels(col), c("secure", "anxious"))
  expect_equal(sum(!is.na(col)), sum(d$style %in% c("secure", "anxious")))
})

test_that("focal_test names a real model-matrix column of the expanded design", {
  d <- mkcat()
  cc <- dmsa_contrasts(d, des3(d))
  for (dd in cc$designs) {
    sub <- cc$data[!is.na(cc$data[[dd$focal_vars]]), , drop = FALSE]
    X <- model.matrix(reformulate(c(dd$focal, dd$fixed)), sub)
    expect_true(dd$focal_test %in% colnames(X))
  }
})

test_that("the reference level can be chosen and is validated", {
  d <- mkcat()
  cc <- dmsa_contrasts(d, des3(d), reference = "avoidant")
  expect_equal(cc$reference, "avoidant")
  expect_setequal(cc$contrasts, c("secure", "anxious"))
  expect_error(dmsa_contrasts(d, des3(d), reference = "nope"),
               "not among the levels")
})

test_that("a 2-level exposure is refused with the reason", {
  d <- mkcat(k = 2)
  expect_error(dmsa_contrasts(d, des3(d)),
               "already a single directional contrast")
})

test_that("a continuous focal term is refused rather than binned", {
  d <- mkcat(); d$style <- rnorm(nrow(d))
  expect_error(dmsa_contrasts(d, des3(d)), "looks\\s+continuous")
})

test_that("an interaction focal term is refused with instructions", {
  d <- mkcat(); d$time <- rep(0:1, nrow(d) / 2)
  des <- dmsa_design(focal = "time:style", fixed = "age", exchangeable = "cID")
  expect_error(dmsa_contrasts(d, des), "is an interaction")
})

test_that("a thin cell warns", {
  d <- mkcat(n = 200)
  d$style <- factor(c(rep("secure", 100), rep("anxious", 95),
                      rep("avoidant", 5)),
                    levels = c("secure", "anxious", "avoidant"))
  expect_warning(dmsa_contrasts(d, des3(d)), "fewer than 20")
})

test_that("holm is the default correction and it is applied", {
  d <- mkcat()
  cc <- dmsa_contrasts(d, des3(d))
  expect_equal(cc$correct, "holm")
  a <- dmsa_contrast_adjust(cc, c(0.01, 0.04))
  expect_equal(a$p_adj, p.adjust(c(0.01, 0.04), "holm"))
  expect_equal(nrow(a), 2L)
  expect_error(dmsa_contrast_adjust(cc, 0.01), "one p-value per contrast")
})

test_that("p-values supplied by name are reordered to the contrast order", {
  d <- mkcat()
  cc <- dmsa_contrasts(d, des3(d))
  p <- c(avoidant = 0.20, anxious = 0.01)
  a <- dmsa_contrast_adjust(cc, p)
  expect_equal(a$p[1], unname(p[cc$contrasts[1]]))
})

test_that("the omnibus runs, is calibrated, and refuses to claim a direction", {
  d <- mkcat(n = 100)
  set.seed(4); Y <- matrix(plogis(rnorm(100 * 12)), 100, 12)
  o <- dmsa_omnibus(d, Y, des3(d), B = 99, seed = 5)
  expect_false(o$directional)
  expect_true(o$p_perm > 0 && o$p_perm <= 1)
  expect_equal(o$df1, 2L)
  expect_true(grepl("no sign", o$note))
})

test_that("the omnibus detects a factor effect that has no common direction", {
  d <- mkcat(n = 160, seed = 11)
  set.seed(12)
  Y <- matrix(rnorm(160 * 20), 160, 20)
  ## half the probes go up in 'anxious', the other half go down: no single
  ## direction exists, so only a directionless test can see it
  up <- d$style == "anxious"
  Y[up, 1:10]  <- Y[up, 1:10]  + 1.1
  Y[up, 11:20] <- Y[up, 11:20] - 1.1
  o <- dmsa_omnibus(d, plogis(Y), des3(d), B = 199, seed = 13)
  expect_lt(o$p_perm, .05)
})

test_that("the omnibus needs 3+ levels", {
  d <- mkcat(k = 2)
  set.seed(6); Y <- matrix(plogis(rnorm(nrow(d) * 8)), nrow(d), 8)
  expect_error(dmsa_omnibus(d, Y, des3(d), B = 9), "use dmsa_fit")
})
