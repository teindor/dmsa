## Design-layer tests. Several of these encode real failures: a batch random
## effect silently omitted across an analysis family, a covariate that is a
## constant in one build, and a covariate that proxies sex.

mk <- function(n = 200, seed = 1) {
  set.seed(seed)
  fam <- rep(seq_len(n / 2), each = 2)
  data.frame(
    cID = fam,
    chip = rep(seq_len(n / 8), each = 8)[seq_len(n)],
    sex_c = rep(c(-0.5, 0.5), n / 2),
    age = rnorm(n, 40, 5),
    Epi = runif(n),
    Fib = rep(0, n),                      # constant, as in the child build
    sv3 = rnorm(n), sv5 = rnorm(n),
    svSex = rep(c(-0.5, 0.5), n / 2) + rnorm(n, 0, .2),   # proxies sex
    exposure = rnorm(n),
    time = rep(c(-0.5, 0.5), n / 2),
    y = rnorm(n)
  )
}

test_that("dmsa_design expands random specs both ways", {
  a <- dmsa_design("exposure", c("sex_c"), random = c("chip", "cID"), exchangeable = "cID")
  b <- dmsa_design("exposure", c("sex_c"), random = ~ (1 | chip) + (1 | cID), exchangeable = "cID")
  expect_equal(a$random_groups, c("chip", "cID"))
  expect_equal(sort(b$random_groups), sort(c("chip", "cID")))
})

test_that("a never-list term cannot enter as fixed or focal", {
  expect_error(
    dmsa_design("exposure", c("sex_c", "svSex"), forbid = "svSex"),
    "never-list"
  )
  expect_error(
    dmsa_design("svSex", c("sex_c"), forbid = "svSex"),
    "never-list"
  )
})

test_that("interaction focal terms are split into their variables", {
  d <- dmsa_design("time:exposure", "sex_c", exchangeable = "cID")
  expect_setequal(d$focal_vars, c("time", "exposure"))
  expect_true("time" %in% d$vars && "exposure" %in% d$vars)
})

test_that("a constant covariate is caught, not silently absorbed", {
  d <- dmsa_design("exposure", c("sex_c", "Fib"), random = "chip", exchangeable = "cID")
  expect_error(dmsa_check_design(d, mk()), "constant in this subsample")
})

test_that("a covariate that proxies another is caught by VIF", {
  d <- dmsa_design("exposure", c("sex_c", "svSex"), random = "chip", exchangeable = "cID")
  expect_error(dmsa_check_design(d, mk()), "VIF > 5")
})

test_that("declared independence is challenged when rows plainly repeat", {
  d <- dmsa_design("exposure", "sex_c", random = NULL, exchangeable = "cID")
  expect_error(dmsa_check_design(d, mk()), "no random effects declared")
})

test_that("a grouping factor aliased with the focal term is caught", {
  dd <- mk()
  dd$exposure <- dd$chip                     # exposure constant within chip
  d <- dmsa_design("exposure", "sex_c", random = "chip", exchangeable = "cID")
  expect_error(dmsa_check_design(d, dd), "aliased")
})

test_that("missing columns are named", {
  d <- dmsa_design("exposure", c("sex_c", "nope"), exchangeable = "cID")
  expect_error(dmsa_check_design(d, mk()), "nope")
})

test_that("a clean design passes", {
  d <- dmsa_design("exposure", c("sex_c", "age", "Epi", "sv3", "sv5"),
                   random = c("chip", "cID"), exchangeable = "cID")
  res <- dmsa_check_design(d, mk())
  expect_length(res$problems, 0)
  expect_equal(res$n, 200)
  expect_equal(unname(res$blocks["n_blocks"]), 100)
})

test_that("alpha_design carries the real build contracts and never-lists", {
  b1 <- alpha_design(1, focal = "IC_T1_c")
  expect_setequal(b1$fixed, c("sex_c", "age_at_array_T1", "Epi_T1", "Fib_T1",
                              "ctrlSV3_T1", "ctrlSV5_T1"))
  expect_setequal(b1$random_groups, c("chip_T1", "cID"))
  expect_true("ctrlSV4_T1" %in% b1$forbid)

  b4 <- alpha_design(4, focal = "MIBS_c")
  expect_setequal(b4$fixed, c("sex_c", "birth_week", "birth_weight", "Epi_T4"))
  expect_equal(b4$random_groups, "chip_T4")
  expect_true("Fib_T4" %in% b4$forbid)      # constant in every child
  expect_false("cID" %in% b4$random_groups) # build 4 declares no family random term

  b2 <- alpha_design(2, focal = "time:BSI_c")
  expect_true("submission_array" %in% b2$forbid)  # perfectly separates time
  expect_error(alpha_design(9, "x"), "must be 1, 2, 3 or 4")
})

test_that("a declared deviation is recorded, not silent", {
  dd <- alpha_design(1, focal = "IC_T1_c", drop = "Epi_T1")
  expect_false("Epi_T1" %in% dd$fixed)
  expect_equal(dd$dropped, "Epi_T1")
  expect_match(dd$label, "declared deviation")
})

test_that("block permutation keeps blocks whole and rows ordered", {
  blk <- rep(seq_len(50), each = 2)
  idx <- dmsa:::dmsa_block_index(blk, 5)
  for (i in idx) {
    expect_equal(sort(i), seq_along(blk))          # a permutation
    expect_equal(length(unique(blk[i])), 50)
    ## each destination block received exactly one source block
    tab <- tapply(blk[i], blk, function(z) length(unique(z)))
    expect_true(all(tab == 1))
  }
})

test_that("block permutation moves all rows of a longitudinal block together", {
  blk <- rep(seq_len(30), each = 4)                # family x person x time
  idx <- dmsa:::dmsa_block_index(blk, 3)
  for (i in idx) expect_true(all(tapply(blk[i], blk, function(z) length(unique(z))) == 1))
})
