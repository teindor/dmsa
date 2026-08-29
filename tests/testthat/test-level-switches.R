## Spec sections 21-24: `system` is a TRUE/FALSE analysis-level switch, never a
## biological selector. The plural arguments select biology.

fake_data <- function(n = 40)
  data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n))

test_that("a system NAME passed to the singular switch hard-errors", {
  d <- fake_data()
  for (v in list("oxytocin", "hpa", c("oxytocin", "sex_steroids"))) {
    e <- expect_error(dmsa_frame(d, outcome = "y", system = v))
    msg <- conditionMessage(e)
    expect_match(msg, "does NOT select", fixed = TRUE)
    expect_match(msg, "systems = ", fixed = TRUE)
    expect_match(msg, "No analysis was run", fixed = TRUE)
  }
})

test_that("each level switch points at its OWN plural argument", {
  d <- fake_data()
  expect_error(dmsa_frame(d, outcome = "y", gene   = "NR3C1"), "genes = ")
  expect_error(dmsa_frame(d, outcome = "y", module = "2.6"),   "modules = ")
  expect_error(dmsa_frame(d, outcome = "y", probe  = "cg1"),   "probes = ")
})

test_that("non-logical or non-scalar switches hard-error", {
  d <- fake_data()
  expect_error(dmsa_frame(d, outcome = "y", system = NA), "single TRUE or FALSE")
  expect_error(dmsa_frame(d, outcome = "y", system = c(TRUE, FALSE)), "length 2")
  expect_error(dmsa_frame(d, outcome = "y", system = 1), "single TRUE or FALSE")
})

test_that("other logical controls are validated too", {
  d <- fake_data()
  expect_error(dmsa_frame(d, outcome = "y", plots = "yes"), "single TRUE or FALSE")
  expect_error(dmsa_frame(d, outcome = "y", moderation = NA), "single TRUE or FALSE")
})

test_that("systems = 'all' cannot be combined with named systems", {
  d <- fake_data()
  expect_error(dmsa_frame(d, outcome = "y", systems = c("all", "hpa")),
               "cannot be combined")
})
