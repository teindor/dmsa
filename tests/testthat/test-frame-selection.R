## Spec 12/13/14: cpg_map default is "full"; "confidence" is opt-in.
## Spec 47: a character `methylation` selection is honoured at any length.

test_that("cpg_map defaults to 'full', not 'confidence'", {
  expect_identical(eval(formals(dmsa_frame)$cpg_map)[1], "full")
  expect_setequal(eval(formals(dmsa_frame)$cpg_map), c("full", "confidence"))
})

.sel_fixture <- function(n = 60, k = 3) {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene == "FKBP5"], k)
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  list(data = d, cols = cols)
}

.run_sel <- function(fx, sel)
  dmsa_frame(fx$data, methylation = sel, direction_source = "bundled",
             outcome = "y", covariates = "cov1",
             blocks = NULL, chip = FALSE, B = 19, plots = FALSE,
             tables = FALSE, summary = FALSE, progress = FALSE, beep = FALSE,
             outdir = tempfile("dmsa_sel"))

test_that("a ONE-column methylation selection analyses exactly that probe", {
  set.seed(1); fx <- .sel_fixture()
  skip_if(length(fx$cols) < 2, "fixture columns unavailable")
  f <- .run_sel(fx, fx$cols[1])
  expect_identical(sort(unique(f$map$column)), sort(fx$cols[1]))
})

test_that("a multi-column selection still works and is not widened", {
  set.seed(1); fx <- .sel_fixture()
  skip_if(length(fx$cols) < 2, "fixture columns unavailable")
  f <- .run_sel(fx, fx$cols[1:2])
  expect_setequal(unique(f$map$column), fx$cols[1:2])
})
