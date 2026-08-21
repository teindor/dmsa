test_that("an engine that reads only its own unit cannot dance", {
  set.seed(1)
  n <- 120L; K <- 40L
  f <- rnorm(n)
  M <- sapply(seq_len(K), function(j) if (j <= 6) 0.7 * f + rnorm(n) else rnorm(n))
  colnames(M) <- sprintf("cg%04d", seq_len(K))
  y <- 0.6 * f + rnorm(n)
  focal <- colnames(M)[1:6]

  competitive <- function(u) {
    rest <- setdiff(u, focal)
    if (length(rest) < 2) return(NA_real_)
    t.test(abs(cor(M[, focal], y)), abs(cor(M[, rest], y)))$p.value
  }
  unit_only <- function(u)
    summary(lm(rowMeans(M[, focal]) ~ y))$coefficients["y", 4]

  pd <- dmsa_pdance(list(competitive = competitive, unit_only = unit_only),
                    set = colnames(M), focal = focal,
                    grid = c(.1, .3, .5), R = 20, seed = 1)

  spread <- tapply(pd$p, pd$engine, function(x) diff(range(x, na.rm = TRUE)))
  ## the whole claim: unit-only is EXACTLY invariant, not merely stabler
  expect_identical(unname(spread[["unit_only"]]), 0)
  expect_gt(spread[["competitive"]], 0)
  expect_equal(sum(pd$flip[pd$engine == "unit_only"]), 0)
  expect_gt(sum(pd$flip[pd$engine == "competitive"]), 0)
})

test_that("the addition form moves a competitive engine too", {
  set.seed(2)
  n <- 100L
  f <- rnorm(n); y <- 0.6 * f + rnorm(n)
  M <- cbind(sapply(1:6, function(j) 0.7 * f + rnorm(n)),
             matrix(rnorm(n * 34), n, 34))
  colnames(M) <- sprintf("cg%03d", seq_len(40))
  extra <- matrix(rnorm(n * 40), n, 40,
                  dimnames = list(NULL, sprintf("zz%03d", seq_len(40))))
  A <- cbind(M, extra)
  focal <- colnames(M)[1:6]
  comp <- function(u) {
    rest <- setdiff(u, focal)
    if (length(rest) < 2) return(NA_real_)
    t.test(abs(cor(A[, focal], y)), abs(cor(A[, rest], y)))$p.value
  }
  pa <- dmsa_pdance(list(competitive = comp), set = colnames(M), focal = focal,
                    pool = colnames(extra), form = "addition",
                    grid = c(5, 20), R = 10, seed = 3)
  expect_true(all(pa$form == "addition"))
  expect_gt(diff(range(pa$p, na.rm = TRUE)), 0)
})

test_that("the contract is enforced, and the caller's RNG survives", {
  eng <- list(e = function(u) 0.5)
  expect_error(dmsa_pdance(eng, set = letters[1:5], focal = "z"), "subset")
  expect_error(dmsa_pdance(eng, set = letters[1:5], focal = letters[1:5]),
               "nothing to perturb")
  expect_error(dmsa_pdance(list(function(u) .5), set = letters[1:5],
                           focal = "a"), "named")
  expect_error(dmsa_pdance(eng, set = letters[1:5], focal = "a",
                           form = "addition"), "pool")
  expect_error(dmsa_pdance(eng, set = letters[1:5], focal = "a",
                           pool = c("e", "f"), form = "addition"), "disjoint")

  set.seed(99); before <- .Random.seed
  invisible(dmsa_pdance(eng, set = letters[1:9], focal = "a", grid = .2,
                        R = 3, seed = 7))
  expect_identical(before, .Random.seed)
})
