## Cascade tests. The substantive claims are checked in simulation elsewhere;
## these pin the machinery: level-specific engines, the sparse companion, the
## gate actually gating, and the local-fdr estimator refusing to guess.

mk_tree <- function(nsys = 6, nmod = 3, ngene = 4, nprobe = 6) {
  tr <- expand.grid(probe = seq_len(nprobe), gene = seq_len(ngene),
                    module = seq_len(nmod), system = seq_len(nsys))
  tr <- tr[order(tr$system, tr$module, tr$gene, tr$probe), ]
  data.frame(system = tr$system, module = tr$module, gene = tr$gene)
}
mk_data <- function(tree, effect_rows = integer(0), h = 0.4, seed = 1) {
  set.seed(seed); P <- nrow(tree)
  b <- rnorm(P, 0, .05); b[effect_rows] <- b[effect_rows] + h
  se <- rep(.05, P)
  al <- data.frame(probe = seq_len(P), gene = tree$gene,
                   s = 1L, p_s_plus = .95)
  list(b = b, se = se, al = al)
}
mk_nulls <- function(sd_by = c(system = 1, module = 1, gene = 1),
                     sparse_sd = 1.6, n = 3000) {
  set.seed(2)
  list(system = rnorm(n, 0, sd_by[["system"]]),
       module = rnorm(n, 0, sd_by[["module"]]),
       gene   = rnorm(n, 0, sd_by[["gene"]]),
       system_sparse = abs(rnorm(n, 0, sparse_sd)),
       module_sparse = abs(rnorm(n, 0, sparse_sd)))
}

test_that("dmsa_lfdr refuses to estimate from too few units", {
  expect_null(dmsa_lfdr(rnorm(10)))
  expect_null(dmsa_lfdr(rep(0, 200)))
  lf <- dmsa_lfdr(c(rnorm(400), rnorm(40, 5)))
  expect_length(lf, 440)
  expect_true(all(lf >= 0 & lf <= 1))
  ## the extreme units must get the lower local fdr
  expect_lt(mean(lf[401:440]), mean(lf[1:400]))
})

test_that("prior-odds multiplier lowers local fdr monotonically", {
  z <- c(rnorm(400), rnorm(40, 4))
  a <- dmsa_lfdr(z, 1); b <- dmsa_lfdr(z, 4)
  expect_true(mean(b) <= mean(a) + 1e-9)
})

test_that("dmsa_bfdr_select respects the mean-lfdr budget", {
  lf <- c(rep(.001, 10), rep(.5, 10), rep(.99, 10))
  s <- dmsa_bfdr_select(lf, .05)
  expect_true(mean(lf[s]) <= .05)
  expect_gte(length(s), 10)
  expect_length(dmsa_bfdr_select(rep(.9, 20), .05), 0)
})

test_that("the cascade finds a planted unit and stops elsewhere", {
  tree <- mk_tree()
  rows <- which(tree$system == 1 & tree$module == 1 & tree$gene <= 2)
  d <- mk_data(tree, rows, h = 0.5)
  cs <- dmsa_cascade(d$b, d$se, d$al, tree, mk_nulls(),
                     engine = c(system = "freq", module = "freq",
                                gene = "freq", probe = "freq"))
  expect_s3_class(cs, "dmsa_cascade")
  sys <- cs$tables$system
  expect_true(sys$selected[sys$unit == "1"])
  expect_equal(sum(sys$selected), 1L)
  expect_true(length(cs$selected_probes) > 0)
  expect_true(all(cs$selected_probes %in% rows))
})

test_that("nothing is selected under a flat null", {
  tree <- mk_tree()
  d <- mk_data(tree, integer(0))
  cs <- dmsa_cascade(d$b, d$se, d$al, tree, mk_nulls(),
                     engine = c(system = "freq", module = "freq",
                                gene = "freq", probe = "freq"))
  expect_equal(sum(cs$tables$system$selected), 0L)
  expect_length(cs$selected_probes, 0)
})

test_that("the gate gates: no descent without a rejected parent", {
  tree <- mk_tree()
  d <- mk_data(tree, integer(0))
  cs <- dmsa_cascade(d$b, d$se, d$al, tree, mk_nulls(),
                     engine = c(system = "freq", module = "freq",
                                gene = "freq", probe = "freq"))
  expect_true(is.null(cs$tables$module) || nrow(cs$tables$module) == 0)
})

test_that("the sparse companion is computed and can open a gate the dense one misses", {
  tree <- mk_tree()
  ## signal confined to one gene of twelve in system 1 -> dense pool diluted
  rows <- which(tree$system == 1 & tree$module == 1 & tree$gene == 1)
  d <- mk_data(tree, rows, h = 1.2, seed = 7)
  nl <- mk_nulls(sparse_sd = 1.2)
  dense <- dmsa_cascade(d$b, d$se, d$al, tree, nl, gate = "dense",
                        engine = c(system = "freq", module = "freq",
                                   gene = "freq", probe = "freq"))
  both  <- dmsa_cascade(d$b, d$se, d$al, tree, nl, gate = "both",
                        engine = c(system = "freq", module = "freq",
                                   gene = "freq", probe = "freq"))
  expect_true(all(!is.na(both$tables$system$sparse_z)))
  expect_gte(sum(both$tables$system$selected), sum(dense$tables$system$selected))
})

test_that("engine defaults are frequentist everywhere (EB is opt-in)", {
  tree <- mk_tree()
  d <- mk_data(tree, integer(0))
  cs <- dmsa_cascade(d$b, d$se, d$al, tree, mk_nulls())
  expect_equal(unname(cs$engine[c("system", "module")]), c("freq", "freq"))
  expect_equal(unname(cs$engine[["gene"]]), "freq")
  expect_equal(unname(cs$engine[["probe"]]), "freq")
})

test_that("an EB level with too few units falls back rather than guessing", {
  tree <- mk_tree(nsys = 3, nmod = 2, ngene = 2, nprobe = 4)
  rows <- which(tree$system == 1 & tree$module == 1)
  d <- mk_data(tree, rows, h = 0.6)
  cs <- dmsa_cascade(d$b, d$se, d$al, tree, mk_nulls())
  expect_s3_class(cs, "dmsa_cascade")           # ran despite tiny EB levels
  expect_true(!is.null(cs$tables$gene))
})

test_that("input mismatches are caught", {
  tree <- mk_tree(); d <- mk_data(tree, integer(0))
  expect_error(dmsa_cascade(d$b[-1], d$se, d$al, tree, mk_nulls()),
               "one entry per probe")
})

test_that("dmsa_cascade_null returns a null per level plus sparse companions", {
  tree <- mk_tree(nsys = 3, nmod = 2, ngene = 2, nprobe = 4)
  set.seed(3); n <- 60; P <- nrow(tree)
  M <- matrix(rnorm(n * P), n, P); x <- rnorm(n)
  al <- data.frame(probe = seq_len(P), gene = tree$gene, s = 1L, p_s_plus = .95)
  nl <- dmsa_cascade_null(rep(0, P), rep(1, P), al, tree, x, M, B = 12L)
  expect_true(all(c("system", "module", "gene") %in% names(nl)))
  expect_true("system_sparse" %in% names(nl))
  expect_true(all(vapply(nl, length, integer(1)) > 0))
})

test_that("a unit with no usable probe does not poison the whole level", {
  set.seed(77); n <- 200; K <- 40
  gene <- rep(paste0("g", 1:8), each = 5)
  tree <- data.frame(system = "S", gene = gene, stringsAsFactors = FALSE)
  d <- rep(c(1, -1), length.out = K)
  ## two genes get NO direction call at all -> their pooled z is NA, and before
  ## the fix that single NA turned every gene's calibrated p into NA and the
  ## level silently selected nothing
  d[gene %in% c("g7", "g8")] <- NA
  al <- dmsa_align(data.frame(cpg = paste0("p", 1:K), d = d,
                              p_plus = ifelse(is.na(d), NA, ifelse(d > 0, .9, .1))),
                   genes = gene, level = "gene")
  x <- rnorm(n)
  M <- matrix(rnorm(n * K), n, K)
  M[, 1:5] <- M[, 1:5] + outer(x, rep(1.2, 5)) * d[1]
  X <- cbind(1, x); XtXi <- solve(crossprod(X))
  bh <- XtXi %*% (t(X) %*% M); r <- M - X %*% bh
  se <- sqrt(colSums(r^2) / (n - 2) * XtXi[2, 2]); b <- bh[2, ]
  nl <- dmsa_cascade_null(b, se, al, tree, exposure = x, M = M, B = 99)
  expect_true(any(!is.finite(nl$gene)))          # the null really does carry NA
  cs <- dmsa_cascade(b, se, al, tree, nl, q = .05)
  gt <- cs$tables$gene
  expect_true(any(is.finite(gt$stat)))           # testable genes get p-values
  expect_true(all(is.na(gt$stat[gt$unit %in% c("S|g7", "S|g8")])))
  expect_true(any(gt$selected))                  # and the real gene is found
})

test_that("BH at a level uses only the testable units", {
  expect_equal(dmsa:::.bh_sel(c(0.01, NA, 0.9), .05), 1L)
  expect_equal(dmsa:::.bh_sel(c(NA, NA), .05), integer(0))
  expect_equal(dmsa:::.bh_sel(c(0.5, 0.6), .05), integer(0))
})

test_that("calibrated p drops NA null entries instead of propagating them", {
  expect_equal(dmsa:::.cal_p(3, c(1, 2, NaN)), 1 / 3)
  expect_true(is.na(dmsa:::.cal_p(NA, c(1, 2))))
  expect_true(is.na(dmsa:::.cal_p(3, c(NA, NA))))
})
