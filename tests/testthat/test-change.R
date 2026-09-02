# dmsa_change (the 2-wave identity) and the mod2 three-way through the frame.

.sim_change <- function(n = 160, seed = 5) {
  set.seed(seed)
  n <- 2L * (n %/% 2L)
  G <- 4L; K <- c(6L, 4L, 3L, 5L)
  gene <- rep(paste0("g", 1:G), K); P <- length(gene)
  d <- sample(c(-1, 1), P, TRUE)
  E <- rnorm(n); y <- rnorm(n); cid <- rep(seq_len(n / 2), each = 2)
  M0 <- matrix(rnorm(n * P), n); M1 <- M0 * .6 + matrix(rnorm(n * P), n) * .8
  ## plant: for g1, the CHANGE tracks y x E through the alignment
  idx <- which(gene == "g1")
  for (j in idx) M1[, j] <- M1[, j] + .35 * y * E * d[j]
  al <- dmsa_align(data.frame(cpg = sprintf("cg%04d", 1:P), d = d,
                              p_plus = ifelse(d > 0, .9, .1)),
                   genes = gene, level = "gene")
  list(M0 = M0, M1 = M1, gene = gene, al = al,
       dat = data.frame(y = y, E = E, m2 = rnorm(n), cv = rnorm(n)),
       cid = cid)
}

test_that("dmsa_change finds the planted dS x E unit and adjusts monotonely", {
  s <- .sim_change()
  r <- suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, outcome = "y", exposure = "E",
                units = s$gene, alignment = s$al, covariates = "cv",
                block = s$cid, B = 199, seed = 1))
  expect_s3_class(r, "dmsa_change")
  expect_true(all(is.finite(r$p)))
  g1 <- r[r$unit == "g1", ]
  expect_true(g1$p_adj < .05)
  expect_true(all(r$p_adj[r$unit != "g1"] > g1$p_adj))
  expect_true(all(diff(r$p_adj[order(r$p)]) >= -1e-12))  # step-down monotone
  expect_output(print(r), "2-wave identity")
})

test_that("dmsa_change runs the outcome role and the mod2 three-way", {
  s <- .sim_change()
  r2 <- suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, exposure = "E", role = "outcome",
                units = s$gene, alignment = s$al, covariates = "cv",
                block = s$cid, B = 99, seed = 1))
  expect_true(all(is.finite(r2$p)))
  r3 <- suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, outcome = "y", exposure = "E", mod2 = "m2",
                units = s$gene, alignment = s$al, covariates = "cv",
                block = s$cid, B = 99, seed = 1))
  expect_match(attr(r3, "term"), "mod2")
  expect_true(all(is.finite(r3$p)))
  expect_error(suppressMessages(
    dmsa_change(s$M0, s$M1, s$dat, outcome = "y", exposure = "nope",
                units = s$gene, alignment = s$al, B = 99)), "exposure")
})

test_that("the frame's mod2 three-way detects a planted interaction", {
  set.seed(3)
  n <- 160
  map <- data.frame(gene = rep(c("G1", "G2", "G3"), each = 4))
  map$probe <- sprintf("cg%06d", 1:12)
  map$column <- paste0(map$probe, "_", map$gene)
  map$system_id <- 1L; map$system <- "SimSys"
  map$best_direction <- rep(c(-1, 1), 6)
  map$p_plus <- ifelse(map$best_direction > 0, .9, .1)
  map$best_tier <- "A"; map$smr_tier <- ""
  d <- data.frame(out1 = rnorm(n), m1 = rnorm(n), m2 = rnorm(n),
                  cov1 = rnorm(n), cID = rep(1:(n / 2), each = 2))
  y3 <- scale(d$m1)[, 1] * scale(d$m2)[, 1]
  for (i in seq_len(nrow(map))) {
    base <- rnorm(n)
    if (map$gene[i] == "G1")
      base <- base + .5 * d$out1 * y3 * sign(map$best_direction[i])
    d[[map$column[i]]] <- stats::plogis(base)
  }
  fr <- dmsa_frame(d, map = map, outcome = "out1", covariates = "cov1",
                   blocks = "cID", moderation = TRUE, mod = "m1",
                   mod2 = "m2", B = 199, seed = 1,
                   outdir = file.path(tempdir(), "mod2out"))
  r <- dmsa_report(fr)
  M <- r$moderation
  g1 <- M[M$level == "gene" & M$unit == "G1", ]
  expect_true(nrow(g1) == 1 && g1$p_composite_adj < .05)
  expect_true(all(M[M$level == "gene" & M$unit != "G1", "p_composite"] > .05))
})
