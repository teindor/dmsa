
test_that("the moderated report block states the design, not the main-effect one", {
  MOD <- data.frame(
    unit = c("FKBP5", "NR3C1"), n_probes = c(5L, 7L),
    b = c(-0.171, 0.204), t = c(-2.09, 3.41),
    p_composite = c(.040, .001), p_composite_adj = c(.265, .008),
    family = "HPA", level = "gene", outcome = "anx",
    selected = c(FALSE, TRUE), stringsAsFactors = FALSE)
  cfg <- list(mod = "ACE", mod2 = "", correction = "maxT", B = 1999,
              weighting = "flat", cpg_map = "full", frame_role = "predictor",
              block_cols = "cID")
  txt <- paste(.rp_mod_block(cfg, MOD, 396L), collapse = " ")
  expect_match(txt, "subject-level aligned tone score")
  expect_match(txt, "composite lens")
  expect_match(txt, "do not carry through a product")
  expect_match(txt, "NR3C1")
  ## the survivor sentence must read as a slope-on-a-slope
  expect_match(txt, "changed that slope by \\+0\\.204")

  ## two moderators: the coefficient is the THREE-way term and must not be
  ## described as a change in the simple slope
  cfg2 <- cfg; cfg2$mod2 <- "sex"
  t2 <- paste(.rp_mod_block(cfg2, MOD, 396L), collapse = " ")
  expect_match(t2, "tone score x ACE x sex")
  expect_match(t2, "depends on sex")
  expect_false(grepl("changed that slope", t2))

  ## nothing survives -> the calibrated-null wording, still with the design
  MOD$selected <- FALSE
  t3 <- paste(.rp_mod_block(cfg, MOD, 396L), collapse = " ")
  expect_match(t3, "calibrated null")
  expect_match(t3, "No unit survived")
})

test_that("the moderated Methods paragraph enumerates the model it fitted", {
  fr2 <- list(mod = "ACE", mod2 = "", frame_role = "predictor",
              outcome = c("y1", "y2"), covariates = c("sex_c", "age"),
              chip = "chip_f", block_cols = "cID")
  w <- .rp_mod_model_words(fr2, "y1")
  expect_match(w$model, "the tone score \\+ ACE \\+ the tone score x ACE")
  expect_match(w$model, "both main effects")
  expect_match(w$covs, "sex_c, age, chip_f")
  ## the multi-outcome cross-adjustment must be stated, not left implicit
  expect_match(w$covs, "each was additionally adjusted for the other")
  ## and the block must not be passed off as a random effect
  expect_match(w$covs, "rather than by a random intercept")
  expect_false(grepl("random effect[^s]", w$covs))

  ## three-way: every lower-order term named, three-way flagged as the test
  fr3 <- fr2; fr3$mod2 <- "sex"
  w3 <- .rp_mod_model_words(fr3, "y1")
  for (tm in c("the tone score x ACE", "the tone score x sex", "ACE x sex",
               "the tone score x ACE x sex"))
    expect_true(grepl(tm, w3$model, fixed = TRUE))
  expect_match(w3$model, "all three main effects and all three two-way")
  expect_match(w3$model, "three-way product")
})

test_that("a non-linear run surfaces its shape scan instead of discarding it", {
  ## Regression guard for the 1.8.9 bug: shape_rows was assembled and never
  ## read, so type = "non-linear" returned a linear result silently.
  skip_if_not(exists("alpha_reference", mode = "function"))
  set.seed(4)
  n <- 160
  d <- data.frame(y = rnorm(n), cID = rep(seq_len(n / 2), each = 2),
                  age = rnorm(n))
  for (nm in paste0("cg", 1:6)) d[[nm]] <- stats::runif(n, .2, .8)
  al <- dmsa_align(data.frame(cpg = paste0("cg", 1:6), d = c(1, 1, -1, 1, -1, 1),
                              p_plus = rep(.9, 6)),
                   genes = rep(c("G1", "G2"), each = 3), level = "gene")
  sc <- dmsa_scores(as.matrix(d[, paste0("cg", 1:6)]), al, weighting = "flat")
  expect_true(is.list(sc))
  ## the omnibus helper must expose per-arm p-values, which is what the report
  ## needs and what .report_shape used to drop on the floor
  arms <- data.frame(lin = as.numeric(scale(sc$aligned)),
                     quad = as.numeric(scale(sc$aligned))^2 - 1,
                     thr = as.numeric(scale(sc$aligned) > 0))
  o <- dmsa_model_omnibus(arms, y ~ S + age, d, "S", block = d$cID, B = 99,
                          seed = 1)
  expect_true(all(c("lin", "quad", "thr") %in% names(o$p_by_flavour)))
  expect_true(is.finite(o$p))
})

test_that("a monotone convex relationship is NOT reported as a U", {
  ## The Lind & Mehlum / Simonsohn failure mode, as a regression test. y grows
  ## throughout but with increasing steepness, so the quadratic term is real
  ## while the shape is monotone. A significant squared term must NOT be allowed
  ## to license a "U-shaped" claim.
  set.seed(11)
  n <- 400
  x <- stats::runif(n, 0, 3)                  # note: x >= 0, curve rises always
  y <- 0.5 * x + 0.8 * x^2 + stats::rnorm(n, sd = 1.2)
  d <- data.frame(y = y, S = as.numeric(scale(x)))
  d$.sq <- d$S^2 - mean(d$S^2)
  m <- stats::lm(y ~ S + .sq, d)
  cf <- stats::coef(m); V <- stats::vcov(m); dfr <- stats::df.residual(m)
  ## the squared term IS significant - that is the whole point of the example
  expect_lt(summary(m)$coefficients[".sq", "Pr(>|t|)"], .01)
  ## but the slope has the SAME sign at both ends, so the shape test must fail
  cv <- function(z) { v <- stats::setNames(numeric(length(cf)), names(cf))
                      v["S"] <- 1; v[".sq"] <- 2 * z; v }
  sl <- function(z) sum(cv(z) * cf)
  se <- function(z) { v <- cv(z); sqrt(drop(v %*% V %*% v)) }
  xL <- min(d$S); xH <- max(d$S)
  expect_gt(sl(xL) * sl(xH), 0)               # same sign at both ends
  pu  <- max(stats::pt(sl(xL) / se(xL), dfr), stats::pt(-sl(xH) / se(xH), dfr))
  piu <- max(stats::pt(-sl(xL) / se(xL), dfr), stats::pt(sl(xH) / se(xH), dfr))
  expect_gt(min(pu, piu), .05)                # neither U nor inverted-U
  ## and the turning point sits outside the observed range
  xs <- -cf["S"] / (2 * cf[".sq"])
  expect_true(xs < xL || xs > xH)
})

test_that("a genuine U is detected by the same slope-sign machinery", {
  set.seed(12)
  n <- 400
  x <- stats::runif(n, -3, 3)
  y <- 0.9 * x^2 + stats::rnorm(n, sd = 1.2)   # true minimum at 0, interior
  d <- data.frame(y = y, S = as.numeric(scale(x)))
  d$.sq <- d$S^2 - mean(d$S^2)
  m <- stats::lm(y ~ S + .sq, d)
  cf <- stats::coef(m); V <- stats::vcov(m); dfr <- stats::df.residual(m)
  cv <- function(z) { v <- stats::setNames(numeric(length(cf)), names(cf))
                      v["S"] <- 1; v[".sq"] <- 2 * z; v }
  sl <- function(z) sum(cv(z) * cf)
  se <- function(z) { v <- cv(z); sqrt(drop(v %*% V %*% v)) }
  xL <- min(d$S); xH <- max(d$S)
  expect_lt(sl(xL) * sl(xH), 0)                # opposite signs
  pu <- max(stats::pt(sl(xL) / se(xL), dfr), stats::pt(-sl(xH) / se(xH), dfr))
  expect_lt(pu, .05)
  xs <- -cf["S"] / (2 * cf[".sq"])
  expect_true(xs > xL && xs < xH)              # turning point interior
})

test_that("the unit builder covers the full declared hierarchy", {
  ## Regression guard: .rp_units had no probe branch, so the shape scan and the
  ## moderated path stopped at gene while the design said system > module >
  ## gene > probe.
  expect_true("probe" %in% names(formals(dmsa_frame)))
  b <- body(.rp_units)
  txt <- paste(deparse(b), collapse = " ")
  for (lv in c("system", "gene", "module", "probe"))
    expect_true(grepl(paste0('level == "', lv, '"'), txt, fixed = TRUE))
})

test_that("the moderated non-linear term is resolved from the design, not guessed", {
  ## Regression guard: `S * .mc + .sq * .mc` names the column `.mc:.sq`, not
  ## `.sq:.mc`. Guessing wrong made dmsa_model() throw, and a tryCatch turned
  ## that into a silent column of NAs.
  d <- data.frame(y = rnorm(40), S = rnorm(40), age = rnorm(40))
  d$.mc <- rnorm(40); d$.sq <- d$S^2 - mean(d$S^2)
  ff <- y ~ S * .mc + .sq * .mc + age
  cn <- colnames(model.matrix(ff, model.frame(ff, d)))
  expect_false(".sq:.mc" %in% cn)
  expect_true(".mc:.sq" %in% cn)
  ## the resolver must find it whichever way round R writes it
  parts <- strsplit(cn, ":", fixed = TRUE)
  ok <- vapply(parts, function(z) length(z) == 2L &&
                 setequal(z, c(".sq", ".mc")), logical(1))
  expect_equal(sum(ok), 1L)
  ## and every lower-order term must be present (marginality)
  for (tm in c("S", ".mc", ".sq", "S:.mc")) expect_true(tm %in% cn)
})

# ---------------------------------------------------------------------------
# Title blocks: the figure must never draw above its own top margin.
#
# Regression for a shipped defect. The module figure fixed mar[3] at 4.0 lines
# and drew the first of its wrapped title lines at line 4.2, so a title that
# needed three lines lost line one off the top of the device. What went missing
# was the outcome name - the one part of the title a reader cannot reconstruct
# from what survived. Nothing errored and no text overflowed horizontally, so
# neither .rp_text_fits() nor the eye caught it.
test_that("the title block fits inside the top margin it asks for", {
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 8.4 * 200, height = 5 * 200, res = 200)
  on.exit({
    if (grDevices::dev.cur() > 1L) grDevices::dev.off()
    unlink(f)
  }, add = TRUE)

  labs <- c(
    "X",
    "Attachment anxiety (T1)",
    "Attachment_Anxiety_General_T1",
    paste("Change in maternal attachment anxiety, general relational form,",
          "measured at T1 and again at T4 with the 36-item ECR-R"))
  mars <- list(c(4.2, 17, 4.0, 6.5), c(4.2, 11, 5.2, 5.5))

  for (m in mars) for (lab in labs) {
    h <- dmsa:::.rp_head(
      m,
      sprintf("%s - module level, three lenses, maxT within each system", lab),
      paste("Each block is one family: a module is corrected only against the",
            "other modules of its own system. Bold = survives."))
    e <- environment(h)
    top <- graphics::par("mar")[3]
    ## every title line lands strictly inside the margin that was reserved
    expect_true(max(e$ttl_at) < top,
                info = paste("title overflows top margin:", lab))
    ## and the subtitle stays in the margin rather than dropping into the panel
    expect_true(min(e$sub_at) > 0,
                info = paste("subtitle falls into the panel:", lab))
    ## title above subtitle, never overlapping
    expect_true(min(e$ttl_at) > max(e$sub_at))
  }
})

test_that("locus panels carry the outcome, and survive a label of any length", {
  pr <- data.frame(probe = sprintf("cg%08d", 1:8),
                   b = seq(.01, .03, length.out = 8), se = rep(.008, 8),
                   d = rep(-1, 8), chr = "20",
                   pos = 3084690 + c(0, 40, 70, 130, 160, 220, 250, 255),
                   stringsAsFactors = FALSE)
  for (ctx in c("", "Attachment anxiety (T1)",
                paste(rep("a very long outcome label indeed", 6), collapse = " "))) {
    f <- tempfile(fileext = ".png")
    expect_silent(suppressMessages(
      dmsa_plot_locus(pr, gene = "AVP", context = ctx, file = f)))
    expect_true(file.exists(f) && file.size(f) > 0)
    unlink(f)
  }
})
