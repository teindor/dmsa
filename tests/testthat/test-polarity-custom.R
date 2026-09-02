## Spec 67: a user polarity table, end to end.
##
## Validation (malformed tables error) and the engine-level merge (E4) were
## pinned before. These run a user cascade + user polarity THROUGH
## dmsa_report() and pin the properties a user relies on:
##   * antisymmetry - flipping every gene's polarity mirrors the system level
##     exactly (a RULE ZERO invariant: same magnitude, opposite direction,
##     identical p under the same seed)
##   * locality - polarity touches ONLY the system level; gene and probe
##     levels align on d alone and are identical between runs
##   * the frame carries the user's table gene for gene, and says so
##   * missing_polarity policies at the engine level do what they say

.pc_sim <- function(n = 120, seed = 42) {
  set.seed(seed)
  map <- expand.grid(gene = paste0("G", 1:6), k = 1:4, stringsAsFactors = FALSE)
  map$system_id <- ifelse(map$gene %in% paste0("G", 1:3), 1L, 2L)
  map$system <- ifelse(map$system_id == 1, "Sim system one", "Sim system two")
  map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
  map$column <- paste0(map$probe, "_", map$gene)
  map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
  map$p_plus <- ifelse(map$best_direction > 0, .9, .1)
  map$best_tier <- "A"; map$smr_tier <- ""
  d <- data.frame(out1 = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cov2 = stats::rnorm(n), cID = rep(seq_len(n / 2), each = 2))
  for (i in seq_len(nrow(map))) {
    sig <- if (map$gene[i] == "G1")
      0.35 * d$out1 * sign(map$best_direction[i]) else 0
    d[[map$column[i]]] <- stats::plogis(stats::rnorm(n, sd = 1) + sig)
  }
  cas <- unique(map[c("system_id", "system", "gene")])
  cas$module_id <- paste0(cas$system_id, ".1"); cas$module <- "M"
  cas$cpg <- map$probe[match(cas$gene, map$gene)]
  list(data = d, map = map, cascade = cas)
}

.pc_report <- function(s, w_g, od = tempfile("dmsa_pol67")) {
  pol <- data.frame(system_id = s$cascade$system_id, gene = s$cascade$gene,
                    w_g = w_g, stringsAsFactors = FALSE)
  st <- dmsa_sets(s$cascade, polarity = pol)
  fr <- suppressMessages(dmsa_frame(
    s$data, map = s$map, sets = st, outcome = "out1",
    covariates = c("cov1", "cov2"), blocks = "cID",
    B = 99, seed = 1, outdir = od, plot_type = "png",
    progress = FALSE, beep = FALSE))
  list(frame = fr,
       report = suppressWarnings(suppressMessages(dmsa_report(fr))))
}

.pc_lvl <- function(r, lv) {
  x <- r$results[r$results$level == lv, , drop = FALSE]
  x[order(x$unit), , drop = FALSE]
}

test_that("spec 67: flipping every gene's polarity mirrors the system level exactly", {
  s <- .pc_sim()
  w <- c(1, -1, 1, 1, 1, -1)
  a <- .pc_report(s, w)
  b <- .pc_report(s, -w)
  sa <- .pc_lvl(a$report, "system"); sb <- .pc_lvl(b$report, "system")
  expect_equal(nrow(sa), 2L); expect_equal(nrow(sb), 2L)
  ## opposite direction, identical evidence
  expect_equal(sb$direction, -sa$direction)
  pcols <- c("p_coherence", "p_composite", "p_diffuse", "p_coherence_adj",
             "p_composite_adj", "p_diffuse_adj", "p_omnibus", "p_union_exact")
  for (pc in pcols) expect_equal(sb[[pc]], sa[[pc]], info = pc)
  expect_equal(sb$concordance, sa$concordance)
  expect_equal(sb$n_probes, sa$n_probes)
  ## the signal system is the one that is named, either way round
  expect_true(sa$selected[sa$unit == "Sim system one"])
  expect_true(sb$selected[sb$unit == "Sim system one"])
})

test_that("spec 67: polarity touches ONLY the system level", {
  s <- .pc_sim()
  w <- c(1, -1, 1, 1, 1, -1)
  a <- .pc_report(s, w)
  w2 <- w; w2[1] <- -1                          # flip G1 alone (the signal gene)
  b <- .pc_report(s, w2)
  lvls <- intersect(c("gene", "module", "probe"), unique(a$report$results$level))
  expect_true("gene" %in% lvls)
  for (lv in lvls) {
    ga <- .pc_lvl(a$report, lv); gb <- .pc_lvl(b$report, lv)
    keep <- intersect(names(ga), c("unit", "direction", "p_coherence",
                                   "p_composite", "p_diffuse", "p_coherence_adj",
                                   "p_composite_adj", "p_diffuse_adj", "selected"))
    expect_equal(gb[keep], ga[keep], info = lv)
  }
  ## and system one, whose only signal gene was flipped, flips with it
  sa <- .pc_lvl(a$report, "system"); sb <- .pc_lvl(b$report, "system")
  expect_equal(sb$direction[sb$unit == "Sim system one"],
               -sa$direction[sa$unit == "Sim system one"])
})

test_that("spec 67: the frame carries the user's polarity gene for gene, and says so", {
  s <- .pc_sim()
  w <- c(1, -1, 1, 1, 1, -1)
  a <- .pc_report(s, w)
  pt <- a$frame$polarity_table
  expect_s3_class(pt, "dmsa_polarity")
  p <- pt$polarity
  expect_equal(p$w_g[match(s$cascade$gene, p$gene)], w)
  expect_identical(pt$name, "user data.frame")
  out <- paste(utils::capture.output(print(a$frame)), collapse = "\n")
  expect_match(out, "polarity", ignore.case = TRUE)
  ## and the system line in summary.md is keyed to the SIGNED score
  sm <- readLines(file.path(a$frame$outdir, "summary.md"))
  expect_true(any(grepl("Sim system one", sm)))
})

test_that("spec 67: missing_polarity policies do what they say (engine level)", {
  probes <- c("cgA1", "cgA2", "cgB1")
  dir <- data.frame(cpg = probes, d = c(1, -1, 1), p_plus = c(.9, .1, .9))
  pol <- data.frame(gene = "GENEA", w_g = 1)          # GENEB has no entry
  ## error: names the gene and how to fix it
  e <- expect_error(
    suppressMessages(dmsa_align(dir, genes = c("GENEA", "GENEA", "GENEB"),
                                level = "system", polarity = pol,
                                missing_polarity = "error")),
    "No polarity")
  expect_match(conditionMessage(e), "GENEB")
  expect_match(conditionMessage(e), "polarity = data.frame", fixed = TRUE)
  ## zero: the gene stays, off-axis, with a warning naming it
  expect_warning(
    z <- suppressMessages(dmsa_align(dir, genes = c("GENEA", "GENEA", "GENEB"),
                                     level = "system", polarity = pol,
                                     missing_polarity = "zero")),
    "GENEB")
  expect_equal(z$w_g[z$gene == "GENEB"], 0)
  expect_equal(z$s[z$gene == "GENEB"], 0)
  expect_false(any(z$usable[z$gene == "GENEB"]))
  expect_true(all(z$usable[z$gene == "GENEA"]))
  ## drop: the gene's rows are gone, the others untouched
  dr <- suppressMessages(dmsa_align(dir, genes = c("GENEA", "GENEA", "GENEB"),
                                    level = "system", polarity = pol,
                                    missing_polarity = "drop"))
  expect_false("GENEB" %in% dr$gene)
  expect_equal(sum(dr$gene == "GENEA"), 2L)
  expect_equal(dr$s[dr$gene == "GENEA"], c(1, -1))
})

test_that("spec 67: a user table naming one gene twice is refused", {
  dir <- data.frame(cpg = c("cgA1", "cgA2"), d = c(1, -1), p_plus = c(.9, .1))
  expect_error(
    suppressMessages(dmsa_align(dir, genes = c("GENEA", "GENEA"), level = "system",
                                polarity = data.frame(gene = c("GENEA", "GENEA"),
                                                      w_g = c(1, -1)))),
    "more than once")
})
