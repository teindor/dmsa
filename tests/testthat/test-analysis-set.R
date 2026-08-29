## dmsa_save_analysis_set(): the probes a frame actually tested, one row per
## CpG x gene pair, annotated, as a publication/validation CSV (PI,
## 2026-08-29). Both paths are covered: the pair path carries evidence tiers,
## the bundled/custom-map path carries what it has.

.as_pairs <- function() data.frame(
  cpg_id = c("cg0000001", "cg0000001", "cg0000002"),
  target_gene = c("OXTR", "AVP", "NR3C1"),
  best_direction = c(-1, 1, -1),
  probability_plus1 = c(.05, .9, .08),
  usable = TRUE, abstain_reason = NA_character_,
  best_evidence = c("measured", "smr_weak", "smr_high"),
  direction_tier = c("M", "S3", "S1"), stringsAsFactors = FALSE)

.as_frame_pairs <- function() {
  set.seed(3)
  cols <- c("cg0000001", "cg0000002")
  d <- data.frame(y = rnorm(60), cov1 = rnorm(60),
                  cID = rep(1:30, each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(rnorm(60))
  old <- options(dmsa.pair_table = .as_pairs()); on.exit(options(old))
  dmsa_frame(d, methylation = cols, direction_source = "cpgdirection",
             outcomes = "y", covariates = "cov1", random_effects = "cID",
             chip = FALSE, B = 19, plots = FALSE, tables = FALSE,
             summary = FALSE, progress = FALSE, beep = FALSE,
             outdir = tempfile("dmsa_as"))
}

test_that("the analysis set has one annotated row per tested pair", {
  f <- .as_frame_pairs()
  a <- dmsa_save_analysis_set(f, file = NULL)
  expect_equal(nrow(a), 3L)                     # co-effect kept as two rows
  expect_setequal(a$gene, c("OXTR", "AVP", "NR3C1"))
  expect_true(all(c("probe", "gene", "system", "direction", "p_plus1",
                    "direction_tier", "evidence", "co_effect",
                    "genes_for_this_cpg", "direction_source", "tissue",
                    "reference") %in% names(a)))
  expect_identical(unique(a$direction_source), "cpgdirection")
  expect_identical(unique(a$tissue), "blood")
  ## the co-effect probe is marked and counts both of its genes
  co <- a[a$probe == "cg0000001", ]
  expect_true(all(co$co_effect))
  expect_true(all(co$genes_for_this_cpg == 2L))
})

test_that("the analysis set is written as a readable CSV", {
  f <- .as_frame_pairs()
  fp <- tempfile(fileext = ".csv")
  r <- dmsa_save_analysis_set(f, file = fp)
  expect_true(file.exists(fp))
  back <- utils::read.csv(fp, stringsAsFactors = FALSE)
  expect_equal(nrow(back), 3L)
  expect_identical(back$probe, r$probe)
})

test_that("a bundled/custom-map frame gets an analysis set too", {
  set.seed(4)
  map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
                    probe = c("cg01", "cg02"), column = c("cg01", "cg02"),
                    best_direction = c(-1, 1), p_plus = c(.1, .9))
  dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40),
                    cID = rep(1:20, each = 2),
                    cg01 = stats::plogis(rnorm(40)),
                    cg02 = stats::plogis(rnorm(40)))
  fr <- dmsa_frame(dat, map = map, outcomes = "anx", covariates = "cov1",
                   B = 19, plots = FALSE, tables = FALSE, summary = FALSE,
                   progress = FALSE, beep = FALSE, outdir = tempfile())
  a <- dmsa_save_analysis_set(fr, file = NULL)
  expect_equal(nrow(a), 2L)
  expect_true(all(c("probe", "gene", "system", "direction") %in% names(a)))
  ## no evidence tiers on this path - and no empty placeholder columns either
  expect_false("direction_tier" %in% names(a))
  expect_false("tissue" %in% names(a))
})

test_that("only a dmsa_frame is accepted", {
  expect_error(dmsa_save_analysis_set(list()), "dmsa_frame")
})
