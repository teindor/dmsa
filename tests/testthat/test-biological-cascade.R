## Spec section 4: the BIOLOGICAL cascade is system > module > gene.
## A CpG/probe column is optional provenance, never part of the definition.

bio_cascade <- function()
  data.frame(
    system_id = c(1, 1, 1, 2, 2),
    system    = c(rep("Sys A", 3), rep("Sys B", 2)),
    module_id = c("1.1", "1.1", "1.2", "2.1", "2.1"),
    module    = c("Mod A1", "Mod A1", "Mod A2", "Mod B1", "Mod B1"),
    gene      = c("G1", "G2", "G3", "G4", "G5"),
    stringsAsFactors = FALSE)

test_that("a cascade with no CpG column is accepted", {
  s <- dmsa_sets(bio_cascade(), name = "biological")
  expect_s3_class(s, "dmsa_sets")
  expect_equal(nrow(s$systems), 2L)
  expect_equal(length(unique(s$cascade$gene)), 5L)
  expect_false("cpg" %in% names(s$cascade))
  expect_true(all(is.na(s$systems$n_cpgs)))
})

test_that("dmsa_sets_check passes a biological cascade", {
  chk <- dmsa_sets_check(bio_cascade(), verbose = FALSE)
  expect_true(chk$ok)
  expect_equal(chk$shared_cpgs, 0L)
})

test_that("printing a biological cascade does not fabricate CpG counts", {
  out <- paste(utils::capture.output(print(dmsa_sets(bio_cascade()))),
               collapse = "\n")
  expect_match(out, "no CpG column")
  expect_false(grepl("NA cpg", out, fixed = TRUE))
})

test_that("gene-level selection works without a CpG column", {
  s   <- dmsa_sets(bio_cascade())
  sel <- dmsa_select(s, systems = 1, genes = c("G1", "G2"))
  expect_setequal(unique(sel$cascade$gene), c("G1", "G2"))
})

test_that("asking for probes on a biological cascade is an explicit error", {
  s <- dmsa_sets(bio_cascade())
  expect_error(dmsa_select(s, systems = 1, probes = "cg123"),
               "carries no CpG column")
  ## the default must NOT trip it
  expect_silent(invisible(dmsa_select(s, systems = 1)))
})

test_that("a cascade WITH CpGs is unaffected", {
  a <- dmsa_sets("alpha")
  expect_true("cpg" %in% names(a$cascade))
  expect_equal(nrow(a$systems), 30L)
  expect_true(all(a$systems$n_cpgs > 0))
})
