test_that("the bundled map is present, versioned, and self-describing", {
  m <- dmsa_directions()
  meta <- attr(m, "dmsa_direction_map")
  expect_true(nrow(m) > 4e5)
  expect_setequal(names(m), c("probe", "gene", "d", "p_plus", "tier"))
  ## an abstention is not a direction: nothing unusable should have shipped
  expect_true(all(m$d %in% c(-1L, 1L)))
  expect_false(any(is.na(m$p_plus)))
  expect_true(all(m$p_plus >= 0 & m$p_plus <= 1))
  ## provenance is the point, not decoration
  expect_true(all(c("tissue", "version", "source", "source_doi") %in%
                    names(meta)))
  expect_identical(meta$tissue, "blood")
})

test_that("lookup by probe, by gene, and by pinned pair each behave", {
  by_gene <- dmsa_directions(genes = "AVP")
  expect_true(nrow(by_gene) > 0)
  expect_true(all(by_gene$gene == "AVP"))

  by_probe <- dmsa_directions(c("cg00052046", "not_a_real_probe"))
  expect_true(all(by_probe$probe == "cg00052046"))
  expect_equal(nrow(dmsa_directions("not_a_real_probe")), 0)

  ## same length probes + genes selects those specific pairs, not the union
  pinned <- dmsa_directions(c("cg00052046", "cg00176879"), c("AVP", "AVP"))
  expect_equal(nrow(pinned), 2)
  expect_true(all(pinned$gene == "AVP"))
})

test_that("dmsa_align takes probe IDs alone and records which map it used", {
  al <- suppressMessages(
    dmsa_align(c("cg00052046", "cg00176879", "cg00308631")))
  expect_s3_class(al, "dmsa_alignment")
  expect_true(nrow(al) >= 3)
  expect_true(all(al$d %in% c(-1, 1)))
  expect_match(attr(al, "direction_map"), "bundled blood map v")

  ## pinning genes restricts to those pairs
  pinned <- suppressMessages(
    dmsa_align(c("cg00052046", "cg00176879"), genes = c("AVP", "AVP")))
  expect_equal(nrow(pinned), 2)
  expect_true(all(pinned$gene == "AVP"))

  ## a probe set with no calls is an error that names the alternative
  expect_error(suppressMessages(dmsa_align(c("nope1", "nope2"))),
               "cpgdirection")
})

test_that("a user-supplied table is untouched by the bundled map", {
  dd <- data.frame(cpg = c("p1", "p2"), d = c(-1, -1), p_plus = c(.2, .2))
  al <- dmsa_align(dd, genes = c("CRH", "NR3C1"), level = "system",
                   polarity = data.frame(gene = c("CRH", "NR3C1"),
                                         w_g = c(1, -1)))
  expect_equal(al$s, c(-1, 1))
  ## and it says so, rather than claiming provenance it does not have
  expect_identical(attr(al, "direction_map"), "user-supplied direction table")
  ## genes remains required for a table
  expect_error(dmsa_align(dd), "`genes` is required")
})
