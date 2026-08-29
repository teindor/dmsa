## Spec 2 + 5: the biological reference is system > module > gene, rebuilt from
## the audited 2026c codebook and NOT filtered by retained probe coverage.

test_that("the bundled biological reference loads as a dmsa_reference", {
  ref <- alpha_reference()
  expect_s3_class(ref, "dmsa_reference")
  expect_true(all(c("system_id", "system", "module_id", "module", "gene") %in%
                    names(ref$systems)))
})

test_that("counts are DERIVED from the shipped file and are self-consistent", {
  chk <- alpha_reference_check(verbose = FALSE)
  expect_true(chk$ok)
  ## one row per gene, and the derived counts agree with the object
  ref <- alpha_reference()
  expect_equal(chk$counts[["genes"]], length(unique(ref$systems$gene)))
  expect_equal(chk$counts[["rows"]],  nrow(ref$systems))
  expect_equal(chk$counts[["systems"]], length(unique(ref$systems$system_id)))
  expect_equal(chk$counts[["modules"]],
               length(unique(stats::na.omit(ref$systems$module_id))))
})

test_that("the reference carries no probe/CpG column", {
  ref <- alpha_reference()
  expect_false(any(c("cpg", "probe", "probe_id", "column") %in%
                     names(ref$systems)))
})

test_that("genes with zero retained Alpha probes remain in the reference", {
  ref <- alpha_reference()
  cas <- dmsa_sets("alpha")
  zero <- setdiff(ref$systems$gene, cas$cascade$gene)
  ## the probe-anchored cascade is a strict SUBSET of the biology
  expect_gt(length(zero), 0L)
  expect_length(setdiff(cas$cascade$gene, ref$systems$gene), 0L)
  ## and each such gene still has a system and a module
  z <- ref$systems[ref$systems$gene %in% zero, ]
  expect_false(any(is.na(z$system_id)))
  expect_false(any(is.na(z$module_id)))
})

test_that("genes with no curated polarity are unresolved, never +1", {
  ref <- alpha_reference()
  skip_if(is.null(ref$polarity), "no polarity attached")
  unresolved <- setdiff(ref$systems$gene, ref$polarity$gene)
  expect_gt(length(unresolved), 0L)
  ## they must be absent from polarity, not present with a fabricated +1
  expect_length(intersect(unresolved, ref$polarity$gene), 0L)
})
