wcsv <- function(d) { f <- tempfile(fileext = ".csv")
  utils::write.csv(d, f, row.names = FALSE); f }

test_that("the template round-trips into a usable bundle", {
  f <- tempfile(fileext = ".csv")
  expect_message(dmsa_reference_template(f), "template written")
  r <- dmsa_reference_csv(f, quiet = TRUE)
  expect_s3_class(r, "dmsa_reference")
  expect_equal(sort(unique(r$systems$system_id)), c("HPA axis", "Oxytocin"))
  expect_true(any(!is.na(r$systems$module_id)))
  expect_equal(r$anchor_method, "user")
  expect_setequal(r$anchors$gene, c("CRH", "OXT"))
  expect_true(any(r$polarity$w_g < 0))
})

test_that("column names are matched loosely", {
  d <- data.frame(`System ID` = c("S1","S1","S1"), SubSystem = c("a","a","b"),
                  Gene_Symbol = c("A","B","C"), WG = c(1,-1,0),
                  check.names = FALSE, stringsAsFactors = FALSE)
  r <- suppressWarnings(dmsa_reference_csv(wcsv(d), quiet = TRUE))
  expect_equal(unique(r$systems$system_id), "S1")
  expect_equal(nrow(r$polarity), 3L)
  expect_true(any(!is.na(r$systems$module_id)))
})

test_that("a missing system or gene column is an actionable error", {
  d <- data.frame(pathway = "S1", w_g = 1, stringsAsFactors = FALSE)
  expect_error(dmsa_reference_csv(wcsv(d)), "no column found for the gene")
  d2 <- data.frame(gene = "A", w_g = 1, stringsAsFactors = FALSE)
  expect_error(dmsa_reference_csv(wcsv(d2)), "no column found for the system")
})

test_that("no w_g column means gene-level only, and it says so", {
  d <- data.frame(system = rep("S1", 4), gene = c("A","B","C","D"),
                  stringsAsFactors = FALSE)
  expect_message(r <- dmsa_reference_csv(wcsv(d), quiet = FALSE),
                 "only GENE-level DMSA is valid")
  expect_null(r$polarity)
  expect_equal(r$anchor_method, "none")
})

test_that("modules are namespaced by system so ids cannot collide", {
  d <- data.frame(system = c("S1","S1","S2","S2"), module = c("a","a","a","a"),
                  gene = c("A","B","C","D"), stringsAsFactors = FALSE)
  r <- dmsa_reference_csv(wcsv(d), min_genes = 2L, quiet = TRUE)
  expect_equal(length(unique(r$systems$module_id)), 2L)
  expect_true(all(grepl("::", r$systems$module_id)))
})

test_that("a non-numeric or out-of-range w_g is refused with the offending value", {
  d <- data.frame(system = rep("S1", 3), gene = c("A","B","C"),
                  w_g = c("1", "up", "0"), stringsAsFactors = FALSE)
  expect_error(dmsa_reference_csv(wcsv(d)), "not numeric")
  expect_error(dmsa_reference_csv(wcsv(d)), "'up'")
  d2 <- data.frame(system = rep("S1", 3), gene = c("A","B","C"),
                   w_g = c(1, 5, 0), stringsAsFactors = FALSE)
  expect_error(dmsa_reference_csv(wcsv(d2)), "\\[-1, 1\\]")
})

test_that("anchors come from the anchor column, else role, else are inferred", {
  base <- data.frame(system = rep("S1", 3), gene = c("A","B","C"),
                     w_g = c(1, .8, -1), stringsAsFactors = FALSE)
  d1 <- transform(base, anchor = c(FALSE, TRUE, FALSE))
  r1 <- dmsa_reference_csv(wcsv(d1), quiet = TRUE)
  expect_equal(r1$anchors$gene, "B"); expect_equal(r1$anchor_method, "user")

  d2 <- transform(base, role = c("readout", "driver", "brake"))
  r2 <- dmsa_reference_csv(wcsv(d2), quiet = TRUE)
  ## provenance is "user": these anchors are derived from the USER's role
  ## column, not from package curation ("curated" was a misstatement)
  expect_equal(r2$anchors$gene, "B"); expect_equal(r2$anchor_method, "user")
  ## and role matching is EXACT - "brake-of-driver" must never anchor
  d2b <- transform(base, role = c("driver", "brake-of-driver", "readout"))
  r2b <- dmsa_reference_csv(wcsv(d2b), quiet = TRUE)
  expect_equal(r2b$anchors$gene, "A")

  expect_warning(r3 <- dmsa_reference_csv(wcsv(base), quiet = TRUE),
                 "anchors inferred")
  expect_equal(r3$anchors$gene, "A")            # the most positive w_g
})

test_that("an all-negative system yields no inferred anchor rather than a wrong one", {
  d <- data.frame(system = rep("S1", 3), gene = c("A","B","C"),
                  w_g = c(-1, -.5, -.2), stringsAsFactors = FALSE)
  ## the warning is correct behaviour - an all-negative system leaves the
  ## system-level sign undefined - so capture it rather than letting it surface
  ## as an uncaptured warning in the suite
  expect_warning(r <- dmsa_reference_csv(wcsv(d), quiet = TRUE),
                 "anchor_method")
  expect_equal(r$anchor_method, "none")
  expect_null(r$anchors)
})

test_that("tiny systems are dropped and reported", {
  d <- data.frame(system = c("Big","Big","Big","Tiny","Tiny"),
                  gene = c("A","B","C","D","E"), stringsAsFactors = FALSE)
  expect_message(r <- dmsa_reference_csv(wcsv(d), min_genes = 3L, quiet = FALSE),
                 "Tiny")
  expect_equal(unique(r$systems$system_id), "Big")
  r2 <- dmsa_reference_csv(wcsv(d), min_genes = 1L, quiet = TRUE)
  expect_setequal(unique(r2$systems$system_id), c("Big", "Tiny"))
})

test_that("blank rows are dropped, not silently kept", {
  d <- data.frame(system = c("S1","S1","","S1"), gene = c("A","B","C",""),
                  stringsAsFactors = FALSE)
  expect_message(dmsa_reference_csv(wcsv(d), min_genes = 1L, quiet = TRUE),
                 "dropping 2 row")
})

test_that("a csv bundle drives dmsa_tree and dmsa_align", {
  f <- tempfile(fileext = ".csv"); dmsa_reference_template(f)
  r <- dmsa_reference_csv(f, quiet = TRUE)
  genes <- c("CRH", "NR3C1", "POMC", "NOTINSET")
  tr <- dmsa_tree(r, genes, system_id = "HPA axis")
  expect_equal(names(tr), c("system", "module", "gene"))
  expect_true(is.na(tr$system[4]))
  pol <- dmsa_polarity_for(r, "HPA axis")
  al <- dmsa_align(data.frame(cpg = paste0("cg", 1:3), d = c(1, 1, -1),
                              p_plus = c(.9, .9, .1)),
                   genes = c("CRH", "NR3C1", "POMC"), level = "system",
                   polarity = pol, missing_polarity = "zero")
  expect_equal(al$w_g, c(1, -1, 1))
  expect_equal(al$s, c(1, -1, -1))
})

test_that("dmsa_reference_write round-trips through dmsa_reference_read", {
  f <- tempfile(fileext = ".csv"); dmsa_reference_template(f)
  r <- dmsa_reference_csv(f, quiet = TRUE)
  d <- file.path(tempdir(), paste0("bundle", sample.int(1e6, 1)))
  dmsa_reference_write(r, d)
  expect_true(all(file.exists(file.path(d, c("systems.csv", "polarity.csv",
                                             "anchors.csv", "manifest.dcf")))))
  r2 <- dmsa_reference_read(d)
  expect_equal(r2$anchor_method, r$anchor_method)
  expect_equal(nrow(r2$systems), nrow(r$systems))
  expect_equal(sort(r2$polarity$w_g), sort(r$polarity$w_g))
})
