mk_ref <- function(with_modules = TRUE, with_polarity = TRUE) {
  s <- data.frame(
    gene = c("A", "B", "C", "D", "E"),
    system_id = c("S1", "S1", "S1", "S2", "S2"),
    system = c("sys one", "sys one", "sys one", "sys two", "sys two"),
    stringsAsFactors = FALSE)
  if (with_modules) {
    s$module_id <- c("S1.a", "S1.a", "S1.b", "S2.a", "S2.a")
    s$module <- c("m a", "m a", "m b", "m a", "m a")
  }
  p <- if (with_polarity) data.frame(
    gene = c("A", "B", "C", "D", "E"),
    system_id = c("S1", "S1", "S1", "S2", "S2"),
    w_g = c(1, -1, 0, 1, 0.4),
    role = c("driver", "brake", "readout", "driver", "modulator"),
    stringsAsFactors = FALSE) else NULL
  a <- data.frame(system_id = c("S1", "S2"), gene = c("A", "D"),
                  stringsAsFactors = FALSE)
  dmsa_reference(s, p, a, anchor_method = "curated", name = "toy")
}

test_that("a bundle validates its inputs", {
  r <- mk_ref()
  expect_s3_class(r, "dmsa_reference")
  expect_equal(r$anchor_method, "curated")
  bad <- data.frame(gene = "A", system = "x", stringsAsFactors = FALSE)
  expect_error(dmsa_reference(bad), "system_id")
  s <- mk_ref()$systems
  expect_error(dmsa_reference(s, data.frame(gene = "A", system_id = "S1", w_g = 5)),
               "\\[-1, 1\\]")
  expect_error(dmsa_reference(s, data.frame(gene = c("A","A"), system_id = c("S1","S1"),
                                            w_g = c(1, -1))),
               "more than one w_g")
})

test_that("polarity rows outside the system map are flagged", {
  s <- mk_ref()$systems
  expect_warning(
    dmsa_reference(s, data.frame(gene = "ZZZ", system_id = "S1", w_g = 1)),
    "absent from")
})

test_that("continuous and signed polarity are both accepted and distinguished", {
  r <- mk_ref()
  expect_true(any(!r$polarity$w_g %in% c(-1, 0, 1)))
  out <- capture.output(print(r))
  expect_true(any(grepl("continuous", out)))
  r2 <- mk_ref(); r2$polarity$w_g <- c(1, -1, 0, 1, 0)
  out2 <- capture.output(print(r2))
  expect_true(any(grepl("signed", out2)))
})

test_that("a bundle with no polarity says gene-level only", {
  r <- dmsa_reference(mk_ref()$systems, NULL, anchor_method = "none")
  out <- capture.output(print(r))
  expect_true(any(grepl("gene-level analysis only", out)))
})

test_that("anchor_method='none' with signed polarity warns", {
  s <- mk_ref()$systems; p <- mk_ref()$polarity
  expect_warning(dmsa_reference(s, p, anchor_method = "none"), "undefined")
})

test_that("dmsa_tree builds the columns dmsa_cascade expects", {
  r <- mk_ref()
  probes_genes <- c("A", "A", "B", "C", "ZZZ")
  tr <- dmsa_tree(r, probes_genes, system_id = "S1")
  expect_equal(names(tr), c("system", "module", "gene"))
  expect_equal(nrow(tr), 5L)
  expect_equal(tr$system[1:4], rep("S1", 4))
  expect_equal(tr$module[1:4], c("S1.a", "S1.a", "S1.a", "S1.b"))
  expect_true(is.na(tr$system[5]))          # unannotated probe -> flat arm
})

test_that("dmsa_tree omits the module level when the bundle has none", {
  r <- mk_ref(with_modules = FALSE)
  tr <- dmsa_tree(r, c("A", "B"), system_id = "S1")
  expect_equal(names(tr), c("system", "gene"))
})

test_that("a gene in several systems warns unless system_id is given", {
  s <- data.frame(gene = c("A", "A"), system_id = c("S1", "S2"),
                  system = c("one", "two"), stringsAsFactors = FALSE)
  r <- dmsa_reference(s, anchor_method = "none")
  expect_warning(dmsa_tree(r, c("A")), "several systems")
  expect_silent(dmsa_tree(r, c("A"), system_id = "S1"))
})

test_that("dmsa_polarity_for returns what dmsa_align consumes", {
  r <- mk_ref()
  p <- dmsa_polarity_for(r, "S1")
  expect_equal(names(p), c("gene", "w_g"))
  expect_equal(nrow(p), 3L)
  expect_null(dmsa_polarity_for(dmsa_reference(mk_ref()$systems,
                                               anchor_method = "none"), "S1"))
  expect_null(dmsa_polarity_for(r, "NOPE"))
})

test_that("a reference polarity drives dmsa_align end to end", {
  r <- mk_ref()
  dir <- data.frame(cpg = paste0("cg", 1:3), d = c(1, -1, 1), p_plus = c(.9, .1, .9))
  al <- dmsa_align(dir, genes = c("A", "B", "C"), level = "system",
                   polarity = dmsa_polarity_for(r, "S1"),
                   missing_polarity = "zero")
  expect_equal(al$w_g, c(1, -1, 0))
  expect_equal(al$s, c(1, 1, 0))            # d * w_g
  expect_equal(sum(al$usable), 2L)
})

test_that("alpha_reference wraps the bundled tables and finds its anchors", {
  expect_warning(r <- alpha_reference(), "never be used")
  expect_s3_class(r, "dmsa_reference")
  expect_equal(r$anchor_method, "curated")
  expect_gt(nrow(r$systems), 400)
  expect_gt(nrow(r$polarity), 100)
  expect_gt(nrow(r$anchors), 10)
  expect_match(r$notes, "DRAFT")
})

test_that("alpha_reference accepts an added module layer", {
  mods <- data.frame(system_id = "2", gene = c("CRH", "NR3C1"),
                     module_id = c("H", "F"), module = c("hypothalamic", "feedback"),
                     stringsAsFactors = FALSE)
  expect_warning(r <- alpha_reference(modules = mods), "never be used")
  s <- r$systems
  expect_true(any(!is.na(s$module_id)))
  expect_equal(s$module_id[s$gene == "CRH" & s$system_id == "2"], "H")
})

test_that("a bundle round-trips through CSV", {
  d <- file.path(tempdir(), "refbundle"); dir.create(d, showWarnings = FALSE)
  r <- mk_ref()
  utils::write.csv(r$systems, file.path(d, "systems.csv"), row.names = FALSE)
  utils::write.csv(r$polarity, file.path(d, "polarity.csv"), row.names = FALSE)
  utils::write.csv(r$anchors, file.path(d, "anchors.csv"), row.names = FALSE)
  writeLines(c("name: toy-from-disk", "version: 1", "anchor_method: graph_sink"),
             file.path(d, "manifest.dcf"))
  r2 <- dmsa_reference_read(d)
  expect_equal(r2$name, "toy-from-disk")
  expect_equal(r2$anchor_method, "graph_sink")
  expect_equal(nrow(r2$systems), nrow(r$systems))
  expect_equal(r2$polarity$w_g, r$polarity$w_g)
  expect_error(dmsa_reference_read(tempfile()), "no systems.csv")
})
