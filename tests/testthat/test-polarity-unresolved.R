## Spec 40 + E13: unresolved polarity is COUNTED and SAID, and a user
## polarity table on the bundled cascade MERGES instead of replacing.
##
## A gene whose activation sign is unresolved is weighted 0 in the system
## score (never assumed +1, spec 39/40) - which means it silently drops out
## of the system tone. That is now stated at the frame, in summary.md, in
## tables/polarity_audit.csv and in dmsa_coverage(). And dmsa_sets("alpha",
## polarity = <table>) now behaves like dmsa_align(polarity = ): user rows
## override gene for gene, the curation fills the rest (E13).

.pu_sim <- function(n = 120, seed = 42) {
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

test_that("E13: a user polarity table on the bundled cascade MERGES gene for gene", {
  bnd <- dmsa_polarity("alpha")$polarity
  key <- paste(bnd$system_id, bnd$gene)
  ## override one curated gene: table size unchanged, that row is the user's
  expect_message(
    st <- dmsa_sets("alpha", polarity = data.frame(system_id = "2", gene = "NR3C1",
                                                    w_g = 1)),
    "polarity merged")
  p <- st$polarity$polarity
  expect_equal(nrow(p), nrow(bnd))
  i <- p$system_id == "2" & p$gene == "NR3C1"
  expect_equal(sum(i), 1L)
  expect_equal(p$w_g[i], 1)
  expect_identical(p$grade[i], "user")
  ## every OTHER gene keeps its bundled sign
  j <- !(paste(p$system_id, p$gene) %in% "2 NR3C1")
  m <- match(paste(p$system_id, p$gene)[j], key)
  expect_equal(p$w_g[j], bnd$w_g[m])
  expect_match(st$polarity$name, "user overrides")
  ## a gene the curation does not know is ADDED, not dropped
  st2 <- suppressMessages(dmsa_sets("alpha", polarity = data.frame(
    system_id = "2", gene = "NOTAGENE", w_g = -1)))
  p2 <- st2$polarity$polarity
  expect_equal(nrow(p2), nrow(bnd) + 1L)
  expect_equal(p2$w_g[p2$gene == "NOTAGENE"], -1)
})

test_that("E13: a USER cascade keeps its own table as the polarity (no merge)", {
  s <- .pu_sim()
  pol <- data.frame(system_id = s$cascade$system_id, gene = s$cascade$gene,
                    w_g = c(1, -1, 1, 1, 1, -1))
  expect_no_message(st <- dmsa_sets(s$cascade, polarity = pol),
                    message = "polarity merged")
  expect_equal(nrow(st$polarity$polarity), 6L)
})

test_that("spec 40: unresolved polarity is counted at the frame, in summary.md, in the audit table", {
  s <- .pu_sim()
  pol <- data.frame(system_id = s$cascade$system_id, gene = s$cascade$gene,
                    w_g = c(1, -1, 1, 1, 1, -1))
  pol <- pol[pol$gene != "G3", ]                      # G3 has no sign
  st <- dmsa_sets(s$cascade, polarity = pol)
  od <- tempfile("dmsa_pol40")
  fr <- suppressMessages(dmsa_frame(
    s$data, map = s$map, sets = st, outcome = "out1",
    covariates = c("cov1", "cov2"), blocks = "cID",
    B = 99, seed = 1, outdir = od, plot_type = "png",
    progress = FALSE, beep = FALSE))
  ## the frame says it, per system, naming the gene
  out <- paste(utils::capture.output(print(fr)), collapse = "\n")
  expect_match(out, "polarity unresolved: Sim system one 1 of 3 testable gene\\(s\\) \\(G3\\)")
  expect_match(out, "weighted 0 in the system score")
  ## the audit itself
  a <- dmsa:::.rp_polarity_audit(fr)
  expect_true(isTRUE(attr(a, "has_polarity")))
  expect_equal(a$n_genes_testable, c(3L, 3L))
  expect_equal(a$n_polarity_signed, c(2L, 3L))
  expect_equal(a$n_polarity_unresolved, c(1L, 0L))
  expect_identical(a$genes_unresolved, c("G3", ""))
  ## the report: summary sentence + table
  r <- suppressWarnings(suppressMessages(dmsa_report(fr)))
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  line <- grep("System-level composition", sm, value = TRUE)
  expect_length(line, 1L)
  expect_match(line, "Sim system one - tone built from 2 of 3 testable gene\\(s\\), 1 with UNRESOLVED polarity \\(G3\\)")
  expect_match(line, "Sim system two - tone built from 3 of 3 testable gene\\(s\\)")
  expect_match(line, "NOT assumed \\+1")
  tb <- utils::read.csv(file.path(od, "tables", "polarity_audit.csv"),
                        stringsAsFactors = FALSE)
  expect_equal(tb$n_polarity_unresolved, c(1L, 0L))
  expect_true(any(grepl("tables/polarity_audit", r$tables)))
  ## and the system level really did run on 2 genes: G3's probes unusable
  sf <- dmsa:::.rp_system_frame(fr)
  al <- sf$alignment
  expect_false(any(al$usable[al$gene == "G3"]))
  expect_true(all(al$usable[al$gene %in% c("G1", "G2")]))
})

test_that("spec 40: a fully resolved table reports no unresolved genes and changes nothing", {
  s <- .pu_sim()
  pol <- data.frame(system_id = s$cascade$system_id, gene = s$cascade$gene,
                    w_g = c(1, -1, 1, 1, 1, -1))
  st <- dmsa_sets(s$cascade, polarity = pol)
  od <- tempfile("dmsa_pol40b")
  fr <- suppressMessages(dmsa_frame(
    s$data, map = s$map, sets = st, outcome = "out1",
    covariates = c("cov1", "cov2"), blocks = "cID",
    B = 99, seed = 1, outdir = od, plot_type = "png",
    progress = FALSE, beep = FALSE))
  out <- paste(utils::capture.output(print(fr)), collapse = "\n")
  expect_false(grepl("polarity unresolved", out))
  a <- dmsa:::.rp_polarity_audit(fr)
  expect_equal(a$n_polarity_unresolved, c(0L, 0L))
  r <- suppressWarnings(suppressMessages(dmsa_report(fr)))
  sm <- readLines(file.path(od, "summary.md"), warn = FALSE)
  expect_true(any(grepl("tone built from 3 of 3 testable gene\\(s\\); Sim system two - tone built from 3 of 3", sm)))
})
