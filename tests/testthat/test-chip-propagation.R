## Spec 44: a declared chip reaches every level, not only the finer ones.
## `.report_gene_level()` had always forwarded ri_group, so gene, module and
## probe were fitted with the chip random intercept while the SYSTEM level -
## the headline - was fitted without it, from the same frame and the same
## declaration. This pins the argument in place at every call site.

.chip_src <- function() {
  for (p in c(file.path("..", "..", "R"), "R", file.path("..", "..", "..", "R")))
    if (dir.exists(p) && length(Sys.glob(file.path(p, "*.R")))) return(p)
  NULL
}

test_that("every dmsa_triangulate CALL forwards ri_group", {
  dir <- .chip_src()
  skip_if(is.null(dir), "package sources not visible from the test directory")
  txt <- paste(unlist(lapply(Sys.glob(file.path(dir, "*.R")), readLines,
                             warn = FALSE)), collapse = "\n")
  st <- gregexpr("dmsa_triangulate\\(", txt)[[1]]
  st <- st[st > 0]
  expect_gte(length(st), 2L)
  checked <- 0L
  for (a in st) {
    rest <- substring(txt, a)
    ch <- strsplit(rest, "")[[1]]
    depth <- 0L; end <- NA_integer_
    for (i in seq_along(ch)) {
      if (ch[i] == "(") depth <- depth + 1L
      if (ch[i] == ")") { depth <- depth - 1L; if (depth == 0L) { end <- i; break } }
    }
    if (is.na(end)) next
    one <- substring(rest, 1, end)
    ## the definition itself, and bare prose mentions, are not call sites
    if (grepl("^dmsa_triangulate\\(M, data", one)) next
    if (grepl("^dmsa_triangulate\\(\\s*\\)$", one)) next
    if (!grepl(",", one)) next
    if (grepl("#'", one, fixed = TRUE)) next      # a roxygen example, not code
    checked <- checked + 1L
    expect_match(one, "ri_group", info = gsub("\\s+", " ", substr(one, 1, 90)))
  }
  expect_gte(checked, 2L)   # at least the system-level and gene-level calls
})

test_that("a declared random chip is retained on the frame", {
  mp <- utils::read.csv(system.file("extdata", "coverage_v4_full.csv",
                                    package = "dmsa"))
  mp <- mp[is.finite(mp$best_direction), ]
  cols <- utils::head(mp$column[mp$gene %in% c("FKBP5", "NR3C1")], 6)
  skip_if(length(cols) < 4)
  set.seed(1); n <- 60
  d <- data.frame(y = stats::rnorm(n), cov1 = stats::rnorm(n),
                  slide = factor(rep(1:4, length.out = n)),
                  cID = rep(seq_len(n / 2), each = 2))
  for (cl in cols) d[[cl]] <- stats::plogis(stats::rnorm(n))
  f <- dmsa_frame(d, methylation = cols, direction_source = "bundled",
                  outcome = "y", covariates = "cov1",
                  blocks = "cID", chip = "slide",
                  chip_effect = "random", B = 19, plots = FALSE,
                  tables = FALSE, summary = FALSE, progress = FALSE,
                  beep = FALSE, outdir = tempfile("dmsa_chip"))
  expect_true(nzchar(f$chip_random))
  expect_true(f$chip_random %in% names(f$data))
  expect_gt(nlevels(factor(f$data[[f$chip_random]])), 1L)
})
