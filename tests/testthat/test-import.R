## dmsa_import(): the pipeline adapter. Orientation, value scale and sample
## alignment must each be either unambiguous or an error - never a guess.

.imp_beta <- function(np = 8, ns = 10, seed = 7) {
  set.seed(seed)
  matrix(runif(np * ns), nrow = np,
         dimnames = list(sprintf("cg%08d", seq_len(np)),
                         paste0("S", seq_len(ns))))
}

test_that("probes-in-rows matrix (the pipeline shape) is transposed", {
  B <- .imp_beta()
  imp <- dmsa_import(B)
  expect_s3_class(imp, "dmsa_import")
  expect_equal(dim(imp$methylation), c(10, 8))
  expect_true(all(grepl("^cg", colnames(imp$methylation))))
  expect_identical(rownames(imp$methylation), paste0("S", 1:10))
  expect_identical(imp$values, "betas")
})

test_that("samples-in-rows matrix passes through unchanged", {
  B <- t(.imp_beta())
  imp <- dmsa_import(B)
  expect_equal(dim(imp$methylation), c(10, 8))
})

test_that("ambiguous orientation is an error, and a declaration resolves it", {
  B <- .imp_beta()
  dimnames(B) <- NULL
  expect_error(dmsa_import(B), "cannot tell which axis")
  imp <- dmsa_import(B, orientation = "probes_rows",
                     pheno = data.frame(x = 1:10), align_by = ".position")
  expect_equal(dim(imp$methylation), c(10, 8))
})

test_that("probe-looking names on BOTH axes is an error", {
  B <- .imp_beta(6, 6)
  colnames(B) <- sprintf("cg%08d", 101:106)
  expect_error(dmsa_import(B), "BOTH")
})

test_that("M-values are detected and converted to betas", {
  B <- .imp_beta()
  M <- log2(B / (1 - B))
  imp <- dmsa_import(M)
  expect_match(imp$values, "converted")
  expect_true(all(imp$methylation >= 0 & imp$methylation <= 1))
  expect_equal(unname(imp$methylation), unname(t(B)), tolerance = 1e-10)
})

test_that("pheno rows are matched and reordered to the sample order", {
  B <- .imp_beta()
  ph <- data.frame(Sample_Name = paste0("S", 10:1), age = 10:1)
  imp <- dmsa_import(B, pheno = ph)
  expect_identical(imp$data$Sample_Name, paste0("S", 1:10))
  expect_identical(imp$data$age, 1:10)
  expect_match(imp$align, "Sample_Name")
})

test_that("minfi-style Basename paths and Sentrix pairs are matched", {
  B <- .imp_beta()
  ph1 <- data.frame(Basename = file.path("/idats/run1", paste0("S", 1:10)),
                    grp = rep(1:2, 5))
  expect_match(dmsa_import(B, pheno = ph1)$align, "basename")
  B2 <- .imp_beta()
  colnames(B2) <- paste0("205060", 1:10, "_R0", rep(1:5, 2), "C01")
  ph2 <- data.frame(Sentrix_ID = paste0("205060", 1:10),
                    Sentrix_Position = paste0("R0", rep(1:5, 2), "C01"))
  expect_match(dmsa_import(B2, pheno = ph2)$align, "Sentrix")
})

test_that("an unmatchable pheno is an error naming what was tried", {
  B <- .imp_beta()
  ph <- data.frame(who = paste0("person", 1:10))
  expect_error(dmsa_import(B, pheno = ph), "align_by")
  expect_error(dmsa_import(B, pheno = ph), "who")
})

test_that("positional alignment must be declared, and checks row counts", {
  B <- .imp_beta()
  ph <- data.frame(y = rnorm(10))
  imp <- dmsa_import(B, pheno = ph, align_by = ".position")
  expect_match(imp$align, "positional")
  expect_error(dmsa_import(B, pheno = data.frame(y = rnorm(9)),
                           align_by = ".position"), "same number of rows")
})

test_that("no pheno at all still yields a working object", {
  imp <- dmsa_import(.imp_beta())
  expect_identical(imp$data$.sample, paste0("S", 1:10))
})

test_that("data.frame exports: probe id column found, annotation dropped", {
  B <- .imp_beta()
  df <- data.frame(IlmnID = rownames(B), CHR = "5", as.data.frame(B),
                   check.names = FALSE)
  imp <- dmsa_import(df)
  expect_equal(dim(imp$methylation), c(10, 8))
  expect_true(any(grepl("non-numeric", imp$notes)))
})

test_that("csv and rds paths import", {
  B <- .imp_beta()
  df <- data.frame(Probe_ID = rownames(B), as.data.frame(B),
                   check.names = FALSE)
  f <- tempfile(fileext = ".csv"); utils::write.csv(df, f, row.names = FALSE)
  expect_equal(dim(dmsa_import(f)$methylation), c(10, 8))
  r <- tempfile(fileext = ".rds"); saveRDS(B, r)
  expect_equal(dim(dmsa_import(r)$methylation), c(10, 8))
  expect_error(dmsa_import(tempfile(fileext = ".xlsx")), "file not found")
})

test_that("SummarizedExperiment: Beta assay, colData as sheet", {
  B <- .imp_beta()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Beta = B),
    colData = data.frame(age = 1:10, row.names = colnames(B)))
  imp <- dmsa_import(se)
  expect_equal(dim(imp$methylation), c(10, 8))
  expect_identical(imp$data$age, 1:10)
  expect_match(imp$source, "Beta")
})

test_that("SummarizedExperiment: M assay converted, Meth/Unmeth computed", {
  B <- .imp_beta()
  M <- log2(B / (1 - B))
  se_m <- SummarizedExperiment::SummarizedExperiment(assays = list(M = M))
  imp_m <- dmsa_import(se_m)
  expect_equal(unname(imp_m$methylation), unname(t(B)), tolerance = 1e-10)
  Meth <- B * 1000; Unmeth <- (1 - B) * 1000
  se_mu <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Meth = Meth, Unmeth = Unmeth))
  imp_mu <- dmsa_import(se_mu)
  expect_equal(unname(imp_mu$methylation),
               unname(t(Meth / (Meth + Unmeth + 100))), tolerance = 1e-10)
})

test_that("raw channel containers are refused with a preprocessing recipe", {
  B <- .imp_beta()
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Green = B, Red = B))
  expect_error(dmsa_import(se), "preprocess")
  fake_rg <- structure(list(), class = "RGChannelSet")
  expect_error(dmsa_import(fake_rg), "RGChannelSet")
})

test_that("sesame SigDF lists and RnBeads objects get recipes, not guesses", {
  sig <- data.frame(Probe_ID = "cg00000001", MG = 1, MR = 1, UG = 1, UR = 1)
  expect_error(dmsa_import(list(sig, sig)), "openSesame")
  fake_rnb <- structure(list(), class = "RnBeadSet")
  expect_error(dmsa_import(fake_rnb), "RnBeads")
  expect_error(dmsa_import(list(1, 2)), "SigDF")
})

test_that("print method reports shape, source and alignment", {
  imp <- dmsa_import(.imp_beta())
  expect_output(print(imp), "10 samples x 8 probes")
  expect_output(print(imp), "dmsa_frame")
})

test_that("dmsa_frame() accepts a dmsa_import object end to end", {
  set.seed(42)
  map <- data.frame(gene = rep(paste0("G", 1:3), each = 3), system_id = 1L,
                    system = "Sim system")
  map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
  map$column <- paste0(map$probe, "_", map$gene)
  map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
  map$p_plus <- ifelse(map$best_direction > 0, 0.9, 0.1)
  n <- 60
  ph <- data.frame(out1 = rnorm(n), cov1 = rnorm(n),
                   cID = rep(1:30, each = 2),
                   Sample_Name = sprintf("S%02d", 1:n))
  sig <- 0.5 * outer(ph$out1, map$best_direction * (map$gene == "G1"))
  B <- t(plogis(matrix(rnorm(n * nrow(map)), n) + sig))   # probes x samples
  dimnames(B) <- list(map$probe, ph$Sample_Name)
  imp <- dmsa_import(B, pheno = ph)
  expect_equal(dim(imp$methylation), c(n, nrow(map)))
  fr <- dmsa_frame(imp, map = map, outcome = "out1", covariates = "cov1",
                   random_effects = "cID", B = 99, plots = FALSE,
                   tables = FALSE, summary = FALSE,
                   outdir = tempfile("dmsa_imp"))
  expect_s3_class(fr, "dmsa_frame")
})
