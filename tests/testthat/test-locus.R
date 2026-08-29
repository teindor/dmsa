test_that("dmsa_plot_locus works with no coordinates at all", {
  set.seed(1)
  p <- data.frame(probe = paste0("cg", 1:6), b = rnorm(6, .02, .01),
                  se = rep(.01, 6), d = c(1, 1, -1, -1, 1, -1),
                  stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png")
  out <- dmsa_plot_locus(p, gene = "AVP", file = f)
  expect_true(file.exists(f))
  expect_gt(file.size(f), 1000)
  ## evenly spaced in the order supplied
  expect_identical(out$x, seq_len(6))
  expect_identical(out$probe, p$probe)
  unlink(f)
})

test_that("gene_region orders probes 5' to 3' when pos is absent", {
  p <- data.frame(probe = paste0("cg", 1:5), b = rep(.02, 5), se = rep(.01, 5),
                  d = c(1, -1, 1, -1, 1),
                  gene_region = c("Body", "TSS200", "3'UTR", "TSS1500",
                                  "1stExon"),
                  stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png")
  out <- dmsa_plot_locus(p, gene = "G", file = f)
  expect_identical(out$gene_region,
                   c("TSS1500", "TSS200", "1stExon", "Body", "3'UTR"))
  expect_identical(out$x, seq_len(5))
  unlink(f)
})

test_that("pos is used when present and forces genomic order", {
  p <- data.frame(probe = paste0("cg", 1:4), b = rep(.02, 4), se = rep(.01, 4),
                  d = c(1, -1, 1, -1), pos = c(300, 100, 400, 200),
                  chr = rep("20", 4), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png")
  out <- dmsa_plot_locus(p, gene = "AVP", file = f)
  expect_identical(out$pos, c(100, 200, 300, 400))
  expect_identical(out$x, out$pos)
  expect_true(file.exists(f))
  unlink(f)
})

test_that("partial coordinates never drop a probe from the figure", {
  p <- data.frame(probe = paste0("cg", 1:4), b = rep(.02, 4), se = rep(.01, 4),
                  d = c(1, -1, 1, -1), pos = c(300, NA, 400, 200),
                  stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png")
  expect_message(out <- dmsa_plot_locus(p, file = f), "spacing all of them")
  expect_equal(nrow(out), 4L)          # nothing dropped
  expect_identical(out$x, seq_len(4))  # even spacing, not the genomic axis
  ## forcing the genomic axis is allowed, and then the unplaceable one goes
  out2 <- dmsa_plot_locus(p, order_by = "pos", file = f)
  expect_equal(nrow(out2), 3L)
  unlink(f)
})

test_that("order_by errors rather than silently falling back", {
  p <- data.frame(probe = "cg1", b = .02, se = .01, d = 1,
                  stringsAsFactors = FALSE)
  expect_error(dmsa_plot_locus(p, order_by = "pos"), "no usable 'pos'")
  expect_error(dmsa_plot_locus(p, order_by = "region"), "gene_region")
})

test_that("invert reflects the right probes regardless of layout mode", {
  p <- data.frame(probe = c("a", "b"), b = c(.02, .02), se = c(.01, .01),
                  d = c(1, -1), stringsAsFactors = FALSE)
  f <- tempfile(fileext = ".png")
  o1 <- dmsa_plot_locus(p, invert = "+1", file = f)
  expect_equal(o1$shown, c(-.02, .02))
  o2 <- dmsa_plot_locus(p, invert = "none", file = f)
  expect_equal(o2$shown, c(.02, .02))
  unlink(f)
})

test_that("dmsa_probe_coords refuses to guess and reports the columns", {
  man <- tempfile(fileext = ".tsv")
  write.table(data.frame(weird = "cg1", other = 1),
              man, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(dmsa_probe_coords("cg1", file = man), "weird")
  unlink(man)
})

test_that("dmsa_probe_coords reads a plain manifest and keeps probe order", {
  man <- tempfile(fileext = ".tsv")
  write.table(data.frame(probeID = c("cg2", "cg1", "cg3"),
                         CpG_chrm = c("chr20", "chr20", "chr20"),
                         CpG_beg = c(200, 100, 300)),
              man, sep = "\t", row.names = FALSE, quote = FALSE)
  co <- dmsa_probe_coords(c("cg1", "cg3", "cgX"), file = man)
  expect_identical(co$probe, c("cg1", "cg3", "cgX"))
  ## E10 (PI-approved): CpG_beg is BED-convention 0-based and is now shifted
  ## to the 1-based scale everything else in the locus panel uses
  expect_equal(co$pos, c(101, 301, NA_real_))
  expect_identical(co$chr[1], "20")
  unlink(man)
})

## --- the gene model in strip B ---------------------------------------------

.mini_model <- function(strand = "-") {
  gff <- paste(
    paste("20", "hav", "gene", "3082556", "3084724", ".", strand, ".",
          "gene_id=ENSG1;canonical_transcript=ENST1.1", sep = "\t"),
    paste("20", "hav", "exon", "3084555", "3084724", ".", strand, ".",
          "Parent=transcript:ENST1;rank=1", sep = "\t"),
    paste("20", "hav", "exon", "3082977", "3083178", ".", strand, ".",
          "Parent=transcript:ENST1;rank=2", sep = "\t"),
    paste("20", "hav", "exon", "3082556", "3082802", ".", strand, ".",
          "Parent=transcript:ENST1;rank=3", sep = "\t"),
    paste("20", "hav", "CDS", "3084555", "3084664", ".", strand, "0",
          "Parent=transcript:ENST1", sep = "\t"),
    paste("20", "hav", "CDS", "3082700", "3082802", ".", strand, "0",
          "Parent=transcript:ENST1", sep = "\t"), sep = "\n")
  .gm_from_gff(gff, "AVP", genome = "hg38", text = TRUE)
}

.mini_probes <- function() {
  data.frame(probe = paste0("cg", 1:5), chr = "chr20",
             pos = c(3082600, 3082750, 3083050, 3084600, 3084700),
             b = c(-.04, .02, .03, -.05, -.06), se = rep(.015, 5),
             d = c(-1, 1, 1, -1, -1), p = c(.01, .3, .04, .002, .001),
             stringsAsFactors = FALSE)
}

test_that("a gene model is drawn and the panel says whose it is", {
  f <- file.path(tempdir(), "locus_model.png"); on.exit(unlink(f), add = TRUE)
  out <- dmsa_plot_locus(.mini_probes(), gene = "AVP",
                         gene_model = .mini_model(), file = f)
  expect_true(file.exists(f) && file.size(f) > 5000)
  expect_equal(nrow(out), 5L)
  ## the note names the source, the genome and the transcript actually drawn
  note <- .locus_model_note(.mini_model(), "canonical")
  expect_match(note, "hg38")
  expect_match(note, "ENST1")
  expect_match(note, "chevrons")
})

test_that("the window opens to the whole gene, not just the probes", {
  ## Cropping to the probes would cut off the exons that explain where those
  ## probes sit - the reason for drawing the model at all.
  f <- file.path(tempdir(), "locus_win.png"); on.exit(unlink(f), add = TRUE)
  gm <- .mini_model()
  pr <- .mini_probes()[1:2, ]          # both probes at the far 3' end
  out <- dmsa_plot_locus(pr, gene = "AVP", gene_model = gm, file = f)
  expect_true(file.exists(f))
  expect_equal(nrow(out), 2L)
})

test_that("a model for the wrong chromosome is refused, not drawn", {
  f <- file.path(tempdir(), "locus_wrongchr.png"); on.exit(unlink(f), add = TRUE)
  gm <- .mini_model(); gm$chr <- "chr7"
  expect_message(dmsa_plot_locus(.mini_probes(), gene = "AVP", gene_model = gm,
                                 file = f), "not drawing it")
  expect_true(file.exists(f))
})

test_that("a model without probe positions is refused rather than faked", {
  f <- file.path(tempdir(), "locus_nopos.png"); on.exit(unlink(f), add = TRUE)
  pr <- .mini_probes(); pr$pos <- NULL; pr$chr <- NULL
  expect_message(
    suppressMessages(dmsa_plot_locus(pr, gene = "AVP", gene_model = .mini_model(),
                                     file = f)) ,
    NA)
  expect_true(file.exists(f))
})

test_that("the evenly-spaced fallback announces itself", {
  ## It must never be a silent default: it looks like a map and is not one.
  f <- file.path(tempdir(), "locus_fallback.png"); on.exit(unlink(f), add = TRUE)
  pr <- .mini_probes(); pr$pos <- NULL; pr$chr <- NULL
  expect_message(dmsa_plot_locus(pr, gene = "AVP", file = f),
                 "cannot be drawn to scale")
})

test_that("both strands put the TSS at the right end", {
  ## The minus-strand TSS is the HIGH coordinate. This is the check that would
  ## have caught a mirrored gene track.
  m <- .mini_model("-"); p <- .mini_model("+")
  expect_equal(unique(m$strand), "-")
  expect_equal(unique(p$strand), "+")
  expect_gt(m$start[m$feature == "utr5"], m$start[m$feature == "utr3"])
  expect_lt(p$start[p$feature == "utr5"], p$start[p$feature == "utr3"])
})

test_that("all transcripts can be stacked", {
  f <- file.path(tempdir(), "locus_iso.png"); on.exit(unlink(f), add = TRUE)
  gm <- .mini_model()
  gm2 <- gm; gm2$transcript <- "ENST2"; gm2$canonical <- FALSE
  gm2 <- gm2[gm2$feature != "utr3", ]
  both <- rbind(gm, gm2)
  out <- dmsa_plot_locus(.mini_probes(), gene = "AVP", gene_model = both,
                         transcripts = "all", file = f)
  expect_true(file.exists(f) && file.size(f) > 5000)
  ## and "canonical" narrows it back to one
  expect_match(.locus_model_note(both, "canonical"), "canonical transcript ENST1")
  expect_match(.locus_model_note(both, "all"), "all 2 transcripts")
})

test_that("a genome-build mismatch is caught, not drawn silently", {
  ## The real one: cpgdirection ships hg19 EPIC coordinates, the cascade and
  ## Ensembl are hg38, and AVP moves 19,353 bp between them. Both layers are
  ## individually correct, so nothing errors - the CpGs simply land 19 kb from
  ## their own gene. This is the check that catches it.
  f <- file.path(tempdir(), "locus_build.png"); on.exit(unlink(f), add = TRUE)
  gm <- .mini_model()                                   # hg38, chr20 3082556-3084724
  pr <- .mini_probes()
  pr$pos <- pr$pos - 19353                              # the hg19 positions
  expect_message(dmsa_plot_locus(pr, gene = "AVP", gene_model = gm, file = f),
                 "GENOME BUILD MISMATCH")
  expect_message(dmsa_plot_locus(pr, gene = "AVP", gene_model = gm, file = f),
                 "hg38")
})

test_that("real promoter probes outside the gene are NOT flagged", {
  ## All 8 of AVP's probes sit upstream of the TSS, outside the gene span. A
  ## naive "is every probe inside the gene" rule would cry wolf on the very
  ## figure this feature exists to draw, so the rule is distance vs gene length.
  f <- file.path(tempdir(), "locus_prom.png"); on.exit(unlink(f), add = TRUE)
  gm <- .mini_model()
  pr <- data.frame(probe = paste0("cg", 1:8), chr = "chr20",
                   pos = c(3084696, 3084756, 3084776, 3084826,
                           3084841, 3084912, 3084935, 3084937),
                   b = rep(.02, 8), se = rep(.01, 8), d = rep(-1, 8),
                   stringsAsFactors = FALSE)
  expect_silent(dmsa_plot_locus(pr, gene = "AVP", gene_model = gm, file = f))
})
