# The gene model decides what a reader believes about WHERE a probe sits, and
# DMSA's whole premise is that where a probe sits determines which way its
# methylation points. These tests protect the two things that would be wrong
# silently: the 5'/3' UTR assignment on the minus strand, and features from a
# neighbouring gene surviving a region query.

## A real Ensembl GFF3 response for AVP (chr20, minus strand, 3 exons), kept
## verbatim so the parser is tested against the format it will actually meet.
avp_gff <- paste(
  "##gff-version 3",
  "##sequence-region   20 3082556 3084724",
  paste("20", "ensembl_havana", "gene", "3082556", "3084724", ".", "-", ".",
        "gene_id=ENSG00000101200;Name=AVP;canonical_transcript=ENST00000380293.4",
        sep = "\t"),
  paste("20", "ensembl_havana", "exon", "3084555", "3084724", ".", "-", ".",
        "Parent=transcript:ENST00000380293;rank=1;exon_id=ENSE00001484506",
        sep = "\t"),
  paste("20", "ensembl_havana", "exon", "3082977", "3083178", ".", "-", ".",
        "Parent=transcript:ENST00000380293;rank=2", sep = "\t"),
  paste("20", "ensembl_havana", "exon", "3082556", "3082802", ".", "-", ".",
        "Parent=transcript:ENST00000380293;rank=3", sep = "\t"),
  paste("20", "ensembl_havana", "CDS", "3084555", "3084664", ".", "-", "0",
        "Parent=transcript:ENST00000380293", sep = "\t"),
  paste("20", "ensembl_havana", "CDS", "3082700", "3082802", ".", "-", "0",
        "Parent=transcript:ENST00000380293", sep = "\t"),
  sep = "\n")

test_that("Ensembl GFF3 parses into the canonical gene-model shape", {
  gm <- .gm_from_gff(avp_gff, "AVP", genome = "hg38", text = TRUE)
  expect_true(all(.GM_COLS %in% names(gm)))
  expect_equal(unique(gm$gene), "AVP")
  expect_equal(unique(gm$chr), "chr20")
  expect_equal(unique(gm$strand), "-")
  expect_equal(unique(gm$transcript), "ENST00000380293")
  expect_true(all(gm$canonical))
  expect_equal(unique(gm$gene_start), 3082556)
  expect_equal(unique(gm$gene_end), 3084724)
  ## exon rank survives the CDS split, so a reader can name the exon a probe is in
  expect_equal(sort(unique(gm$exon_rank)), c(1, 2, 3))
})

test_that("on the minus strand the 5' UTR is the HIGH coordinate", {
  ## The trap. Exon 1 of AVP is 3084555-3084724 with coding stopping at
  ## 3084664, so 3084665-3084724 is the 5' UTR - the piece with the LARGER
  ## coordinate, because the gene reads right to left. Getting this backwards
  ## puts the promoter at the wrong end of every minus-strand gene, which is
  ## half of them.
  gm <- .gm_from_gff(avp_gff, "AVP", genome = "hg38", text = TRUE)
  u5 <- gm[gm$feature == "utr5", ]
  u3 <- gm[gm$feature == "utr3", ]
  expect_equal(nrow(u5), 1L); expect_equal(nrow(u3), 1L)
  expect_equal(u5$start, 3084665); expect_equal(u5$end, 3084724)
  expect_equal(u3$start, 3082556); expect_equal(u3$end, 3082699)
  expect_gt(u5$start, u3$start)
  ## the coding span runs from the first CDS base to the last, and no UTR may
  ## fall inside it - that is what makes the envelope test valid per exon
  cds <- gm[gm$feature == "cds", ]
  expect_equal(min(cds$start), 3082700)
  expect_equal(max(cds$end), 3084664)
  expect_true(all(u5$start > max(cds$end)))
  expect_true(all(u3$end < min(cds$start)))
})

test_that("the same exons on the PLUS strand give the mirror answer", {
  plus <- gsub("\t-\t", "\t+\t", avp_gff, fixed = TRUE)
  gm <- .gm_from_gff(plus, "FAKE", genome = "hg38", text = TRUE)
  expect_equal(gm$start[gm$feature == "utr5"], 3082556)
  expect_equal(gm$start[gm$feature == "utr3"], 3084665)
})

test_that("an exon with no CDS is reported as unresolved, not as coding", {
  ## Silence here would draw a full-height box and assert coding status the
  ## annotation never gave.
  no_cds <- paste(grep("CDS", strsplit(avp_gff, "\n")[[1]], invert = TRUE,
                       value = TRUE), collapse = "\n")
  gm <- .gm_from_gff(no_cds, "AVP", genome = "hg38", text = TRUE)
  expect_setequal(unique(gm$feature), "exon")
  expect_equal(nrow(gm), 3L)
  chk <- dmsa_gene_model_check(gm, verbose = FALSE)
  expect_true(chk$ok)
})

test_that("GTF attribute syntax is read as well as GFF3", {
  gtf <- paste(
    paste("20", "havana", "exon", "3084555", "3084724", ".", "-", ".",
          'gene_id "ENSG00000101200"; transcript_id "ENST00000380293"; exon_number "1";',
          sep = "\t"),
    paste("20", "havana", "exon", "3082556", "3082802", ".", "-", ".",
          'gene_id "ENSG00000101200"; transcript_id "ENST00000380293"; exon_number "3";',
          sep = "\t"), sep = "\n")
  gm <- .gm_from_gff(gtf, "AVP", genome = "hg38", text = TRUE)
  expect_equal(unique(gm$transcript), "ENST00000380293")
  expect_equal(nrow(gm), 2L)
})

test_that("the validator catches a model built from the wrong gene", {
  ## A region query returns every gene overlapping the span. If the neighbour's
  ## exons are not filtered out they are drawn as this gene's - the failure
  ## mode this check exists for.
  gm <- .gm_from_gff(avp_gff, "AVP", genome = "hg38", text = TRUE)
  bad <- gm
  bad$start[1] <- bad$gene_start[1] - 10000
  bad$end[1] <- bad$gene_start[1] - 9000
  chk <- dmsa_gene_model_check(bad, verbose = FALSE)
  expect_false(chk$ok)
  expect_false(chk$checks[["every feature inside the gene span"]]$ok)
  ## and an unknown feature word cannot slip in
  bad2 <- gm; bad2$feature[1] <- "promoter"
  expect_false(dmsa_gene_model_check(bad2, verbose = FALSE)$ok)
})

test_that("a model that cannot be obtained is empty, never invented", {
  gm <- dmsa_gene_model("NOT_A_REAL_GENE_XYZ", source = "ensembl",
                        host = NULL, cache = NULL, quiet = TRUE)
  expect_s3_class(gm, "dmsa_gene_model")
  expect_equal(nrow(gm), 0L)
  expect_output(print(gm), "empty")
})

test_that("source = 'table' round-trips and rejects a wrong shape", {
  gm <- .gm_from_gff(avp_gff, "AVP", genome = "hg38", text = TRUE)
  again <- dmsa_gene_model("AVP", source = "table", table = gm)
  expect_equal(nrow(again), nrow(gm))
  expect_s3_class(again, "dmsa_gene_model")
  expect_error(dmsa_gene_model("AVP", source = "table",
                               table = gm[, c("gene", "start")]),
               "missing column")
})

test_that("the Gviz engine refuses clearly rather than half-drawing", {
  skip_if(requireNamespace("Gviz", quietly = TRUE), "Gviz is installed")
  pr <- data.frame(probe = "cg1", pos = 3084600, b = .01,
                   stringsAsFactors = FALSE)
  expect_error(dmsa_plot_locus_gviz(pr, gene = "AVP"), "BiocManager::install")
  ## and it names the dependency-free alternative in the same breath
  expect_error(dmsa_plot_locus_gviz(pr, gene = "AVP"), "dmsa_plot_locus")
})

## ---- Ensembl lookup is JSON, parsed without a JSON dependency -------------
## lookup/symbol serves json/xml/jsonp ONLY; the old request for text/x-gff3
## was an unconditional HTTP 400 from Ensembl, on every machine (found
## 2026-08-29 when gene_models = "auto" first exercised the path live).

test_that(".gm_json_field reads the flat fields of a real lookup response", {
  loc <- paste0('{"source":"ensembl_havana","object_type":"Gene",',
                '"logic_name":"ensembl_havana_gene_homo_sapiens",',
                '"species":"homo_sapiens",',
                '"description":"prolactin receptor [Source:HGNC Symbol;Acc:HGNC:9446]",',
                '"display_name":"PRLR","assembly_name":"GRCh38",',
                '"biotype":"protein_coding","end":35230612,',
                '"seq_region_name":"5","db_type":"core","strand":-1,',
                '"id":"ENSG00000113494","start":35048829,"version":21,',
                '"canonical_transcript":"ENST00000618457.4"}')
  expect_identical(.gm_json_field(loc, "seq_region_name"), "5")
  expect_identical(.gm_json_field(loc, "start"), "35048829")
  expect_identical(.gm_json_field(loc, "end"), "35230612")
  expect_identical(.gm_json_field(loc, "id"), "ENSG00000113494")
  expect_identical(.gm_json_field(loc, "canonical_transcript"),
                   "ENST00000618457.4")
  expect_true(is.na(.gm_json_field(loc, "not_a_key")))
})
