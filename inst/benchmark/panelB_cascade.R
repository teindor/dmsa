# ============================================================================
# PANEL B - where the signal sits, and whether a method can say where.
#   Rscript bench/panelB_cascade.R <reps>
#
# Detection is only half the question. Once a set responds, the analyst wants to
# know WHICH genes and WHICH probes. Panel B holds the total signal constant and
# moves it from spread-across-all-genes to concentrated-in-one, then scores:
#
#   detection    did the method reject the set at all           (all engines)
#   gene layer   sensitivity and FDR over the 12 genes          (cascade vs BH)
#   probe layer  sensitivity and FDR over the 60 probes         (cascade vs BH)
#
# The competitors have no localization layer at all: the honest comparator is
# what an analyst actually does after a significant set - Benjamini-Hochberg on
# the probes, either genome-wide (strict, standard) or within the set (lenient,
# common but rarely defended).
# ============================================================================
source("/home/claude/bench/engines.R")
source("/home/claude/bench/generator.R")
suppressPackageStartupMessages(library(data.table))

A <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(A) >= 1) as.integer(A[1]) else 250L
ACC <- .85; NROT <- 200L; Q <- .05
K <- 60L; PER_GENE <- 5L; G <- K %/% PER_GENE          # 12 genes
PER_MODULE <- 3L                                        # 4 modules of 3 genes

GENE <- paste0("g", sprintf("%02d", rep(seq_len(G), each = PER_GENE)))
MOD  <- paste0("m", sprintf("%02d", (rep(seq_len(G), each = PER_GENE) - 1L) %/%
                              PER_MODULE + 1L))
TREE <- data.frame(system = "S", module = MOD, gene = GENE,
                   stringsAsFactors = FALSE)

CONC <- c(12L, 4L, 1L)                 # dense -> concentrated, matched signal
H    <- c(0, .03, .05)
GRID <- CJ(conc = CONC, h = H)
OUT  <- "/home/claude/outputs/bench_panelB.csv"
first <- !file.exists(OUT)

cat("calibrating flat null ...\n")
CRIT <- calibrate(make_gen0(K = K), NC = 1000, acc = ACC, nrot = NROT)

## the cascade needs its own permutation null; it depends only on the tree and
## the correlation structure, not on where the signal is, so one is enough.
cat("calibrating cascade null (this is the slow part) ...\n")
g0 <- make_gen0(K = K)()
f0 <- allfit(cbind(g0$M, g0$B), g0$x)
dc0 <- g0$d
AL0 <- dmsa_align(
  data.frame(cpg_id = paste0("p", seq_len(K)), d = dc0,
             p_plus = ifelse(dc0 > 0, ACC, 1 - ACC)),
  genes = GENE, level = "gene")
NULLS <- dmsa_cascade_null(f0$b[1:K], f0$se[1:K], AL0, TREE,
                           exposure = g0$x, M = g0$M, B = 299L)
cat("  cascade null ready\n")

## ---- helpers ---------------------------------------------------------------
gene_of_probe <- GENE
sens_fdr <- function(called, truth) {
  s <- if (!length(truth)) NA_real_ else mean(truth %in% called)
  f <- if (!length(called)) 0 else mean(!called %in% truth)
  c(sens = s, fdr = f)
}

for (i in seq_len(nrow(GRID))) {
  cfg <- GRID[i]
  gen <- make_gen(K = K, h = cfg$h, p_plus = .30, shape = "linear",
                  conc = cfg$conc)
  HIT <- matrix(NA_real_, REPS, length(ENGINES), dimnames = list(NULL, ENGINES))
  LOC <- matrix(NA_real_, REPS, 12, dimnames = list(NULL, c(
    "casc_sys", "casc_gene_sens", "casc_gene_fdr", "casc_probe_sens",
    "casc_probe_fdr", "bhgw_gene_sens", "bhgw_gene_fdr", "bhgw_probe_sens",
    "bhset_gene_sens", "bhset_gene_fdr", "bhset_probe_sens",
    "bhset_probe_fdr")))

  for (r in seq_len(REPS)) {
    g <- gen()
    e <- run_engines(g$M, g$B, g$x, g$d, g$gid, g$gidB, ACC, CRIT, nrot = NROT)
    HIT[r, ] <- e$hit
    true_genes  <- unique(GENE[g$signal_probe])
    true_probes <- g$signal_probe

    ## ---- cascade -----------------------------------------------------------
    dc <- sign(e$m_expected)
    al <- dmsa_align(
      data.frame(cpg_id = paste0("p", seq_len(K)), d = dc,
                 p_plus = ifelse(dc > 0, ACC, 1 - ACC)),
      genes = GENE, level = "gene")
    cs <- tryCatch(dmsa_cascade(e$b, e$se, al, TREE, NULLS, q = Q),
                   error = function(err) NULL)
    if (!is.null(cs)) {
      st <- cs$tables$system
      LOC[r, "casc_sys"] <- as.numeric(!is.null(st) && any(st$selected))
      gt <- cs$tables$gene
      cg <- if (is.null(gt)) character(0) else
        vapply(strsplit(gt$unit[gt$selected], "|", fixed = TRUE),
               function(z) z[length(z)], character(1))
      sf <- sens_fdr(cg, true_genes)
      LOC[r, "casc_gene_sens"] <- sf[["sens"]]
      LOC[r, "casc_gene_fdr"]  <- sf[["fdr"]]
      sfp <- sens_fdr(as.integer(cs$selected_probes), true_probes)
      LOC[r, "casc_probe_sens"] <- sfp[["sens"]]
      LOC[r, "casc_probe_fdr"]  <- sfp[["fdr"]]
    }

    ## ---- what an analyst does instead --------------------------------------
    ## genome-wide BH over set + background, gene called if any probe survives
    pgw <- stats::p.adjust(e$p, "BH")[seq_len(K)]
    gw <- which(pgw < Q)
    LOC[r, "bhgw_probe_sens"] <- sens_fdr(gw, true_probes)[["sens"]]
    sfg <- sens_fdr(unique(GENE[gw]), true_genes)
    LOC[r, "bhgw_gene_sens"] <- sfg[["sens"]]
    LOC[r, "bhgw_gene_fdr"]  <- sfg[["fdr"]]
    ## BH within the set only - lenient, common, rarely defended
    pin <- stats::p.adjust(e$p[1:K], "BH")
    ws <- which(pin < Q)
    sfs <- sens_fdr(unique(GENE[ws]), true_genes)
    LOC[r, "bhset_gene_sens"] <- sfs[["sens"]]
    LOC[r, "bhset_gene_fdr"]  <- sfs[["fdr"]]
    sfw <- sens_fdr(ws, true_probes)
    LOC[r, "bhset_probe_sens"] <- sfw[["sens"]]
    LOC[r, "bhset_probe_fdr"]  <- sfw[["fdr"]]
  }

  row <- data.table(panel = "B", shape = "linear", K = K, h = cfg$h,
                    conc = cfg$conc, p_plus = .30, reps = REPS)
  for (e in ENGINES) row[[paste0("pow_", e)]] <- mean(HIT[, e], na.rm = TRUE)
  for (cn in colnames(LOC)) row[[cn]] <- mean(LOC[, cn], na.rm = TRUE)
  fwrite(row, OUT, append = !first); first <- FALSE
  cat(sprintf("conc=%2d h=%.2f | DMSA %.2f gt %.2f cam %.2f ewas %.2f | casc sys %.2f gene %.2f/%.3f probe %.2f/%.3f | BHgw gene %.2f/%.3f | BHset gene %.2f/%.3f\n",
    cfg$conc, cfg$h, row$pow_dmsa_expected, row$pow_globaltest, row$pow_camera,
    row$pow_ewas_bh, row$casc_sys, row$casc_gene_sens, row$casc_gene_fdr,
    row$casc_probe_sens, row$casc_probe_fdr, row$bhgw_gene_sens,
    row$bhgw_gene_fdr, row$bhset_gene_sens, row$bhset_gene_fdr))
}
cat("DONE-PANEL-B\n")
