# ============================================================================
# PANEL A - the direction-ratio mechanism, and what it does to every method.
#   Rscript bench/panelA_ratio.R <reps> <chunk>     chunk in {1,2}
#
# Sweeps p_plus across the range the eQTM literature actually reports
# (promoter-rich .30, whole-blood .41, Alpha .26, body-rich .74) plus the
# 50/50 tissue where every unaligned method must be null by construction.
# ============================================================================
source("/home/claude/bench/engines.R")
source("/home/claude/bench/generator.R")
suppressPackageStartupMessages(library(data.table))

A <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(A) >= 1) as.integer(A[1]) else 300L
CHUNK <- if (length(A) >= 2) as.integer(A[2]) else 1L
ACC <- .85; NROT <- 200L

P_PLUS <- c(.20, .26, .30, .41, .50, .60, .74)
H      <- c(0, .02, .03, .04)
KGRID  <- 60L
GRID <- CJ(p_plus = P_PLUS, h = H, K = KGRID)
GRID <- GRID[seq(CHUNK, .N, by = 2)]

OUT <- sprintf("/home/claude/outputs/bench_panelA_%d.csv", CHUNK)
first <- !file.exists(OUT)

## the null is the same for every p_plus (no effect => no ratio), so calibrate
## once per K and reuse
cat("calibrating null critical values ...\n")
CRIT <- calibrate(make_gen0(K = KGRID), NC = 1000, acc = ACC, nrot = NROT)
print(round(CRIT, 3))

for (i in seq_len(nrow(GRID))) {
  cfg <- GRID[i]
  gen <- make_gen(K = cfg$K, h = cfg$h, p_plus = cfg$p_plus, shape = "linear")
  HIT <- matrix(NA_real_, REPS, length(ENGINES),
                dimnames = list(NULL, ENGINES))
  SGN <- HIT
  for (r in seq_len(REPS)) {
    g <- gen()
    e <- run_engines(g$M, g$B, g$x, g$d, g$gid, g$gidB, ACC, CRIT, nrot = NROT)
    HIT[r, ] <- e$hit; SGN[r, ] <- e$rep_sign
  }
  ## two truths, deliberately different (see engines.R)
  truth_meth <- sign(2 * cfg$p_plus - 1)     # what an unaligned method can see
  truth_bio  <- 1                            # what the exposure actually does
  cond_sign <- function(eng, truth) {
    k <- which(HIT[, eng] == 1)
    if (!length(k) || all(is.na(SGN[k, eng]))) return(NA_real_)
    mean(SGN[k, eng] == truth, na.rm = TRUE)
  }
  row <- data.table(panel = "A", shape = "linear", K = cfg$K, h = cfg$h,
                    p_plus = cfg$p_plus, reps = REPS)
  for (e in ENGINES) row[[paste0("pow_", e)]] <- mean(HIT[, e], na.rm = TRUE)
  for (e in DIRECTIONAL) {
    row[[paste0("bio_", e)]]  <- cond_sign(e, truth_bio)
    row[[paste0("meth_", e)]] <- cond_sign(e, truth_meth)
  }
  fwrite(row, OUT, append = !first); first <- FALSE
  cat(sprintf("p+=%.2f h=%.2f | DMSA %.2f nDMSA %.2f blind %.2f | cam %.2f fry %.2f roastD %.2f roastM %.2f gt %.2f gsea %.2f RRA %.2f ora %.2f ewas %.2f | bioDir DMSA %s cam %s\n",
    cfg$p_plus, cfg$h, row$pow_dmsa_expected, row$pow_ndmsa,
    row$pow_dmsa_signblind, row$pow_camera, row$pow_fry, row$pow_roast_dir,
    row$pow_roast_mixed, row$pow_globaltest, row$pow_gsea_mcsea,
    row$pow_methylRRA, row$pow_ora, row$pow_ewas_bh,
    ifelse(is.na(row$bio_dmsa_expected), "--", sprintf("%.2f", row$bio_dmsa_expected)),
    ifelse(is.na(row$bio_camera), "--", sprintf("%.2f", row$bio_camera))))
}
cat("DONE-PANEL-A chunk", CHUNK, "\n")
