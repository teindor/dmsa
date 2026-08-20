# ============================================================================
# PANEL C - what happens when the exposure does not act linearly.
#   Rscript bench/panelC_shape.R <reps>
#
# Every method in the roster except nDMSA tests a LINEAR association between the
# exposure and methylation. Two shapes break that assumption in ways the
# substantive literature has reason to expect:
#
#   quadratic  an optimum in mid-range - too little and too much both cost. Its
#              linear component is close to zero by construction, so a linear
#              test is not merely weaker, it is testing for something that is
#              not there.
#   threshold  nothing happens until a level is crossed, then a step. Allostatic
#              load: sustained for a long time, then something gives.
#
# nDMSA pays a real multiplicity price for looking in three places at once, and
# the linear row is where that price is charged. Reporting that cost honestly is
# the point of including the linear shape here.
# ============================================================================
source("/home/claude/bench/engines.R")
source("/home/claude/bench/generator.R")
suppressPackageStartupMessages(library(data.table))

A <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(A) >= 1) as.integer(A[1]) else 300L
ACC <- .85; NROT <- 200L; K <- 60L

SHAPES <- c("linear", "quadratic", "threshold")
H      <- c(0, .03, .05, .08)
GRID <- CJ(shape = SHAPES, h = H, sorted = FALSE)
OUT  <- "/home/claude/outputs/bench_panelC.csv"
first <- !file.exists(OUT)
## RESUME (container restarts are frequent; this grid is longer than the
## interval between them, so without this the run never reaches the end).
DONE <- if (file.exists(OUT)) {
  dn <- tryCatch(fread(OUT), error = function(e) NULL)
  if (is.null(dn) || !nrow(dn)) character(0) else paste(dn$shape, dn$h)
} else character(0)
if (length(DONE)) cat("resume: skipping", length(DONE), "completed panelC cells\n")

## One null per shape: under h = 0 the shape never enters the data, so the null
## is genuinely shared. Calibrating once keeps every row on the same threshold.
CRITF <- "/home/claude/outputs/bench_panelC_crit.rds"
if (file.exists(CRITF)) {
  CRIT <- readRDS(CRITF); cat("null calibration: loaded from cache\n")
} else {
  cat("calibrating null ...\n")
  CRIT <- calibrate(make_gen0(K = K), NC = 1000, acc = ACC, nrot = NROT)
  saveRDS(CRIT, CRITF)
}
print(round(CRIT, 3))

for (i in seq_len(nrow(GRID))) {
  if (paste(GRID$shape[i], GRID$h[i]) %in% DONE) next
  cfg <- GRID[i]
  gen <- make_gen(K = K, h = cfg$h, p_plus = .30, shape = cfg$shape)
  HIT <- matrix(NA_real_, REPS, length(ENGINES), dimnames = list(NULL, ENGINES))
  ARM <- matrix(NA_real_, REPS, 3,
                dimnames = list(NULL, c("linear", "quadratic", "threshold")))
  for (r in seq_len(REPS)) {
    g <- gen()
    e <- run_engines(g$M, g$B, g$x, g$d, g$gid, g$gidB, ACC, CRIT, nrot = NROT)
    HIT[r, ] <- e$hit
    ARM[r, ] <- e$nd[c("linear", "quadratic", "threshold")]
  }
  ## which nDMSA arm won, among replicates where nDMSA fired - this is the
  ## diagnostic that tells the analyst WHAT SHAPE the association took
  fired <- which(HIT[, "ndmsa"] == 1)
  won <- if (length(fired))
    table(factor(colnames(ARM)[max.col(ARM[fired, , drop = FALSE])],
                 levels = colnames(ARM))) / length(fired) else
    stats::setNames(rep(NA_real_, 3), colnames(ARM))

  row <- data.table(panel = "C", shape = cfg$shape, K = K, h = cfg$h,
                    p_plus = .30, reps = REPS)
  for (e in ENGINES) row[[paste0("pow_", e)]] <- mean(HIT[, e], na.rm = TRUE)
  row$arm_linear <- as.numeric(won[["linear"]])
  row$arm_quadratic <- as.numeric(won[["quadratic"]])
  row$arm_threshold <- as.numeric(won[["threshold"]])
  fwrite(row, OUT, append = !first); first <- FALSE
  cat(sprintf("%-10s h=%.2f | nDMSA %.2f DMSA %.2f blind %.2f | gt %.2f cam %.2f fry %.2f roastD %.2f gsea %.2f RRA %.2f ora %.2f ewas %.2f | arm lin/quad/thr %.2f/%.2f/%.2f\n",
    cfg$shape, cfg$h, row$pow_ndmsa, row$pow_dmsa_expected,
    row$pow_dmsa_signblind, row$pow_globaltest, row$pow_camera, row$pow_fry,
    row$pow_roast_dir, row$pow_gsea_mcsea, row$pow_methylRRA, row$pow_ora,
    row$pow_ewas_bh, row$arm_linear, row$arm_quadratic, row$arm_threshold))
}
cat("DONE-PANEL-C\n")
