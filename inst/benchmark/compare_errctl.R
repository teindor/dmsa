# ============================================================================
# HEAD-TO-HEAD: which error-control procedure should DMSA ship?
#   Rscript bench/compare_errctl.R <reps>
#
# Same data, same permutations, five procedures. Scored on the two things that
# decide it: can it NAME the gene carrying the signal, and does it name genes
# that carry none. Run at three concentrations because the whole argument is
# about a signal sitting in one gene rather than spread over twelve.
# ============================================================================
source("/home/claude/bench/generator.R")
source("/home/claude/bench/errctl.R")
suppressPackageStartupMessages(library(data.table))
A <- commandArgs(trailingOnly = TRUE)
REPS <- if (length(A) >= 1) as.integer(A[1]) else 200L
B <- 499L; ACC <- .85; K <- 60L; PER <- 5L; G <- K %/% PER
GID <- rep(seq_len(G), each = PER)
GRID <- CJ(conc = c(1L, 4L, 12L), h = c(0, .05))
OUT <- "/home/claude/outputs/bench_errctl.csv"; first <- TRUE

pool_z <- function(b, se, m) {
  w <- 1 / se^2; den <- sum(m^2 * w, na.rm = TRUE)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(m * b * w, na.rm = TRUE) / sqrt(den)
}
fit <- function(Y, x) {
  X <- cbind(1, x); XtXi <- solve(crossprod(X))
  bh <- XtXi %*% (t(X) %*% Y); r <- Y - X %*% bh
  se <- sqrt(colSums(r^2) / (nrow(Y) - 2) * XtXi[2, 2])
  list(b = bh[2, ], se = se)
}

for (i in seq_len(nrow(GRID))) {
  cfg <- GRID[i]
  gen <- make_gen(K = K, h = cfg$h, p_plus = .30, shape = "linear",
                  conc = cfg$conc, nbg = 5L)
  acc <- matrix(NA_real_, REPS, 12,
    dimnames = list(NULL, c("sys_gate", "sys_artp",
      "sens_single", "fwer_single", "sens_step", "fwer_step",
      "sens_bh", "fdr_bh", "tdp_true", "tdp_top1", "tdp_all", "tdp_wrong")))
  for (r in seq_len(REPS)) {
    g <- gen(); Y <- g$M; x <- g$x
    truth <- unique(GID[g$signal_probe])
    dc <- g$d * ifelse(stats::runif(K) < ACC, 1L, -1L)
    m <- dc * (2 * ACC - 1)
    L <- fit(Y, x)
    gz <- vapply(seq_len(G), function(j) { k <- which(GID == j)
      pool_z(L$b[k], L$se[k], m[k]) }, numeric(1))
    sysz <- pool_z(L$b, L$se, m)
    ## permutation null: gene-level z, system z, and probe-level t
    NG <- matrix(NA_real_, B, G); NS <- numeric(B)
    NT <- matrix(NA_real_, B, K)
    for (bb in seq_len(B)) {
      Lp <- fit(Y, x[sample.int(length(x))])
      NG[bb, ] <- vapply(seq_len(G), function(j) { k <- which(GID == j)
        pool_z(Lp$b[k], Lp$se[k], m[k]) }, numeric(1))
      NS[bb] <- pool_z(Lp$b, Lp$se, m)
      NT[bb, ] <- m * Lp$b / Lp$se          # aligned probe t, a sum statistic
    }
    ## ---- system level: current gate vs adaptive truncation ----------------
    p_dense  <- (1 + sum(abs(NS) >= abs(sysz))) / (B + 1)
    mx <- apply(abs(NG), 1, max, na.rm = TRUE)
    p_sparse <- (1 + sum(mx >= max(abs(gz), na.rm = TRUE))) / (B + 1)
    acc[r, "sys_gate"] <- as.numeric(gate_both(p_dense, p_sparse) < .05)
    acc[r, "sys_artp"] <- as.numeric(artp(gz, NG)$p < .05)
    ## ---- gene level: three ways to name a gene ---------------------------
    ps <- wy_single(gz, NG); pd <- wy_stepdown(gz, NG); qb <- bh_marginal(gz, NG)
    nm <- function(p) which(is.finite(p) & p < .05)
    for (nmv in list(c("single", "sens_single", "fwer_single"),
                     c("step",   "sens_step",   "fwer_step"))) {
      sel <- nm(if (nmv[1] == "single") ps else pd)
      acc[r, nmv[2]] <- as.numeric(all(truth %in% sel))
      acc[r, nmv[3]] <- as.numeric(any(!sel %in% truth))
    }
    selb <- nm(qb)
    acc[r, "sens_bh"] <- as.numeric(all(truth %in% selb))
    acc[r, "fdr_bh"]  <- if (!length(selb)) 0 else mean(!selb %in% truth)
    ## ---- TDP bounds on probe sets ----------------------------------------
    tobs <- m * L$b / L$se
    acc[r, "tdp_true"]  <- tdp_bound(tobs, NT, g$signal_probe)
    top1 <- which(GID == which.max(abs(gz)))
    acc[r, "tdp_top1"]  <- tdp_bound(tobs, NT, top1)
    acc[r, "tdp_all"]   <- tdp_bound(tobs, NT, seq_len(K))
    wrong <- which(!GID %in% truth)
    acc[r, "tdp_wrong"] <- tdp_bound(tobs, NT, wrong)
  }
  row <- as.data.table(c(list(conc = cfg$conc, h = cfg$h, reps = REPS, B = B),
                         as.list(colMeans(acc, na.rm = TRUE))))
  fwrite(row, OUT, append = !first); first <- FALSE
  cat(sprintf("conc=%2d h=%.2f | sys gate %.2f artp %.2f | NAME the gene: maxT-single %.2f (FWER %.3f)  maxT-step %.2f (FWER %.3f)  BH %.2f (FDR %.3f) | TDP true %.2f top1 %.2f all %.2f WRONG %.3f\n",
    cfg$conc, cfg$h, row$sys_gate, row$sys_artp, row$sens_single, row$fwer_single,
    row$sens_step, row$fwer_step, row$sens_bh, row$fdr_bh,
    row$tdp_true, row$tdp_top1, row$tdp_all, row$tdp_wrong))
}
cat("DONE-ERRCTL\n")
