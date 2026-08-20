# ============================================================================
# COMPETITOR RE-RUN - v6 (Mac). Source this file.
#
# WHY THIS EXISTS. The v5 run's DMSA columns are VALID and are NOT recomputed
# here. Its COMPETITOR columns are not, because of three defects in the v5
# harness - all mine, all in the competitor block, all confirmed by direct
# test. v6 fixes them and re-runs the competitors only, on the same grid and
# the same seeds, then merges with your existing v5 DMSA rows.
#
# (1) camera ran with limma's DEFAULT inter.gene.cor = 0.01 instead of the
#     ESTIMATED value. The manuscript says "camera with estimated inter-gene
#     correlation", and the locked benchmark does pass inter.gene.cor = NA
#     (bench/engines.R). v5 did not. With true within-gene correlation .6,
#     assuming .01 makes camera wildly anti-conservative. Measured null
#     family-wise error, 200 reps, n = 250:
#         inter.gene.cor = 0.01 -> .860 (G = 12), 1.000 (G = 42)
#         inter.gene.cor = NA   -> .060 (G = 12),  .045 (G = 42)
#     This is why v5 reported camera nulls of .65-1.00. Every camera power
#     number and every camera N80 in v5 is void.
#
# (2) ORA used a raw Fisher exact p. Fisher assumes independent probes; within
#     a gene they correlate at .6, so the counts are over-dispersed and the p
#     is anti-conservative - null family-wise .185 (G = 12) and .330 (G = 42).
#     v6 gives ORA a permutation null from the same block-permutation stream
#     as the other calibrated statistics.
#
# (3) globaltest was corrected by Holm over rank p-values whose granularity is
#     1/(B+1). At B = 499 and G = 42 the smallest achievable Holm-adjusted p is
#     42/500 = .084, so rejection was ARITHMETICALLY IMPOSSIBLE - which is why
#     v5 reported globaltest power of exactly 0.000 and "no N suffices" at
#     G = 42. That is the same discreteness defect that zeroed the combined
#     DMSA engine on your K = 42 run, now in a competitor. v6 corrects
#     globaltest by maxT on its own permutation stream (floor 1/(B+1) = .002).
#
# FAIRNESS NOTE. (2) and (3) both make competitors STRONGER, and (1) makes
# camera weaker only by removing an error rate it should never have had. The
# net effect of this fix is to make the comparison harder for DMSA, not easier.
#
# BEFORE RUNNING:
#   1. keep ~/dmsa_nsweep_v5/nsweep_v5_raw_empirical.csv where it is
#   2. source("dmsa_nsweep_v6_competitors_mac.R")
# RUNTIME: much shorter than v5 - no DMSA permutations. Expect ~35-50 min at
# NSIM = 100, B = 499 on 16 cores. Saves after every row; resumes on re-source.
# OUTPUT (~/dmsa_nsweep_v6/): raw / summary / N80 CSVs, DMSA rows merged in.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(parallel); library(limma)})

NSIM  <- 100L
B     <- 499L
CORES <- max(1L, detectCores() - 2L)
if (.Platform$OS.type == "windows") CORES <- 1L
V5RAW  <- path.expand("~/dmsa_nsweep_v5/nsweep_v5_raw_empirical.csv")
OUTDIR <- path.expand("~/dmsa_nsweep_v6"); dir.create(OUTDIR, showWarnings = FALSE)
RAW    <- file.path(OUTDIR, "nsweep_v6_competitors_raw.csv")
NBG    <- 60L
KPOOL  <- c(1,1,2,2,2,3,3,3,4,4,4,5,5,6,7,7,8,8,9,10,12,16)
rmvn <- function(n, K, r = .6) { z <- matrix(rnorm(n*K), n); f <- rnorm(n)
  sqrt(r)*f + sqrt(1-r)*z }

## ---------------------------------------------------------------- one cell
## ALL data generation happens before any inference, so the dataset a
## competitor sees here is the dataset DMSA saw in v5 for the same seed.
one <- function(n, G, beta, shape, p_plus_ratio, seed) {
  set.seed(seed)
  n <- 2L*(n %/% 2L); cid <- rep(seq_len(n/2), each = 2)
  K <- sample(KPOOL, G, TRUE)
  gene <- rep(paste0("g", seq_len(G)), K); P <- length(gene)
  y <- sqrt(.25)*rep(rnorm(n/2), each = 2) + sqrt(.75)*rnorm(n)
  dat <- data.frame(y = y, cv1 = rnorm(n), cv2 = rnorm(n))
  M <- matrix(NA_real_, n, P); st <- c(0L, cumsum(K))
  for (g in seq_len(G)) M[, (st[g]+1):st[g+1]] <- rmvn(n, K[g])
  M <- M + .3*rep(rnorm(n/2), each = 2)
  d_true <- ifelse(runif(P) < p_plus_ratio, 1, -1)
  target <- "g1"; units <- gene
  if (shape == "coherent") {
    if (beta > 0) for (j in (st[1]+1):st[2]) M[, j] <- M[, j] + beta*y*d_true[j]
  } else if (shape == "block") {
    if (beta > 0) { idx <- (st[1]+1):st[2]
      flip <- rep(1, length(idx)); flip[seq_len(round(.2*length(idx)))] <- -1
      for (k in seq_along(idx))
        M[, idx[k]] <- M[, idx[k]] + .85*beta*y*d_true[idx[k]]*flip[k] }
  } else if (shape == "sys3") {
    if (beta > 0) for (j in seq_len(P)) if (runif(1) < .5)
      M[, j] <- M[, j] + .55*beta*y*d_true[j]
    units <- rep("SYS1", P)
    for (s in 2:3) {
      K2 <- sample(KPOOL, G, TRUE); P2 <- sum(K2)
      M2 <- matrix(NA_real_, n, P2); s2 <- c(0L, cumsum(K2))
      for (g in seq_len(G)) M2[, (s2[g]+1):s2[g+1]] <- rmvn(n, K2[g])
      M2 <- M2 + .3*rep(rnorm(n/2), each = 2)
      M <- cbind(M, M2); units <- c(units, rep(paste0("SYS", s), P2))
    }
    P <- length(units); target <- "SYS1"
  }
  colnames(M) <- sprintf("cg%05d", seq_len(P))
  ug <- unique(units); Gn <- length(ug); idxL <- split(seq_len(P), units)[ug]
  MB <- cbind(M, matrix(rnorm(n*NBG), n) + .3*rep(rnorm(n/2), each = 2))
  X  <- model.matrix(~ y + cv1 + cv2, dat)
  hit <- function(pv) { v <- pv[match(target, ug)]
    list(hit = isTRUE(is.finite(v) && v < .05),
         false = isTRUE(any(pv[ug != target] < .05, na.rm = TRUE))) }
  out <- list()

  ## ---- camera / fry: native inference, ESTIMATED inter-gene correlation ---
  Yt <- t(scale(MB))
  addc <- function(nm, tb) { if (is.null(tb)) return(invisible())
    p <- tb$PValue[match(ug, rownames(tb))]
    h <- hit(p.adjust(p, "holm")); hb <- hit(p.adjust(p, "BH"))
    out[[nm]] <<- data.table(engine = nm, hit = h$hit, false = h$false,
                             hit_bh = hb$hit) }
  addc("camera", tryCatch(camera(Yt, idxL, X, contrast = 2,
                                 inter.gene.cor = NA), error = function(e) NULL))
  addc("fry",    tryCatch(fry(Yt, idxL, X, contrast = 2), error = function(e) NULL))

  ## ---- ONE shared block-permutation stream for globaltest and ORA ---------
  C <- X[, c(1,3,4)]; PZ <- C %*% solve(crossprod(C), t(C))
  yr <- as.numeric(y - PZ %*% y); Mres <- scale(M)
  gt_stat <- function(v) { s2 <- as.numeric(crossprod(Mres, v))^2
    vapply(idxL, function(j) sum(s2[j])/(length(j)*sum(v^2)), numeric(1)) }
  XtXi <- solve(crossprod(X))
  ora_stat <- function(yv) {                       # -log10 Fisher, permutable
    Xp <- X; Xp[, 2] <- yv
    bh_ <- XtXi %*% t(Xp) %*% MB; R_ <- MB - Xp %*% bh_
    se_ <- sqrt(colSums(R_^2)/(n - ncol(Xp))*XtXi[2,2])
    pp_ <- 2*pt(-abs(bh_[2,]/se_), n - ncol(Xp))
    bg <- mean(pp_[-seq_len(P)] < .05)
    vapply(idxL, function(j) {
      k <- sum(pp_[j] < .05)
      -log10(max(stats::pbinom(k - 1, length(j), max(bg, 1e-6),
                               lower.tail = FALSE), 1e-300)) }, numeric(1)) }
  gt_o <- gt_stat(yr); or_o <- ora_stat(y)
  rows <- split(seq_len(n), cid)
  gt_n <- matrix(NA_real_, B, Gn); or_n <- matrix(NA_real_, B, Gn)
  for (b in seq_len(B)) {
    perm  <- unlist(rows[sample.int(length(rows))], use.names = FALSE)
    gt_n[b, ] <- gt_stat(yr[perm])
    or_n[b, ] <- ora_stat(y[perm])
  }
  ## maxT on each stream: no Holm, so no rank-p x family-size floor
  maxT <- function(obs, nul) {
    mx <- apply(nul, 1, max, na.rm = TRUE)
    vapply(seq_along(obs), function(g)
      (1 + sum(mx >= obs[g]))/(B + 1), numeric(1)) }
  bhp <- function(obs, nul) p.adjust(vapply(seq_along(obs), function(g)
      (1 + sum(nul[, g] >= obs[g]))/(B + 1), numeric(1)), "BH")
  h <- hit(maxT(gt_o, gt_n)); hb <- hit(bhp(gt_o, gt_n))
  out[["globaltest"]] <- data.table(engine = "globaltest", hit = h$hit,
                                    false = h$false, hit_bh = hb$hit)
  h <- hit(maxT(or_o, or_n)); hb <- hit(bhp(or_o, or_n))
  out[["ora"]] <- data.table(engine = "ora", hit = h$hit, false = h$false,
                             hit_bh = hb$hit)

  rbindlist(out)[, `:=`(n = n, G = Gn, beta = beta, shape = shape,
                        p_plus = p_plus_ratio, seed = seed)]
}

NGRID <- c(75L,150L,250L,400L,600L,900L)
GRID <- rbind(
  CJ(n=NGRID, G=c(6L,12L,42L), beta=c(.06,.10,.14), shape="coherent", p_plus_ratio=.30),
  CJ(n=NGRID, G=12L, beta=c(.10,.14), shape="coherent", p_plus_ratio=.50),
  CJ(n=NGRID, G=12L, beta=c(.10,.14), shape="block",    p_plus_ratio=.30),
  CJ(n=NGRID, G=8L,  beta=c(.06,.10,.14), shape="sys3", p_plus_ratio=.30),
  CJ(n=NGRID, G=c(12L,42L), beta=0, shape="coherent",   p_plus_ratio=.30),
  CJ(n=NGRID, G=8L,  beta=0, shape="sys3",              p_plus_ratio=.30))
prev <- if (file.exists(RAW)) fread(RAW) else NULL
done <- vapply(seq_len(nrow(GRID)), function(i) { g <- GRID[i]
  !is.null(prev) && nrow(prev[n==2L*(g$n%/%2L) & beta==g$beta & shape==g$shape &
    p_plus==g$p_plus_ratio & engine=="camera" &
    G==ifelse(g$shape=="sys3",3L,g$G)]) >= NSIM-10L }, logical(1))
cat(sprintf("v6 competitors | grid %d rows (%d done) | NSIM %d | B %d | %d cores\n",
            nrow(GRID), sum(done), NSIM, B, CORES))
acc <- if (is.null(prev)) list() else list(prev); t00 <- Sys.time()
for (i in seq_len(nrow(GRID))) {
  if (done[i]) next
  g <- GRID[i]; t0 <- Sys.time()
  rr <- mclapply(seq_len(NSIM), function(s)
    one(g$n,g$G,g$beta,g$shape,g$p_plus_ratio, seed=940000L+1000L*i+s), mc.cores=CORES)
  acc[[length(acc)+1L]] <- rbindlist(Filter(Negate(is.null), rr), fill=TRUE)
  fwrite(rbindlist(acc, fill=TRUE), RAW)
  el <- as.numeric(Sys.time()-t0, units="mins")
  cat(sprintf("%3d/%3d %-8s n=%3d G=%2d b=%.2f p+=%.2f  %.1f min (~%.0f left)\n",
      i,nrow(GRID),g$shape,g$n,g$G,g$beta,g$p_plus_ratio,el,el*sum(!done[i:nrow(GRID)])))
}

## ---- merge with the v5 DMSA rows (still valid) and summarise -------------
CMP <- rbindlist(acc, fill=TRUE)
DM  <- if (file.exists(V5RAW)) {
  d <- fread(V5RAW); d[grepl("^dmsa_", engine), .(engine, hit, false, hit_bh,
                                                  n, G, beta, shape, p_plus, seed)]
} else { cat("!! v5 raw not found - DMSA rows not merged\n"); NULL }
R <- rbindlist(list(CMP[, .(engine,hit,false,hit_bh,n,G,beta,shape,p_plus,seed)], DM),
               fill = TRUE)
SUM <- R[, .(power=mean(hit), power_bh=mean(hit_bh), false_name=mean(false),
             fw_err=mean(hit | false), sims=.N),
         by=.(shape,G,beta,p_plus,n,engine)]
fwrite(SUM, file.path(OUTDIR, "nsweep_v6_summary.csv"))
cat("\n==== REALISED FAMILY-WISE TYPE I (beta = 0) - ALL ENGINES ====\n")
cat("Every engine here should sit near .05. If camera is not, stop and tell me.\n")
print(dcast(SUM[beta==0], shape+G+n ~ engine, value.var="fw_err"), nrows=40)
cat("\n==== POWER at FWER .05 ====\n")
print(dcast(SUM[beta>0], shape+G+beta+p_plus+n ~ engine, value.var="power"), nrows=100)

n80 <- function(d) {
  if (mean(d$hit) < .02) return(data.table(N80=Inf, lo=NA_real_, hi=NA_real_,
    note="no N suffices (flat power curve)"))
  if (length(unique(d$hit)) < 2) return(data.table(
    N80=if (mean(d$hit)>=.8) as.double(min(d$n)) else Inf,
    lo=NA_real_, hi=NA_real_, note="boundary"))
  f <- suppressWarnings(glm(hit ~ log(n), binomial, d))
  sol <- function(cf) exp((log(.8/.2)-cf[1])/cf[2]); est <- sol(coef(f))
  if (!is.finite(est) || est <= 0) est <- Inf
  bs <- suppressWarnings(replicate(200, { i <- sample(nrow(d), replace=TRUE)
    ff <- try(glm(hit ~ log(n), binomial, d[i]), silent=TRUE)
    if (inherits(ff,"try-error")) NA_real_ else sol(coef(ff)) }))
  data.table(N80=est, lo=quantile(bs,.05,na.rm=TRUE), hi=quantile(bs,.95,na.rm=TRUE),
    note=if (!is.finite(est)) "no N suffices" else
         if (est > max(d$n)) "extrapolated beyond grid" else "within grid")
}
N80 <- R[beta>0, n80(.SD), by=.(shape,G,beta,p_plus,engine)]
N80[, `:=`(N80=ifelse(is.finite(N80), round(N80), Inf), lo=round(lo), hi=round(hi))]
fwrite(N80, file.path(OUTDIR, "nsweep_v6_N80.csv"))
cat("\n==== RECOMMENDED N for 80% power, ALL ENGINES ====\n")
print(N80[note=="within grid"][order(shape,G,beta,p_plus,N80)], nrows=200)
cat("\n(rows outside the grid are in the CSV; treat them as bounds, not estimates)\n")
cat(sprintf("\nDONE in %.1f min. Send back the two CSVs in %s\n",
    as.numeric(Sys.time()-t00, units="mins"), OUTDIR))
