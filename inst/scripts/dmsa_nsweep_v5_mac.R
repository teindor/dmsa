# ============================================================================
# RECOMMENDED N - v5 (Mac). FINAL. Supersedes v2/v3/v4. Source this file.
#
# WHAT CHANGED FROM v3, and why. Two corrections, both from Tsachi's review:
#
# (A) MAP ACCURACY IS NO LONGER AN ARBITRARY CONSTANT. v3 flipped a flat 15%
#     of direction calls (ACC = .85) while telling DMSA every call was 90%
#     reliable - internally inconsistent, and unjustified. The real Alpha map
#     already carries per-probe confidence. Measured on the 620 direction-
#     called probes of coverage_v4_full.csv, the confidence in the SIGNED call
#     - max(p_plus, 1 - p_plus) - has median .850, mean .837, range .655-1.00,
#     i.e. a median implied error of exactly 15%. So .85 was the right centre
#     by accident, but applied wrongly.
#     v5 instead SAMPLES each probe's confidence from that empirical
#     distribution (embedded below as a 41-point quantile grid, so the script
#     is self-contained), and then makes the map SELF-CONSISTENT: a probe that
#     claims confidence c is wrong with probability 1 - c, and DMSA is handed
#     exactly that c as its p_plus. No probe is secretly more or less reliable
#     than it says. This is the honest version of "the direction is validated
#     by cpgdirection": the validation supplies a confidence, and the
#     simulation both honours it and holds the method to it.
#     MAP_MODE lets you vary it: "empirical" (default), "perfect" (c = 1, the
#     ceiling), "degraded" (confidence shrunk 15% toward .5, a stress test).
#     Running all three quantifies how much of DMSA's advantage depends on map
#     quality - which is the question a reviewer will actually ask.
#
# (B) THE PRIMARY RULE IS THE DESIGNED ONE: a unit is named when ANY of its
#     three lens-adjusted p-values clears .05. The three lenses are three
#     views of one problem, and a hit along one path is the architecture, not
#     a loophole (T.E.-D., design decision, 15 Aug 2026). Because the lenses
#     are strongly correlated, the realised family-wise error of that rule is
#     ~.06 at a nominal .05 - a declared and measured property of the design,
#     not a defect, and NOT to be "corrected" away.
#     What this script therefore does: it scores every DMSA engine under the
#     DESIGNED rule (column `hit`), and additionally records the conservative
#     single-statistic rule `p_unit_adj` (column `hit_cons`, realised error
#     .030-.040) purely as a reference. Both are reported; the designed rule
#     is primary.
#     ONE REPORTING CONSEQUENCE, which protects the claim rather than changing
#     it: because the competitors are corrected by Holm - which with
#     correlated tests typically realises well BELOW .05 - the comparison is
#     only readable if each engine's REALISED type I is printed beside its
#     power. The summary below does exactly that, so a reader can see that
#     DMSA is being run at ~.06 and Holm-corrected engines at ~.02-.03, and
#     judge the power difference in that light.
#
# EVERYTHING ELSE AS v3: seven engines on identical data, one p per unit each,
# FWER primary (maxT for DMSA, Holm for competitors) with BH reported
# alongside, native inference for every competitor, curve-fit N80 with
# bootstrap CI, and cells where no N suffices reported as Inf.
#
# BEFORE RUNNING (once):
#   1. R CMD INSTALL dmsa_1.5.0.tar.gz     # RESTART R first if dmsa is loaded
#   2. source("dmsa_nsweep_v5_mac.R")
# Stop v2/v3/v4 if running; v5 writes to its own directory.
#
# RUNTIME - MEASURED, not guessed. Timed on one core at B = 999: 27 s/sim at
# G = 6, 48 s at G = 12, 188 s at G = 42; runtime is driven almost entirely by
# the number of PROBES, barely at all by n (n = 75 and n = 900 differ by ~12%).
# The full 114-row grid therefore costs about 480 CORE-hours at NSIM = 200,
# B = 999 - i.e. roughly 40 h wall-clock on 12 effective cores, NOT the 5-8 h
# an earlier draft of this header claimed. Cost is linear in both NSIM and B:
#     NSIM 200 / B 999  ~ 40 h     NSIM 100 / B 499  ~ 10 h   <- DEFAULT
#     NSIM 100 / B 999  ~ 20 h     NSIM  60 / B 499  ~  6 h   <- first pass
# B = 499 is safe: the discreteness floor is gone (continuous fusion), so B now
# sets only null RESOLUTION - 1/500 = .002, far finer than a .05 threshold -
# and the K = 42 case was verified working at B = 399. NSIM = 100 is enough
# because N80 comes from a curve FIT that pools all six n-points (600 sims per
# cell), not from any single cell.
# The script CALIBRATES ON YOUR MACHINE before it starts and prints a projected
# total; if you do not like the number, Ctrl-C there and lower NSIM.
# Saves after every grid row and resumes on re-source. For the map-quality
# sensitivity, re-source twice more with MAP_MODE <- "perfect" then "degraded".
# OUTPUT (~/dmsa_nsweep_v5/): raw / summary / N80 CSVs - send all back.
# ============================================================================
suppressPackageStartupMessages({
  library(data.table); library(dmsa); library(parallel); library(limma)})
stopifnot(utils::packageVersion("dmsa") >= "1.5.0")

NSIM     <- 100L                    # 60 for a ~6 h first pass; 200 costs 2x
B        <- 499L                    # 999 costs 2x; 499 is ample (see header)
MAP_MODE <- "empirical"             # "empirical" | "perfect" | "degraded"
CORES    <- max(1L, detectCores() - 2L)
if (.Platform$OS.type == "windows") CORES <- 1L
OUTDIR <- path.expand("~/dmsa_nsweep_v5"); dir.create(OUTDIR, showWarnings = FALSE)
RAW <- file.path(OUTDIR, sprintf("nsweep_v5_raw_%s.csv", MAP_MODE))
NBG <- 60L
## Empirical confidence in the signed direction call, measured on the 620
## direction-called probes of the Alpha panel (median .850, mean .837).
CONF_Q <- c(0.655,0.6704,0.6733,0.6844,0.7136,0.7216,0.7358,0.7426,0.76,0.7707,
            0.7855,0.7995,0.8068,0.8148,0.8216,0.8268,0.8372,0.85,0.85,0.85,
            0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.8546,
            0.9183,0.9697,1,1,1,1,1,1)
draw_conf <- function(P) {
  cf <- sample(CONF_Q, P, replace = TRUE)
  switch(MAP_MODE,
    empirical = cf,
    perfect   = rep(1, P),
    degraded  = 0.5 + 0.85*(cf - 0.5))     # shrink 15% toward a coin flip
}
KPOOL <- c(1,1,2,2,2,3,3,3,4,4,4,5,5,6,7,7,8,8,9,10,12,16)
rmvn <- function(n, K, r = .6) { z <- matrix(rnorm(n*K), n); f <- rnorm(n)
  sqrt(r)*f + sqrt(1-r)*z }
ora_p <- function(ps, pb, a = .05) stats::fisher.test(
  matrix(c(sum(ps < a), sum(ps >= a), sum(pb < a), sum(pb >= a)), 2),
  alternative = "greater")$p.value

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
      M <- cbind(M, M2)
      d_true <- c(d_true, ifelse(runif(P2) < p_plus_ratio, 1, -1))
      units <- c(units, rep(paste0("SYS", s), P2))
    }
    P <- length(units); target <- "SYS1"
  }
  colnames(M) <- sprintf("cg%05d", seq_len(P))
  ug <- unique(units); Gn <- length(ug); idxL <- split(seq_len(P), units)[ug]
  ## ---- SELF-CONSISTENT map: confidence c, wrong with probability 1 - c ----
  conf  <- draw_conf(P)
  wrong <- runif(P) > conf                       # <- the only source of error
  d_map <- d_true * ifelse(wrong, -1, 1)
  ## DMSA is handed exactly the confidence the map claims for that probe
  p_plus_j <- ifelse(d_map > 0, conf, 1 - conf)
  al <- dmsa_align(data.frame(cpg = colnames(M), d = d_map, p_plus = p_plus_j),
                   genes = units, level = "gene")
  hit <- function(pv) { v <- pv[match(target, ug)]
    list(hit = isTRUE(is.finite(v) && v < .05),
         false = isTRUE(any(pv[ug != target] < .05, na.rm = TRUE))) }
  out <- list()
  for (eng in c("flat", "reliability", "combined")) {
    r <- tryCatch(suppressWarnings(dmsa_triangulate(
           M, dat, c("y","cv1","cv2"), "y", units, al, block = cid, B = B,
           correction = "maxT", weighting = eng, seed = seed)),
         error = function(e) NULL)
    if (is.null(r)) next
    ## DESIGNED rule: any lens may carry the unit (three views, one problem)
    best <- pmin(r$p_coherence_adj, r$p_composite_adj, r$p_diffuse_adj,
                 na.rm = TRUE)[match(ug, r$unit)]
    h  <- hit(best)
    ## conservative single-statistic reference (not the primary rule)
    hc <- hit(r$p_unit_adj[match(ug, r$unit)])
    hb <- hit(p.adjust(r$p_unit[match(ug, r$unit)], "BH"))
    out[[eng]] <- data.table(engine = paste0("dmsa_", eng), hit = h$hit,
                             false = h$false, hit_cons = hc$hit,
                             false_cons = hc$false, hit_bh = hb$hit)
  }
  MB <- cbind(M, matrix(rnorm(n*NBG), n) + .3*rep(rnorm(n/2), each = 2))
  Yt <- t(scale(MB)); X <- model.matrix(~ y + cv1 + cv2, dat)
  addc <- function(nm, tb) { if (is.null(tb)) return(invisible())
    p <- tb$PValue[match(ug, rownames(tb))]
    h <- hit(p.adjust(p, "holm")); hb <- hit(p.adjust(p, "BH"))
    out[[nm]] <<- data.table(engine = nm, hit = h$hit, false = h$false,
                             hit_cons = h$hit, false_cons = h$false,
                             hit_bh = hb$hit) }
  addc("camera", tryCatch(camera(Yt, idxL, X, contrast = 2), error = function(e) NULL))
  addc("fry",    tryCatch(fry(Yt, idxL, X, contrast = 2),    error = function(e) NULL))
  C <- X[, c(1,3,4)]; PZ <- C %*% solve(crossprod(C), t(C))
  yr <- as.numeric(y - PZ %*% y); Mres <- scale(M)
  qp <- function(v) { s2 <- as.numeric(crossprod(Mres, v))^2
    vapply(idxL, function(j) sum(s2[j])/(length(j)*sum(v^2)), numeric(1)) }
  qo <- qp(yr); rows <- split(seq_len(n), cid); qn <- matrix(NA_real_, B, Gn)
  for (b in seq_len(B)) qn[b, ] <-
    qp(yr[unlist(rows[sample.int(length(rows))], use.names = FALSE)])
  gp <- vapply(seq_len(Gn), function(g) (1+sum(qn[,g] >= qo[g]))/(B+1), numeric(1))
  h <- hit(p.adjust(gp,"holm")); hb <- hit(p.adjust(gp,"BH"))
  out[["globaltest"]] <- data.table(engine="globaltest", hit=h$hit, false=h$false,
                                    hit_cons=h$hit, false_cons=h$false, hit_bh=hb$hit)
  XtXi <- solve(crossprod(X)); bh_ <- XtXi %*% t(X) %*% MB
  R_ <- MB - X %*% bh_; se_ <- sqrt(colSums(R_^2)/(n-ncol(X))*XtXi[2,2])
  pp_ <- 2*pt(-abs(bh_[2,]/se_), n-ncol(X))
  op <- vapply(idxL, function(j) ora_p(pp_[j], pp_[-seq_len(P)]), numeric(1))
  h <- hit(p.adjust(op,"holm")); hb <- hit(p.adjust(op,"BH"))
  out[["ora"]] <- data.table(engine="ora", hit=h$hit, false=h$false,
                             hit_cons=h$hit, false_cons=h$false, hit_bh=hb$hit)
  rbindlist(out)[, `:=`(n=n, G=Gn, beta=beta, shape=shape, p_plus=p_plus_ratio,
                        map=MAP_MODE, mean_conf=mean(conf), seed=seed)]
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
    p_plus==g$p_plus_ratio & engine=="dmsa_flat" &
    G==ifelse(g$shape=="sys3",3L,g$G)]) >= NSIM-10L }, logical(1))
cat(sprintf("grid %d rows (%d done) | NSIM %d | B %d | %d cores | map=%s\n",
            nrow(GRID), sum(done), NSIM, B, CORES, MAP_MODE))

## ---- CALIBRATE ON THIS MACHINE, then project the total -------------------
## Runtime is driven by probe count, so cost per sim is ~linear in the unit
## count (x3 for sys3, which carries three systems). One cheap timed sim at
## B = 99 on the largest cell fixes the constant; everything else scales.
if (any(!done)) {
  .Bkeep <- B; B <- 99L; .calG <- max(GRID$G)
  .t <- Sys.time(); invisible(one(150L, .calG, .10, "coherent", .30, seed = 4242L))
  .sec <- as.numeric(Sys.time() - .t, units = "secs"); B <- .Bkeep
  .rel <- ifelse(GRID$shape == "sys3", 3 * GRID$G, GRID$G) / .calG
  .hrs <- NSIM * sum(.rel[!done]) * .sec * (B / 99) / CORES / 3600
  .fmt <- function(h) if (h < 1) sprintf("%.0f min", h * 60) else sprintf("%.1f h", h)
  cat(sprintf("calibration: %.1f s/sim at G=%d, B=99 -> PROJECTED TOTAL %s",
              .sec, .calG, .fmt(.hrs)))
  cat(sprintf("  (%s at NSIM=200, %s at NSIM=60)\n",
              .fmt(.hrs * 200 / NSIM), .fmt(.hrs * 60 / NSIM)))
  cat("Ctrl-C now if that is too long; lower NSIM and re-source. Resumes cleanly.\n")
  Sys.sleep(10)
}
acc <- if (is.null(prev)) list() else list(prev); t00 <- Sys.time()
for (i in seq_len(nrow(GRID))) {
  if (done[i]) next
  g <- GRID[i]; t0 <- Sys.time()
  rows <- mclapply(seq_len(NSIM), function(s)
    one(g$n,g$G,g$beta,g$shape,g$p_plus_ratio, seed=940000L+1000L*i+s), mc.cores=CORES)
  acc[[length(acc)+1L]] <- rbindlist(Filter(Negate(is.null), rows), fill=TRUE)
  fwrite(rbindlist(acc, fill=TRUE), RAW)
  el <- as.numeric(Sys.time()-t0, units="mins")
  cat(sprintf("%2d/%2d %-8s n=%3d G=%2d b=%.2f p+=%.2f  %.1f min (~%.0f left)\n",
      i,nrow(GRID),g$shape,g$n,g$G,g$beta,g$p_plus_ratio,el,el*sum(!done[i:nrow(GRID)])))
}

R <- rbindlist(acc, fill=TRUE)
SUM <- R[, .(power=mean(hit), power_cons=mean(hit_cons), power_bh=mean(hit_bh),
             false_name=mean(false), false_cons=mean(false_cons),
             ## realised FAMILY-WISE error: ANY unit named. Under beta = 0 the
             ## target is just another null unit, so naming IT is also an error;
             ## P(hit OR false), not mean(false), is the quantity to compare.
             fw_err      = mean(hit      | false),
             fw_err_cons = mean(hit_cons | false_cons),
             sims=.N, mean_conf=round(mean(mean_conf),3)),
         by=.(shape,G,beta,p_plus,n,engine)]
fwrite(SUM, file.path(OUTDIR, sprintf("nsweep_v5_summary_%s.csv", MAP_MODE)))
cat("\n==== REALISED FAMILY-WISE TYPE I (beta = 0) ====\n")
cat("DMSA rows are the DESIGNED any-lens rule (expected ~.06 by design);\n")
cat("competitor rows are Holm (typically well below .05 with correlated tests).\n")
cat("Read every power number below against its engine's realised rate here.\n")
print(dcast(SUM[beta==0], shape+G+n ~ engine, value.var="fw_err"), nrows=40)
cat("\n-- same cells under the conservative single-statistic rule (reference) --\n")
print(dcast(SUM[beta==0], shape+G+n ~ engine, value.var="fw_err_cons"), nrows=40)
cat("\n==== POWER at FWER .05 ====\n")
print(dcast(SUM[beta>0], shape+G+beta+p_plus+n ~ engine, value.var="power"), nrows=100)

n80 <- function(d) {
  if (mean(d$hit) < .02) return(data.table(N80=Inf, lo=NA_real_, hi=NA_real_,
    note="no N suffices (flat power curve)"))
  if (length(unique(d$hit)) < 2) return(data.table(
    N80=if (mean(d$hit)>=.8) as.double(min(d$n)) else Inf,   # as.double: data.table
    lo=NA_real_, hi=NA_real_, note="boundary"))              # needs one type per group
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
fwrite(N80, file.path(OUTDIR, sprintf("nsweep_v5_N80_%s.csv", MAP_MODE)))
cat("\n==== RECOMMENDED N for 80% power, ALL ENGINES (90% bootstrap CI) ====\n")
print(N80[order(shape,G,beta,p_plus,N80)], nrows=200)
cat(sprintf("\nDONE in %.1f min (map = %s). Send back the three CSVs in %s\n",
    as.numeric(Sys.time()-t00, units="mins"), MAP_MODE, OUTDIR))
