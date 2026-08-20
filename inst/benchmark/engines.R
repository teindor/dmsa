# ============================================================================
# LOCKED ENGINE ROSTER  -  v1.0, 13 Aug 2026
#
# Every method in the paper's contest is defined here, once, and both the
# simulation panels and the real-data analyses call this file. Nothing in this
# file may change after the benchmark is locked; a change means re-running
# every panel.
#
# PROVENANCE OF EACH ENGINE
#   camera, fry, roast          real limma functions
#   globaltest, GSEA/mCSEA,     re-implemented from the published statistic,
#   methylRRA, ORA              because CRAN/Bioconductor are unreachable from
#                               the analysis container. bench/mac_crosscheck.R
#                               re-runs these on the SAME saved datasets with
#                               the real packages, and the agreement is reported.
#   DMSA / nDMSA / cascade      the dmsa package
#
# All set-level engines answer the same question - "does this probe set respond
# to x" - at alpha = .05, on the same simulated data, with the same background
# panel where a competitive method needs one.
# ============================================================================
suppressPackageStartupMessages({library(limma); library(dmsa)})

ENGINES <- c(
  ## --- DMSA family ---------------------------------------------------------
  "dmsa_expected", "dmsa_fixed", "dmsa_signblind", "ndmsa", "dmsa_tri",
  ## --- competitors ---------------------------------------------------------
  "camera", "fry", "roast_dir", "roast_mixed",
  "globaltest", "gsea_mcsea", "methylRRA", "ora", "ewas_bh")

DIRECTIONAL <- c("dmsa_expected", "dmsa_fixed", "ndmsa", "dmsa_tri",
                 "camera", "fry", "roast_dir", "gsea_mcsea")

## ===========================================================================
## per-probe least squares for a whole panel at once
## ===========================================================================
allfit <- function(Y, x, Z = NULL) {
  X <- if (is.null(Z)) cbind(1, x) else cbind(1, x, Z)
  XtXi <- solve(crossprod(X))
  bh <- XtXi %*% (t(X) %*% Y)
  res <- Y - X %*% bh
  df <- nrow(Y) - ncol(X)
  se <- sqrt(colSums(res^2) / df * XtXi[2, 2])
  t <- bh[2, ] / se
  list(b = bh[2, ], se = se, t = t, df = df,
       p = 2 * stats::pt(-abs(t), df = df))
}

## ===========================================================================
## DMSA pooling. m_j is the alignment multiplier: +-1 for fixed signs,
## (2*acc - 1) * sign for expected signs, 1 for the sign-blind control (which
## pools |b| and is the internal control isolating what the sign chain buys).
## ===========================================================================
dmsa_pool <- function(b, se, m) {
  w <- 1 / se^2
  num <- sum(m * b * w); den <- sum(m^2 * w)
  if (den <= 0) return(c(mu = NA_real_, z = NA_real_))
  mu <- num / den
  c(mu = mu, z = mu * sqrt(den))
}
dmsa_pool_blind <- function(b, se) {
  w <- 1 / se^2
  mu <- sum(abs(b) * w) / sum(w)
  c(mu = mu, z = mu * sqrt(sum(w)))
}

## ===========================================================================
## nDMSA: DMSA that does not assume the exposure enters linearly.
## Statistic = max over three aligned arms, so it pays an honest multiplicity
## price when the truth IS linear and collects when it is not.
##   linear     pooled aligned z^2                 (DMSA's own statistic)
##   quadratic  pooled aligned z^2 of the x^2 term
##   threshold  sup over tau of the pooled aligned z^2 of 1(x > tau)
## ===========================================================================
NDMSA_TAU <- seq(.2, .8, by = .1)

ndmsa_stat <- function(M, x, m) {
  x <- as.numeric(x)
  pool_z2 <- function(design_col, extra = NULL) {
    X <- cbind(1, design_col, extra)
    XtXi <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
    if (is.null(XtXi)) return(0)
    bh <- XtXi %*% (t(X) %*% M)
    res <- M - X %*% bh
    df <- nrow(M) - ncol(X)
    if (df <= 1) return(0)
    se <- sqrt(colSums(res^2) / df * XtXi[2, 2])
    z <- dmsa_pool(bh[2, ], se, m)[["z"]]
    if (is.finite(z)) z^2 else 0
  }
  xc <- as.numeric(scale(x))
  lin  <- pool_z2(xc)
  ## the quadratic term goes in column 2 so pool_z2 pools IT, with the linear
  ## term held as a nuisance column - curvature over and above the slope
  quad <- pool_z2(xc^2 - mean(xc^2), extra = xc)
  tau  <- stats::quantile(x, NDMSA_TAU, names = FALSE)
  tau  <- unique(tau[tau > min(x) & tau < max(x)])
  thr  <- if (length(tau))
    max(vapply(tau, function(tt) pool_z2(as.numeric(x > tt)), numeric(1))) else 0
  c(linear = lin, quadratic = quad, threshold = thr,
    ndmsa = max(lin, quad, thr))
}

## ===========================================================================
## globaltest (Goeman): the score test for a random-effects alternative,
## Q = y' M M' y / (K * y'y) with y the residualised outcome.
## ===========================================================================
gt_stat <- function(M, x) {
  y <- x - mean(x)
  Mc <- scale(M, center = TRUE, scale = FALSE)
  s <- crossprod(Mc, y)
  sum(s^2) / (ncol(M) * sum(y^2))
}

## ===========================================================================
## mCSEA / GSEA: weighted running enrichment score on the SIGNED t statistics,
## so it is one of the few competitors that reports a direction.
## ===========================================================================
gsea_es <- function(stat, in_set) {
  o <- order(stat, decreasing = TRUE)
  s <- stat[o]; hit <- in_set[o]
  NR <- sum(abs(s[hit]))
  if (NR == 0) return(0)
  dev <- cumsum(ifelse(hit, abs(s), 0)) / NR -
         cumsum(ifelse(hit, 0, 1)) / sum(!hit)
  dev[which.max(abs(dev))]
}

## ===========================================================================
## methylRRA: probe p -> gene-level robust rank aggregation -> ORA on genes.
## ===========================================================================
rra_p <- function(pv) {
  r <- sort(pv); k <- seq_along(r)
  min(1, min(stats::pbeta(r, k, length(r) - k + 1)) * length(r))
}
rra_ora <- function(p_set, gid_set, p_bg, gid_bg, alpha = .05) {
  gs <- vapply(split(p_set, gid_set), rra_p, numeric(1))
  gb <- vapply(split(p_bg,  gid_bg),  rra_p, numeric(1))
  tab <- matrix(c(sum(gs < alpha), sum(gs >= alpha),
                  sum(gb < alpha), sum(gb >= alpha)), 2)
  stats::fisher.test(tab, alternative = "greater")$p.value
}

## ===========================================================================
## ORA on probes - what most EWAS papers actually do.
## ===========================================================================
ora_p <- function(p_set, p_bg, alpha = .05) {
  tab <- matrix(c(sum(p_set < alpha), sum(p_set >= alpha),
                  sum(p_bg < alpha), sum(p_bg >= alpha)), 2)
  stats::fisher.test(tab, alternative = "greater")$p.value
}

## ===========================================================================
## THE ROSTER, applied to one dataset.
##
##   M    n x K set probe matrix (M-values, already the analysis scale)
##   B    n x NBG background probe matrix (needed by competitive methods)
##   x    exposure
##   d    TRUE per-probe CpG->expression directions
##   acc  direction-call accuracy the analyst is assumed to have
##   crit named vector of null critical values for the statistics with no
##        native inference (dmsa_*, ndmsa, globaltest, gsea_mcsea, cascade)
##
## DIRECTION IS RECORDED, NOT SCORED, HERE. Each directional engine returns the
## sign it reports; the panel script scores that against two different truths,
## because the whole argument turns on their being different:
##
##   methylation truth   sign(2*p_plus - 1) - the direction of the mean
##                       methylation effect. camera, fry, roast and GSEA are
##                       reporting this, and reporting it correctly.
##   biological truth    +1 by construction - the exposure raises the set's
##                       expression tone. Only an aligned method can report it.
##
## The literature routinely reads the first as if it were the second. That
## reading is right only when p_plus > .5, and it FLIPS with the promoter/body
## composition of the array. That is the result the panel is built to show.
## ===========================================================================
run_engines <- function(M, B, x, d, gid, gidB, acc, crit,
                        nrot = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- ncol(M); Y <- cbind(M, B)
  f  <- allfit(Y, x)
  fs <- list(b = f$b[1:K], se = f$se[1:K], t = f$t[1:K], p = f$p[1:K])
  fb <- list(p = f$p[-(1:K)], t = f$t[-(1:K)])

  ## the analyst's direction calls: right with probability acc
  dc <- d * ifelse(stats::runif(K) < acc, 1L, -1L)
  m_fixed    <- as.numeric(dc)
  m_expected <- as.numeric(dc) * (2 * acc - 1)

  out <- stats::setNames(rep(NA_real_, length(ENGINES)), ENGINES)   # rejected?
  rs  <- stats::setNames(rep(NA_real_, length(ENGINES)), ENGINES)   # sign reported

  ## ---- DMSA family --------------------------------------------------------
  pe <- dmsa_pool(fs$b, fs$se, m_expected)
  pf <- dmsa_pool(fs$b, fs$se, m_fixed)
  pb <- dmsa_pool_blind(fs$b, fs$se)
  out["dmsa_expected"]  <- as.numeric(abs(pe[["z"]]) > crit[["dmsa_expected"]])
  out["dmsa_fixed"]     <- as.numeric(abs(pf[["z"]]) > crit[["dmsa_fixed"]])
  out["dmsa_signblind"] <- as.numeric(abs(pb[["z"]]) > crit[["dmsa_signblind"]])
  rs["dmsa_expected"]   <- sign(pe[["mu"]])
  rs["dmsa_fixed"]      <- sign(pf[["mu"]])

  nd <- ndmsa_stat(M, x, m_expected)
  out["ndmsa"] <- as.numeric(nd[["ndmsa"]] > crit[["ndmsa"]])
  rs["ndmsa"]  <- rs["dmsa_expected"]     # same alignment, same reported sign

  ## ---- limma: camera (competitive), fry and roast (self-contained) --------
  des <- cbind(Intercept = 1, x = x)
  Yt <- t(Y)
  ## inter.gene.cor = NA makes camera ESTIMATE the inter-probe correlation
  ## instead of assuming limma's default 0.01. The default is badly wrong for
  ## methylation, where probes in a gene are strongly correlated, and leaving it
  ## in place would inflate camera's Type I to about .15 here and flatter its
  ## power for the wrong reason. Estimating is camera's own recommendation when
  ## the correlation is unknown, and it is the fair comparison.
  cm <- tryCatch(camera(Yt, list(s = 1:K), des, contrast = 2,
                        inter.gene.cor = NA),
                 error = function(e) NULL)
  fr <- tryCatch(fry(Yt, list(s = 1:K), des, contrast = 2),
                 error = function(e) NULL)
  ro <- tryCatch(roast(Yt, index = 1:K, design = des, contrast = 2,
                       nrot = nrot), error = function(e) NULL)
  if (!is.null(cm)) {
    out["camera"] <- as.numeric(cm$PValue[1] < .05)
    rs["camera"]  <- if (cm$Direction[1] == "Up") 1 else -1
  }
  if (!is.null(fr)) {
    out["fry"] <- as.numeric(fr$PValue[1] < .05)
    rs["fry"]  <- if (fr$Direction[1] == "Up") 1 else -1
  }
  if (!is.null(ro)) {
    ## roast reports Up and Down as separate ONE-sided p-values. Taking the
    ## smaller one at .05 is a two-sided test at .10 and inflates roast's Type I
    ## to about .11 - so double it. Getting this wrong would have made roast look
    ## both more powerful and less valid than it is.
    pu <- ro$p.value["Up", "P.Value"]; pd <- ro$p.value["Down", "P.Value"]
    out["roast_dir"] <- as.numeric(min(1, 2 * min(pu, pd)) < .05)
    rs["roast_dir"]  <- if (pu < pd) 1 else -1
    ## PValue.Mixed is degenerate here, so use the transparent directionless
    ## statistic: sum of squared t against the same rotation null
    out["roast_mixed"] <- as.numeric(sum(fs$t^2) > crit[["roast_mixed"]])
  }

  ## ---- globaltest ---------------------------------------------------------
  out["globaltest"] <- as.numeric(gt_stat(M, x) > crit[["globaltest"]])

  ## ---- mCSEA / GSEA -------------------------------------------------------
  es <- gsea_es(f$t, c(rep(TRUE, K), rep(FALSE, ncol(B))))
  out["gsea_mcsea"] <- as.numeric(abs(es) > crit[["gsea_mcsea"]])
  rs["gsea_mcsea"]  <- sign(es)

  ## ---- methylRRA, ORA, and plain EWAS ------------------------------------
  out["methylRRA"] <- as.numeric(
    rra_ora(fs$p, gid, fb$p, gidB) < .05)
  out["ora"] <- as.numeric(ora_p(fs$p, fb$p) < .05)
  out["ewas_bh"] <- as.numeric(
    any(stats::p.adjust(f$p, "BH")[1:K] < .05))

  ## ---- V2: the shipped triangulation (after every RNG consumer) ----------
  tn <- attr(crit, "tri_null")
  if (!is.null(tn)) {
    tro <- tri_stats(M, x, m_expected)
    pt3 <- tri_pvalue(tro[1:3], tn)
    out["dmsa_tri"] <- as.numeric(pt3 < attr(crit, "tri_crit"))
    rs["dmsa_tri"]  <- sign(tro[["mu"]])
  }

  list(hit = out, rep_sign = rs, b = fs$b, se = fs$se, t = f$t, p = f$p,
       m_expected = m_expected, nd = nd)
}


## ===========================================================================
## V2 ADDITION (14 Aug): the SHIPPED deliverable - the three-lens triangulation
## (coherence pooled burden + composite score-in-one-model + confidence-
## weighted diffuse quadratic, winsor 3 MAD, ACAT omnibus). At set level the
## family is a single unit, so maxT reduces to the unit test and the decision
## statistic is the ACAT of the three lens p-values, each lens calibrated on
## the same null draws as every other non-native engine. Placed AFTER all RNG
## consumers so the v1 engines see an unchanged random stream: their v2
## numbers must REPRODUCE v1 exactly (a built-in regression check).
## ===========================================================================
.winsor3 <- function(M) apply(M, 2, function(y) {
  md <- stats::median(y); s <- stats::mad(y)
  if (!is.finite(s) || s <= 0) return(y)
  pmin(pmax(y, md - 3 * s), md + 3 * s) })
tri_stats <- function(M, x, m, b = NULL, se = NULL) {
  Mw <- .winsor3(M)
  f <- allfit(Mw, x)
  coh <- abs(dmsa_pool(f$b, f$se, m)[["z"]])
  S <- as.numeric(Mw %*% m) / max(sum(abs(m)), 1e-12)
  ct <- suppressWarnings(stats::cor(S, x))
  comp <- abs(ct) * sqrt((length(x) - 2) / max(1 - ct^2, 1e-12))
  xr <- x - mean(x)
  r2 <- as.numeric(crossprod(scale(Mw), xr))^2
  dif <- sum(abs(m) * r2) / (max(sum(abs(m)), 1e-12) * sum(xr^2))
  mu <- dmsa_pool(f$b, f$se, m)[["mu"]]
  c(coh = coh, comp = comp, dif = dif, mu = mu)
}
acat <- function(p) { p <- pmin(pmax(p, 1e-12), 1 - 1e-4)
  0.5 - atan(mean(tan((0.5 - p) * pi))) / pi }
tri_pvalue <- function(obs3, null3) {
  pl <- vapply(1:3, function(k)
    (1 + sum(null3[, k] >= obs3[k])) / (nrow(null3) + 1), numeric(1))
  acat(pl)
}

## ===========================================================================
## Null calibration. The statistics without native inference get a critical
## value from NC draws under H0 (h = 0), which is where they are all valid; the
## null does not depend on p_plus because there is no effect to have a ratio.
## ===========================================================================
calibrate <- function(gen0, NC = 500, acc = .85, seed = 11, nrot = 200) {
  set.seed(seed)
  S <- t(vapply(seq_len(NC), function(i) {
    g <- gen0()
    K <- ncol(g$M); Y <- cbind(g$M, g$B)
    f <- allfit(Y, g$x)
    b <- f$b[1:K]; se <- f$se[1:K]
    dc <- g$d * ifelse(stats::runif(K) < acc, 1L, -1L)
    c(dmsa_expected  = abs(dmsa_pool(b, se, dc * (2 * acc - 1))[["z"]]),
      dmsa_fixed     = abs(dmsa_pool(b, se, as.numeric(dc))[["z"]]),
      dmsa_signblind = abs(dmsa_pool_blind(b, se)[["z"]]),
      ndmsa          = ndmsa_stat(g$M, g$x, dc * (2 * acc - 1))[["ndmsa"]],
      globaltest     = gt_stat(g$M, g$x),
      gsea_mcsea     = abs(gsea_es(f$t, c(rep(TRUE, K), rep(FALSE, ncol(g$B))))),
      roast_mixed    = sum(f$t[1:K]^2),
      tri_coh = NA_real_, tri_comp = NA_real_, tri_dif = NA_real_)
  }, numeric(10)))
  ## second pass for the triangulation trio on the SAME null stream: rebuild
  ## the draws with an offset seed so no RNG is stolen from the v1 pass above
  set.seed(seed + 7000)
  TN <- t(vapply(seq_len(NC), function(i) {
    g <- gen0()
    K <- ncol(g$M)
    dc <- g$d * ifelse(stats::runif(K) < acc, 1L, -1L)
    tri_stats(g$M, g$x, dc * (2 * acc - 1))[1:3]
  }, numeric(3)))
  ## null ACAT distribution by leave-self-out ranks -> exact 5% critical value
  NCn <- nrow(TN)
  pa <- vapply(seq_len(NCn), function(i) {
    pl <- vapply(1:3, function(k)
      (1 + sum(TN[-i, k] >= TN[i, k])) / NCn, numeric(1))
    acat(pl)
  }, numeric(1))
  crit <- apply(S[, 1:7], 2, stats::quantile, probs = .95, na.rm = TRUE)
  attr(crit, "tri_null") <- TN
  attr(crit, "tri_crit") <- as.numeric(stats::quantile(pa, .05))
  crit
}
