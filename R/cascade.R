# dmsa gated cascade, hybrid inference
# --------------------------------------------------------------------------
# system -> module -> gene -> probe, descending only into rejected parents.
#
# Inference differs by level, because the levels differ in how many units they
# contain:
#   FREQUENTIST at the top. A handful of systems gives an empirical-Bayes
#   procedure nothing to learn a prior from; simulation had EB *losing* to a
#   plain gate there (0.38 vs 0.68 probe sensitivity). The gate combines a DENSE
#   statistic (pooled aligned z - optimal when the whole unit is shifted) with a
#   SPARSE one (max over child units - optimal when one child carries it), which
#   is what fixes dilution: system detection 0.09 -> 0.33 at h = 0.10 for no
#   measurable FDR cost.
#   EMPIRICAL BAYES at the bottom, where hundreds or thousands of units make the
#   mixture estimable and EB was the most sensitive method tested (0.97 vs 0.93).
#   The parent's evidence enters as a prior-odds multiplier for its children -
#   borrowing strength instead of opening or closing a door.

#' Local false discovery rate by central matching
#'
#' Efron-style two-group model on z-values: the null centre and scale are
#' fitted on an iteratively trimmed core (robust to a minority signal slab),
#' \code{pi0} is the observed mass inside the window that holds 50% under
#' the fitted null, and the alternative is a wider normal. Returns one local
#' fdr per unit.
#'
#' Assumptions, stated plainly (calibrated by simulation, 2026-08-29):
#' signal must be a MINORITY of units - power is good to roughly a third
#' signal and collapses conservatively (toward selecting nothing) beyond
#' that; weak signal (|z| ~ 2) is largely invisible to the EB arm (also
#' conservative); and the null must be approximately GAUSSIAN. A null with
#' both tails significantly heavier than the fitted Gaussian - kurtosis,
#' which the two-group model would misread as signal - is detected by a
#' two-radius symmetric-tail-excess guard and the function returns
#' \code{NULL} with a warning, the callers' documented fall-back to the
#' frequentist rule. The guard keys on BOTH tails being in excess, so
#' one-sided signal does not trip it (~3%, same as a pure Gaussian null);
#' symmetric two-sided signal DOES trip it and is handled by the frequentist
#' arm instead - a conservative, never a silent, outcome. Residual risk:
#' with only a few hundred units a heavy-tailed null is detected in about
#' two runs of three; below that scale prefer the frequentist engine.
#'
#' @param z numeric vector of unit-level z statistics (null approximately
#'   standard normal after calibration).
#' @param prior_odds_multiplier Multiplies the estimated \code{pi1}, bounded to
#'   \code{[0.01, 0.9]}. Values above 1 encode evidence inherited from a parent.
#' @param min_n Below this many units the estimate is not attempted and
#'   \code{NULL} is returned, signalling the caller to fall back to a
#'   frequentist rule.
#' @return numeric vector of local fdr values, or NULL if not estimable.
#' @examples
#' set.seed(1)
#' ## 300 null genes and 30 carrying a real aligned effect
#' z <- c(rnorm(300), rnorm(30, mean = 4))
#' lfdr <- dmsa_lfdr(z)
#' round(quantile(lfdr, c(0, .5, 1)), 3)
#' mean(lfdr[301:330]) < mean(lfdr[1:300])
#' ## too few units to fit the two-group model: NULL asks the caller to
#' ## fall back to a frequentist rule rather than guess
#' is.null(dmsa_lfdr(rnorm(10)))
#' @export
dmsa_lfdr <- function(z, prior_odds_multiplier = 1, min_n = 50L) {
  z <- as.numeric(z)
  if (length(z) < min_n || stats::sd(z) == 0) return(NULL)
  ## E5 fix (2026-08-29, PI-approved). The old estimator read the mass
  ## between the sample's OWN quartiles - which is 0.5 by definition of
  ## quartiles - so pi0 was 1.0 whatever the data and the whole EB arm ran
  ## with a hard-coded ~1% signal prior, deaf to the data.
  ##
  ## Efron-style central matching, made robust to a heavy slab: the null's
  ## centre and scale are estimated on an iteratively TRIMMED core (values
  ## within 1.5 fitted SDs of the fitted centre, scale corrected for the
  ## truncation), so a signal cluster cannot drag the fitted null toward
  ## itself the way raw quantiles can. pi0 is then the observed mass inside
  ## the FIXED window delta +/- qnorm(.75)*sigma0 - the interval that holds
  ## exactly 50% under the fitted null - divided by 0.5. Pure noise gives
  ## pi0 ~ 1 (the old behaviour, now earned rather than hard-coded); real
  ## signal pushes mass out of the window and pi0 falls, so selection opens.
  delta <- stats::median(z)
  sigma0 <- max(stats::mad(z), 1e-6)
  .tf <- 1.5
  ## sd of a normal truncated at +/- 1.5 sd understates sigma by this factor
  .corr <- sqrt(1 - 2 * .tf * stats::dnorm(.tf) / (2 * stats::pnorm(.tf) - 1))
  for (.it in 1:3) {
    core <- z[abs(z - delta) <= .tf * sigma0]
    if (length(core) < min_n / 2) break
    delta <- stats::median(core)
    s_new <- stats::sd(core) / .corr
    if (!is.finite(s_new) || s_new <= 0) break
    sigma0 <- s_new
  }
  half <- stats::qnorm(.75) * sigma0
  centre <- mean(z >= delta - half & z <= delta + half)
  pi0 <- min(1, max(.5, centre / .5))
  ## E5 undermining (2026-08-29): the two-group model assumes a GAUSSIAN
  ## null, and a heavy-tailed null (kurtosis, not signal) is read as signal -
  ## under a pure t5 null the selection arm called false units in 80-100% of
  ## simulated runs. Kurtosis has a signature real signal rarely has: BOTH
  ## tails significantly heavier than the fitted Gaussian, in proportion.
  ## One-sided or shouldered signal leaves at most one tail in excess. When
  ## both tails exceed their Poisson 99.5% bound the fit is declared
  ## untrustworthy and NULL is returned - the callers' documented fallback
  ## to the frequentist arm, a safe failure instead of a silent one.
  ## Two radii, OR-ed, each at the .95 Poisson bound per tail: BOTH tails
  ## must be in excess (the symmetry requirement is what protects power -
  ## one-sided signal loads one tail only), and under a true Gaussian null
  ## the joint two-tail false-trip is ~2 (0.05)^2 per radius, well under 1%.
  .trip <- function(r) {
    lam <- length(z) * stats::pnorm(-r)
    cut <- stats::qpois(.95, lam)
    lo <- sum(z < delta - r * sigma0); hi <- sum(z > delta + r * sigma0)
    lo > cut && hi > cut
  }
  if (.trip(2.5) || .trip(3)) {
    .lo <- sum(z < delta - 3 * sigma0); .hi <- sum(z > delta + 3 * sigma0)
    warning("dmsa_lfdr: both tails are heavier than the fitted Gaussian ",
            "null (", .lo, " and ", .hi, " units beyond 3 SD). ",
            "This is the signature of a heavy-tailed null, ",
            "which the two-group model would misread as signal. Falling ",
            "back to the frequentist rule.", call. = FALSE)
    return(NULL)
  }
  pi1 <- min(.9, max(.01, (1 - pi0) * prior_odds_multiplier))
  ## slab width by variance matching
  tau2 <- max((stats::var(z) - sigma0^2) / max(pi1, .01), sigma0^2 * .25)
  f0 <- stats::dnorm(z, delta, sigma0)
  f1 <- stats::dnorm(z, delta, sqrt(sigma0^2 + tau2))
  lf <- (1 - pi1) * f0 / ((1 - pi1) * f0 + pi1 * f1)
  pmin(1, pmax(0, lf))
}

#' Select units by Bayesian FDR
#'
#' Largest set whose MEAN local fdr does not exceed \code{q}.
#' @param lfdr vector of local fdr values
#' @param q target Bayesian FDR
#' @return integer indices of selected units
#' @examples
#' lfdr <- c(rep(0.001, 10), rep(0.5, 10), rep(0.99, 10))
#' sel <- dmsa_bfdr_select(lfdr, q = 0.05)
#' sel
#' mean(lfdr[sel])
#' # nothing is selected when no unit is credible enough to fit the budget
#' dmsa_bfdr_select(rep(0.9, 20), q = 0.05)
#' @export
dmsa_bfdr_select <- function(lfdr, q = 0.05) {
  if (!length(lfdr)) return(integer(0))
  o <- order(lfdr); cm <- cumsum(lfdr[o]) / seq_along(o)
  k <- suppressWarnings(max(which(cm <= q)))
  if (!is.finite(k)) integer(0) else o[seq_len(k)]
}

## BH over the TESTABLE units only. A unit with no usable probe has no p-value;
## counting it in the denominator would make every level needlessly conservative
## in exactly the datasets where coverage is thin.
.bh_sel <- function(p, q) {
  ok <- which(is.finite(p))
  if (!length(ok)) return(integer(0))
  ok[stats::p.adjust(p[ok], "BH") < q]
}

## Calibrated p from a stored null vector, with the +1 correction.
##
## The null vector routinely contains NA: a unit whose probes all carry a zero
## multiplier - no direction call, or a gene the polarity table puts off-axis -
## has no pooled statistic under permutation either. Those entries must be
## DROPPED, not propagated. Without that, a single unusable gene turned every
## p-value at its level into NA and the whole level silently selected nothing,
## which is exactly what happened on the first real-data run.
.cal_p <- function(stat, null, absolute = TRUE) {
  null <- null[is.finite(null)]
  if (!length(null) || !is.finite(stat)) return(NA_real_)
  if (absolute) (1 + sum(abs(null) >= abs(stat))) / (length(null) + 1)
  else (1 + sum(null >= stat)) / (length(null) + 1)
}

#' Gated hierarchical DMSA with level-specific inference
#'
#' @param b,se Probe-level coefficient and standard error, one per probe.
#' @param alignment A \code{dmsa_align()} result covering the same probes in
#'   the same order (supplies the aligned multipliers).
#' @param tree data.frame with one row per probe and one column per level,
#'   ordered outermost first, e.g. \code{data.frame(system=, module=, gene=)}.
#'   Probes are the implicit innermost level.
#' @param nulls Named list of null statistic vectors for calibration, one entry
#'   per level named as in \code{tree}, plus \code{"sparse"} entries named
#'   \code{paste0(level, "_sparse")} where a sparse statistic is used. Generate
#'   these once per design with \code{dmsa_cascade_null()}.
#' @param engine Named character vector giving \code{"freq"} or \code{"eb"} per
#'   level (plus \code{"probe"}).
#'   Default is "freq" everywhere: in a gated tree the surviving unit counts are
#'   too small for an empirical-Bayes mixture to be estimable, and validation showed
#'   EB neutral-to-worse there. EB belongs on an ungated level with many units.
#' @param gate \code{"both"} (default), \code{"dense"} or \code{"sparse"} - the
#'   statistic used at frequentist levels. \code{"both"} takes the smaller of
#'   the two calibrated p-values with a factor-2 correction.
#' @param q Target error rate at each level (BH for \code{"freq"}, Bayesian FDR
#'   for \code{"eb"}).
#' @param method Pooling method passed to \code{dmsa_test()}.
#' @param borrow Upper bound on the prior-odds multiplier passed from a parent
#'   to its children at empirical-Bayes levels. 1 disables borrowing.
#' @return An object of class \code{dmsa_cascade}.
#' @examples
#' set.seed(1)
#' tree <- expand.grid(probe = 1:3, gene = 1:2, system = 1:4)
#' tree <- tree[, c("system", "gene")]      # outermost level first
#' P <- nrow(tree); n <- 60; x <- rnorm(n)
#' M <- matrix(rnorm(n * P), n, P) + outer(x, 0.9 * (tree$system == 1))
#' al <- dmsa_align(data.frame(cpg = paste0("p", 1:P), d = 1, p_plus = 0.95),
#'                  genes = paste0("g", tree$gene), level = "gene")
#' fit <- lm(M ~ x); b <- coef(fit)["x", ]
#' se <- sqrt(colSums(residuals(fit)^2) / (n - 2) / sum((x - mean(x))^2))
#' nulls <- dmsa_cascade_null(b, se, al, tree, exposure = x, M = M, B = 49)
#' cs <- dmsa_cascade(b, se, al, tree, nulls)
#' cs$tables$system
#' @export
dmsa_cascade <- function(b, se, alignment, tree, nulls,
                         engine = NULL, gate = c("both", "dense", "sparse"),
                         q = 0.05, method = c("expected", "fixed"),
                         borrow = 4) {
  gate <- match.arg(gate); method <- match.arg(method)
  tree <- as.data.frame(tree, stringsAsFactors = FALSE)
  al <- as.data.frame(alignment)
  if (nrow(tree) != length(b) || length(b) != length(se) || nrow(al) != length(b))
    stop("b, se, alignment and tree must all have one entry per probe", call. = FALSE)
  levs <- names(tree)
  ## Default: frequentist at every level with the dense+sparse gate. An earlier
  ## exploratory run suggested empirical Bayes should take over at the lower
  ## levels; validating the shipped function showed otherwise. Inside a GATED
  ## tree the surviving unit counts stay small even pooled across parents, and
  ## EB was neutral-to-worse (probe sensitivity 0.16 vs 0.17 and 0.29 vs 0.32,
  ## at equal or higher FDR). EB earns its place on an UNGATED level with
  ## thousands of units - i.e. the flat unannotated-probe arm - not here.
  ## Set engine = c(..., gene = "eb", probe = "eb") to opt in.
  if (is.null(engine))
    engine <- stats::setNames(rep("freq", length(levs) + 1L), c(levs, "probe"))
  ## cumulative unit id at each level, so a child key is unique
  keys <- lapply(seq_along(levs), function(i)
    do.call(paste, c(unname(tree[seq_len(i)]), sep = "|")))
  names(keys) <- levs

  pooled <- function(pos) {
    if (!length(pos)) return(NA_real_)
    dmsa_test(b[pos], se[pos], al[pos, , drop = FALSE], method = method)$z
  }

  out <- list(); parent_evidence <- c(`__root__` = 1)
  alive <- list(`__root__` = seq_len(nrow(tree)))

  for (i in seq_along(levs)) {
    lv <- levs[i]; k <- keys[[lv]]
    rows <- list(); nxt <- list(); ev <- c()
    ## An empirical-Bayes level must see ALL children of ALL surviving parents at
    ## once: inside a single parent there are only a handful of units, far too
    ## few to estimate a mixture. Gather first, fit once, select globally.
    if (engine[[lv]] == "eb") {
      allu <- unlist(lapply(names(alive), function(pk)
        stats::setNames(unique(k[alive[[pk]]]),
                        rep(pk, length(unique(k[alive[[pk]]]))))), use.names = TRUE)
      if (length(allu)) {
        posof <- lapply(seq_along(allu), function(j)
          alive[[names(allu)[j]]][k[alive[[names(allu)[j]]]] == allu[j]])
        zz <- vapply(posof, pooled, numeric(1))
        mult <- vapply(names(allu), function(pk)
          min(borrow, max(1, parent_evidence[[pk]])), numeric(1))
        lf <- dmsa_lfdr(zz, prior_odds_multiplier = stats::median(mult))
        if (is.null(lf)) {
          pv <- vapply(zz, .cal_p, numeric(1), null = nulls[[lv]])
          sel <- .bh_sel(pv, q); stat_p <- pv
          eng_used <- "freq (EB not estimable)"
        } else { sel <- dmsa_bfdr_select(lf, q); stat_p <- lf; eng_used <- "eb" }
        out[[lv]] <- data.frame(
          level = lv, parent = names(allu), unit = unname(allu), z = zz,
          sparse_z = NA_real_, stat = stat_p,
          selected = seq_along(allu) %in% sel,
          n_probes = vapply(posof, length, integer(1)),
          engine = eng_used, stringsAsFactors = FALSE)
        nxt <- stats::setNames(posof[sel], unname(allu[sel]))
        ev <- stats::setNames(
          vapply(sel, function(s) min(borrow, 1 / max(stat_p[s], 1e-3)), numeric(1)),
          unname(allu[sel]))
      } else { out[[lv]] <- NULL }
      alive <- nxt; parent_evidence <- ev
      if (!length(alive)) break
      next
    }
    for (pk in names(alive)) {
      pos <- alive[[pk]]
      units <- unique(k[pos])
      if (!length(units)) next
      zz <- vapply(units, function(u) pooled(pos[k[pos] == u]), numeric(1))
      ## sparse companion: the largest child statistic inside each unit
      sp <- rep(NA_real_, length(units))
      if (gate != "dense" && i < length(levs)) {
        ck <- keys[[levs[i + 1L]]]
        sp <- vapply(units, function(u) {
          p2 <- pos[k[pos] == u]; ch <- unique(ck[p2])
          if (length(ch) < 2) return(NA_real_)
          max(abs(vapply(ch, function(c2) pooled(p2[ck[p2] == c2]), numeric(1))), na.rm = TRUE)
        }, numeric(1))
      }
      if (engine[[lv]] == "freq") {
        pd <- vapply(zz, .cal_p, numeric(1), null = nulls[[lv]])
        pv <- pd
        if (gate != "dense" && any(is.finite(sp)) && !is.null(nulls[[paste0(lv, "_sparse")]])) {
          ps <- vapply(sp, function(s) if (is.na(s)) 1 else
            .cal_p(s, nulls[[paste0(lv, "_sparse")]]), numeric(1))
          pv <- if (gate == "sparse") ps else pmin(1, 2 * pmin(pd, ps))
        }
        sel <- .bh_sel(pv, q)
        stat_p <- pv
      } else {
        lf <- dmsa_lfdr(zz, prior_odds_multiplier =
                          min(borrow, max(1, parent_evidence[[pk]])))
        if (is.null(lf)) {                       # too few units: fall back
          pv <- vapply(zz, .cal_p, numeric(1), null = nulls[[lv]])
          sel <- .bh_sel(pv, q); stat_p <- pv
        } else {
          sel <- dmsa_bfdr_select(lf, q); stat_p <- lf
        }
      }
      rows[[length(rows) + 1L]] <- data.frame(
        level = lv, parent = pk, unit = units, z = zz, sparse_z = sp,
        stat = stat_p, selected = seq_along(units) %in% sel,
        n_probes = vapply(units, function(u) sum(k[pos] == u), integer(1)),
        engine = "freq", stringsAsFactors = FALSE)
      for (s in sel) {
        u <- units[s]
        nxt[[u]] <- pos[k[pos] == u]
        ev[u] <- if (engine[[lv]] == "freq") min(borrow, 1 / max(stat_p[s], 1e-3))
                 else min(borrow, 1 / max(stat_p[s], 1e-3))
      }
    }
    out[[lv]] <- do.call(rbind, rows)
    alive <- nxt; parent_evidence <- ev
    if (!length(alive)) break
  }

  ## probe layer inside surviving innermost units
  probe_sel <- integer(0); probe_tab <- NULL
  if (length(alive)) {
    zp <- b / se
    if (engine[["probe"]] == "eb") {
      ## again: pool every probe under every surviving parent before estimating
      pos <- unlist(alive, use.names = FALSE)
      par <- rep(names(alive), vapply(alive, length, integer(1)))
      mult <- vapply(names(alive), function(pk)
        min(borrow, max(1, parent_evidence[[pk]])), numeric(1))
      lf <- dmsa_lfdr(zp[pos], prior_odds_multiplier = stats::median(mult))
      if (is.null(lf)) {
        pv <- 2 * stats::pnorm(-abs(zp[pos]))
        sel <- .bh_sel(pv, q); stat_p <- pv
        eng <- "freq (EB not estimable)"
      } else { sel <- dmsa_bfdr_select(lf, q); stat_p <- lf; eng <- "eb" }
      probe_sel <- pos[sel]
      probe_tab <- data.frame(level = "probe", parent = par,
        unit = as.character(pos), z = zp[pos], sparse_z = NA_real_,
        stat = stat_p, selected = seq_along(pos) %in% sel, n_probes = 1L,
        engine = eng, stringsAsFactors = FALSE)
    } else {
      rows <- list()
      for (pk in names(alive)) {
        pos <- alive[[pk]]
        pv <- 2 * stats::pnorm(-abs(zp[pos]))
        sel <- .bh_sel(pv, q)
        probe_sel <- c(probe_sel, pos[sel])
        rows[[length(rows) + 1L]] <- data.frame(
          level = "probe", parent = pk, unit = as.character(pos), z = zp[pos],
          sparse_z = NA_real_, stat = pv,
          selected = seq_along(pos) %in% sel, n_probes = 1L,
          engine = "freq", stringsAsFactors = FALSE)
      }
      probe_tab <- do.call(rbind, rows)
    }
  }
  out$probe <- probe_tab

  structure(list(levels = c(levs, "probe"), engine = engine, gate = gate, q = q,
                 method = method, borrow = borrow, tables = out,
                 selected_probes = sort(unique(probe_sel))),
            class = "dmsa_cascade")
}

#' Null statistics for cascade calibration
#'
#' Simulates the pooled aligned statistic under the global null for each unit
#' shape present in the tree, so that \code{dmsa_cascade()} can convert observed
#' statistics into calibrated p-values. Uses the supplied residual covariance
#' structure by resampling the exposure within blocks.
#'
#' @param b,se,alignment,tree as in \code{dmsa_cascade}.
#' @param exposure numeric exposure vector used to fit \code{b}.
#' @param M matrix of probe values (rows = observations), used to refit under
#'   permutation.
#' @param block optional grouping vector defining exchangeable units.
#' @param B number of permutations.
#' @param method pooling method.
#' @return named list suitable for the \code{nulls} argument.
#' @examples
#' set.seed(1)
#' tree <- expand.grid(probe = 1:3, gene = 1:2, system = 1:4)
#' tree <- tree[, c("system", "gene")]      # outermost level first
#' P <- nrow(tree); n <- 60; x <- rnorm(n)
#' M <- matrix(rnorm(n * P), n, P)
#' al <- dmsa_align(data.frame(cpg = paste0("p", 1:P), d = 1, p_plus = 0.95),
#'                  genes = paste0("g", tree$gene), level = "gene")
#' fit <- lm(M ~ x); b <- coef(fit)["x", ]
#' se <- sqrt(colSums(residuals(fit)^2) / (n - 2) / sum((x - mean(x))^2))
#' # B = 49 keeps the example quick; a real calibration wants B >= 999
#' nulls <- dmsa_cascade_null(b, se, al, tree, exposure = x, M = M, B = 49)
#' vapply(nulls, length, integer(1))
#' @export
dmsa_cascade_null <- function(b, se, alignment, tree, exposure, M,
                              block = NULL, B = 999L,
                              method = c("expected", "fixed")) {
  method <- match.arg(method)
  tree <- as.data.frame(tree, stringsAsFactors = FALSE)
  al <- as.data.frame(alignment)
  levs <- names(tree)
  keys <- lapply(seq_along(levs), function(i)
    do.call(paste, c(unname(tree[seq_len(i)]), sep = "|")))
  names(keys) <- levs
  n <- length(exposure)
  idxs <- if (is.null(block)) lapply(seq_len(B), function(i) sample.int(n))
          else dmsa_block_index(block, B)
  Y <- scale(M)
  res <- stats::setNames(vector("list", 2 * length(levs)),
                         c(levs, paste0(levs, "_sparse")))
  for (bi in seq_len(B)) {
    x <- exposure[idxs[[bi]]]
    X <- cbind(1, x); XtXi <- solve(crossprod(X))
    bh <- XtXi %*% (t(X) %*% Y); r <- Y - X %*% bh
    s2 <- colSums(r^2) / (n - 2)
    bb <- bh[2, ]; ss <- sqrt(s2 * XtXi[2, 2])
    for (i in seq_along(levs)) {
      lv <- levs[i]; k <- keys[[lv]]
      us <- unique(k)
      zz <- vapply(us, function(u) {
        pos <- which(k == u)
        dmsa_test(bb[pos], ss[pos], al[pos, , drop = FALSE], method = method)$z
      }, numeric(1))
      res[[lv]] <- c(res[[lv]], zz)
      if (i < length(levs)) {
        ck <- keys[[levs[i + 1L]]]
        sp <- vapply(us, function(u) {
          p2 <- which(k == u); ch <- unique(ck[p2])
          if (length(ch) < 2) return(NA_real_)
          max(abs(vapply(ch, function(c2) {
            q2 <- p2[ck[p2] == c2]
            dmsa_test(bb[q2], ss[q2], al[q2, , drop = FALSE], method = method)$z
          }, numeric(1))), na.rm = TRUE)
        }, numeric(1))
        res[[paste0(lv, "_sparse")]] <- c(res[[paste0(lv, "_sparse")]], sp[!is.na(sp)])
      }
    }
  }
  res[!vapply(res, is.null, logical(1))]
}

#' @export
print.dmsa_cascade <- function(x, ...) {
  cat("DMSA gated cascade (", x$method, "-sign, gate = ", x$gate, ", q = ",
      x$q, ")\n", sep = "")
  cat("  engine:", paste(sprintf("%s=%s", names(x$engine), x$engine),
                         collapse = "  "), "\n")
  for (lv in x$levels) {
    tb <- x$tables[[lv]]
    if (is.null(tb)) next
    cat(sprintf("  %-7s tested %4d   selected %3d\n", lv, nrow(tb), sum(tb$selected)))
  }
  cat("  probes selected:", length(x$selected_probes), "\n")
  invisible(x)
}
