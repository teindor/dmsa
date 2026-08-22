# ============================================================================
# LONGITUDINAL mDMSA - the 2-wave identity.
#
# With two waves, time x S x E is ALGEBRAICALLY the change-score model: per
# person Delta-S = S_after - S_before on the SAME alignment map, and the
# tested term is Delta-S x E (or Delta-S x E x mod2 for the second moderator).
# Time is not a variable to estimate - Delta IS time. Validated on the Oct-7
# before/after cohort, 14 Aug 2026 (n = 172; see DMSA_CHANGE.md).
#
# Change-scale cell-composition covariates (Delta-Epi, Delta-Fib, Delta-SVs)
# are MANDATORY in real analyses - trauma shifts the leukocyte mix, and that
# is the classic artifact on exactly this design. The baseline score S_before
# enters as a covariate by default (regression-to-the-mean guard).
# ============================================================================

#' Two-wave change-score DMSA: time x methylation x factor
#'
#' Tests Delta-S x exposure (x mod2) per unit through the composite machinery
#' (\code{\link{dmsa_model}}, Freedman-Lane, family blocks), with maxT across
#' the units of the level-local family on one shared permutation stream.
#'
#' @param M_before,M_after numeric matrices, same probes in the same order
#'   (rows = the same people in the same order).
#' @param data data.frame with outcome, exposure, moderators and covariates.
#' @param outcome outcome column (role = "predictor").
#' @param exposure the moderating factor E; with two waves this is the `mod`
#'   of the time x S x E three-way (Delta carries the time axis).
#' @param mod2 optional second moderator: tests Delta-S x E x mod2 (moderators
#'   mean-centered, all lower-order terms in the model, highest-order tested).
#' @param units unit label per probe column (genes, modules, or one system).
#' @param alignment alignment as in \code{\link{dmsa_triangulate}}.
#' @param covariates character vector of covariate columns. Supply the
#'   CHANGE-SCALE cell composition here; a message reminds you if empty.
#' @param baseline include the unit's before-wave score as a covariate.
#' @param role "predictor": outcome ~ dS * E (* mod2) + S0 + covariates.
#'   "outcome": dS ~ E (* mod2) + covariates (exposure predicts the change).
#' @param block exchangeable blocks (e.g. couple id).
#' @param B,correction,seed as elsewhere; correction applies across the units
#' @param weighting Character. Probe weighting engine within a unit: \code{"combined"} (default) fuses the flat and reliability statistics on one shared permutation stream, \code{"flat"} weights every usable aligned probe equally, \code{"reliability"} weights each probe by its item-rest correlation with the rest of its unit. Weights are computed from methylation alone, so the permutation null is unaffected.
#' @param w_floor Numeric. Lower bound applied to reliability weights before normalisation, so a single poorly behaved probe cannot be driven to zero influence. Ignored when \code{weighting = "flat"}.
#'   of the family (maxT default).
#' @return data.frame of class "dmsa_change": unit, n, n_probes, b, t, p,
#'   p_adj; attributes term, role, formula.
#' @examples
#' set.seed(5)
#' n <- 80; g <- rep(c("g1", "g2"), each = 4); P <- length(g)
#' d <- rep(c(1, -1), length.out = P); E <- rnorm(n); y <- rnorm(n)
#' M0 <- matrix(rnorm(n * P), n, P)
#' M1 <- 0.6 * M0 + matrix(rnorm(n * P), n, P)
#' # plant a change in g1 that tracks y x E through the alignment
#' M1[, g == "g1"] <- M1[, g == "g1"] + outer(0.5 * y * E, d[g == "g1"])
#' al <- dmsa_align(data.frame(cpg = paste0("cg", 1:P), d = d,
#'                             p_plus = ifelse(d > 0, 0.9, 0.1)), genes = g)
#' dat <- data.frame(y = y, E = E, dEpi = rnorm(n))
#' dmsa_change(M0, M1, dat, outcome = "y", exposure = "E", units = g,
#'             alignment = al, covariates = "dEpi", B = 99, seed = 1)
#' @export
dmsa_change <- function(M_before, M_after, data, outcome = NULL, exposure,
                        mod2 = NULL, units, alignment, covariates = NULL,
                        baseline = TRUE, role = c("predictor", "outcome"),
                        block = NULL, B = 1999,
                        correction = c("maxT", "minP"),
                        weighting = c("combined", "reliability", "flat"), w_floor = 1.5,
                        seed = 1) {
  role <- match.arg(role); correction <- match.arg(correction)
  weighting <- match.arg(weighting)
  M_before <- as.matrix(M_before); M_after <- as.matrix(M_after)
  stopifnot(ncol(M_before) == ncol(M_after),
            nrow(M_before) == nrow(M_after),
            length(units) == ncol(M_before))
  data <- as.data.frame(data)
  stopifnot(nrow(data) == nrow(M_before))
  block <- .dmsa_rows(block, data, "block"); .dmsa_check_block(block)
  if (!exposure %in% names(data))
    stop("exposure column '", exposure, "' not in data", call. = FALSE)
  if (role == "predictor" && (is.null(outcome) || !outcome %in% names(data)))
    stop("role = 'predictor' needs an `outcome` column", call. = FALSE)
  if (!is.null(mod2) && !mod2 %in% names(data))
    stop("mod2 column '", mod2, "' not in data", call. = FALSE)
  if (is.null(covariates))
    message("no covariates supplied - in real analyses the CHANGE-SCALE cell ",
            "composition (Delta-Epi, Delta-Fib, Delta-control-SVs) is mandatory")
  al <- as.data.frame(alignment)

  d0 <- data
  d0$.E <- as.numeric(scale(d0[[exposure]]))
  if (!is.null(mod2)) d0$.M2 <- as.numeric(scale(d0[[mod2]]))

  engs <- if (weighting == "combined") c("flat", "reliability") else weighting
  cauchy <- function(p) tan((0.5 - pmin(pmax(p, 1e-12), 1 - 1e-12)) * pi)
  ug <- unique(units)
  res <- list(); fobs <- list(); fnul <- list()
  for (u in ug) {
    j <- which(units == u)
    aj <- al[j, , drop = FALSE]
    usable <- if (!is.null(aj$usable)) sum(aj$usable & aj$s != 0, na.rm = TRUE)
      else sum(is.finite(aj$s) & aj$s != 0)
    if (usable < 2) next
    mj <- 2 * as.numeric(aj$p_s_plus) - 1
    us <- if (!is.null(aj$usable)) as.logical(aj$usable) else rep(TRUE, length(mj))
    mj[!is.finite(mj) | !us] <- 0
    ## reliability weights on both waves from the stacked (outcome-free) matrix
    Zc <- scale(rbind(as.matrix(M_before[, j, drop = FALSE]),
                      as.matrix(M_after[, j, drop = FALSE])))
    zk <- apply(Zc, 2, function(z) all(is.finite(z)))
    rw_rel <- rep(1, length(mj))
    if (any(zk)) rw_rel[zk] <- dmsa_relweights(Zc[, zk, drop = FALSE],
                            rep(1L, sum(zk)), mj[zk], w_floor = w_floor,
                            weighting = "reliability")
    covs <- c(if (baseline) ".S0", covariates)
    if (role == "predictor") {
      lhs <- outcome
      inter <- if (is.null(mod2)) ".dS * .E" else ".dS * .E * .M2"
      term <- if (is.null(mod2)) ".dS:.E" else ".dS:.E:.M2"
    } else {
      lhs <- ".dS"
      inter <- if (is.null(mod2)) ".E" else ".E * .M2"
      term <- if (is.null(mod2)) ".E" else ".E:.M2"
    }
    ff <- stats::as.formula(paste(lhs, "~", paste(c(inter, covs), collapse = " + ")))
    ## run each engine on the SAME permutation stream (same seed) -> shared null
    fits <- list()
    for (e in engs) {
      meff <- mj * (if (e == "flat") rep(1, length(mj)) else rw_rel)
      S0 <- dmsa_scores(M_before[, j, drop = FALSE], meff, weighting = "flat")$aligned
      S1 <- dmsa_scores(M_after[, j, drop = FALSE], meff, weighting = "flat")$aligned
      d2 <- d0; d2$.dS <- as.numeric(scale(S1 - S0)); d2$.S0 <- as.numeric(scale(S0))
      ft <- tryCatch(dmsa_model(ff, d2, term, block = block, B = B, seed = seed,
                                nulls = TRUE), error = function(e) NULL)
      if (is.null(ft)) { fits <- NULL; break }
      fits[[e]] <- ft
    }
    if (is.null(fits) || !length(fits)) next
    ref <- fits[[1]]                                   # descriptive b/t (flat if combined)
    if (length(fits) == 1L) {
      so <- abs(ref$t); sn <- abs(ref$null_t)
    } else {                                           # Cauchy-fuse the engines
      ## pooled rank-p (obs + nulls together) so obs and null are exchangeable
      rp <- lapply(fits, function(f) {
        pool <- c(abs(f$t), abs(f$null_t)); m <- length(pool)
        r <- (m - data.table::frank(pool, ties.method = "min") + 1) / m
        list(po = r[1], pn = r[-1]) })
      so <- Reduce(`+`, lapply(rp, function(r) cauchy(r$po)))
      sn <- Reduce(`+`, lapply(rp, function(r) cauchy(r$pn)))
    }
    res[[u]] <- data.frame(unit = u, n = ref$n, n_probes = usable,
      b = ref$b, t = ref$t,
      p = (1 + sum(sn >= so)) / (length(sn) + 1), stringsAsFactors = FALSE)
    fobs[[u]] <- so; fnul[[u]] <- sn
  }
  if (!length(res))
    stop("no unit had >= 2 usable direction-called probes", call. = FALSE)
  out <- do.call(rbind, res); rownames(out) <- NULL
  SN <- do.call(cbind, fnul); so_v <- unlist(fobs)
  if (!is.null(SN) && ncol(SN) > 1) {
    if (correction == "minP" && length(engs) == 1L) {
      PN <- apply(SN, 2, function(v) rank(-v, ties.method = "max") / length(v))
      mn <- apply(PN, 1, min, na.rm = TRUE)
      out$p_adj <- vapply(out$p, function(q)
        (1 + sum(mn <= q, na.rm = TRUE)) / (nrow(PN) + 1), numeric(1))
    } else {                                           # maxT on the (fused) statistic
      mx <- apply(SN, 1, max, na.rm = TRUE)
      out$p_adj <- vapply(so_v, function(o)
        (1 + sum(mx >= o, na.rm = TRUE)) / (length(mx) + 1), numeric(1))
    }
    o <- order(out$p); out$p_adj[o] <- cummax(out$p_adj[o])
  } else out$p_adj <- out$p
  attr(out, "term") <- if (role == "predictor") {
    if (is.null(mod2)) "dS x E (time x S x E, 2-wave identity)"
    else "dS x E x mod2"
  } else paste("exposure ->", "dS")
  attr(out, "role") <- role
  attr(out, "correction") <- correction
  class(out) <- c("dmsa_change", "data.frame")
  out[order(out$p), ]
}

#' @export
print.dmsa_change <- function(x, ...) {
  cat("DMSA change-score test -", attr(x, "term"), "\n")
  cat("  family of", nrow(x), "unit(s),", attr(x, "correction"),
      "across the family\n")
  for (i in seq_len(nrow(x)))
    cat(sprintf("  %-28s n %4d  probes %3d  b %+.3f  t %+.2f  p %.4f  adj %.4f%s\n",
                x$unit[i], x$n[i], x$n_probes[i], x$b[i], x$t[i], x$p[i],
                x$p_adj[i], if (x$p_adj[i] < .05) "  *" else ""))
  invisible(x)
}
