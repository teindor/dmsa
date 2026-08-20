# ============================================================================
# RELIABILITY (CENTRALITY) WEIGHTS - the default aligned weighting.
#
# DMSA ships two engines. Both align each probe to its expression consequence;
# they differ only in how much each aligned probe counts toward its unit's
# score:
#
#   "reliability" (DEFAULT) weights probe j by its item-REST correlation - how
#       much it tracks the rest of its own unit - so a probe that reliably
#       reports the unit's shared axis counts more than a noisy one. This is
#       classical reliability weighting (corrected item-total correlation), and
#       in simulation it GAINS power on units whose probe quality is
#       heterogeneous while LOSING nothing when it is uniform (unit weights are
#       already near-optimal there). It is the default because it weakly
#       dominates the flat engine.
#
#   "flat" weights every usable aligned probe equally (the original engine).
#       Retained as an option and as the exact reproducibility reference.
#
# THE CALIBRATION GUARANTEE. The weights are a function of the methylation
# matrix ALONE - never the outcome - so under Freedman-Lane permutation of the
# outcome residuals the weights are frozen and the permutation null is exact.
# Weighting therefore changes power, never type-I error. (Verified in
# tests/testthat/test-weight.R and in the weighted-DMSA simulation.)
#
# THREE GUARDS keep the weights honest:
#   (1) a unit with < 3 usable direction-called probes falls back to flat -
#       the item-rest correlation is too unstable to estimate;
#   (2) a unit whose aligned probes carry no coherent axis - the leading
#       eigenvalue of their correlation matrix is below `w_floor` (default 1.5;
#       independent probes give eigenvalues near 1, so this asks for a shared
#       factor stronger than sampling noise) - falls back to flat; there is
#       nothing to be central to, and reweighting would chase noise;
#   (3) a probe anti-correlated with its unit's axis clamps to weight 0, and if
#       every probe clamps to 0 the unit falls back to flat.
# Surviving weights are renormalised to mean 1 within the unit so the aligned
# score keeps the same scale as the flat engine.
# ============================================================================

## null-coalescing helper (internal; used across the report layer)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Outcome-free reliability weights for aligned probes
#'
#' Item-rest correlation of each aligned probe with the rest of its unit,
#' computed from the (standardised, already winsorised) methylation matrix only.
#' Used internally by \code{\link{dmsa_triangulate}}, \code{\link{dmsa_scores}}
#' and \code{\link{dmsa_change}} when \code{weighting = "reliability"}.
#'
#' @param Z Numeric matrix, samples x probes, standardised methylation (the same
#'   matrix the lenses see; outcome-free).
#' @param units Character/label vector, one per column of \code{Z}.
#' @param mlt Numeric per-probe aligned multiplier (0 for unusable probes).
#' @param w_floor Minimum leading eigenvalue of the aligned probes' correlation
#'   matrix for weighting to engage (independent probes give eigenvalues near 1;
#'   the default 1.5 asks for a shared factor stronger than sampling noise).
#'   Below it the unit is weighted flat.
#' @param weighting \code{"reliability"} (default) or \code{"flat"} (returns all
#'   ones).
#' @return Numeric vector of per-probe weights, length \code{ncol(Z)}, mean 1
#'   within each weighted unit and 1 everywhere a guard fired.
#' @export
dmsa_relweights <- function(Z, units, mlt, w_floor = 1.5,
                            weighting = c("reliability", "flat")) {
  weighting <- match.arg(weighting)
  P <- ncol(Z); relw <- rep(1, P)
  if (weighting == "flat") return(relw)
  for (u in unique(units)) {
    j <- which(units == u & is.finite(mlt) & mlt != 0)
    if (length(j) < 3L) next                              # guard (1): small K
    A <- sweep(Z[, j, drop = FALSE], 2, sign(mlt[j]), "*")  # orient to +axis
    C <- suppressWarnings(stats::cor(A))
    if (any(!is.finite(C))) next
    ev1 <- eigen(C, symmetric = TRUE, only.values = TRUE)$values[1]
    if (!is.finite(ev1) || ev1 < w_floor) next             # guard (2): no axis
    r <- vapply(seq_along(j), function(k) {
      rr <- suppressWarnings(stats::cor(A[, k], rowMeans(A[, -k, drop = FALSE])))
      if (is.finite(rr)) rr else 0 }, numeric(1))
    w <- pmax(0, r)                                        # guard (3): clamp
    if (sum(w > 0) < 2L) next                              # need >=2 for an axis
    relw[j] <- w / mean(w[w > 0])                          # mean 1 within unit
  }
  relw
}
