# ============================================================================
# CANDIDATE ERROR-CONTROL PROCEDURES, one implementation each, so they can be
# compared on the same permutation output instead of argued about.
#
# Every one of these consumes the SAME two objects:
#   z_obs   observed unit-level statistic, one per unit (gene, or probe)
#   Z_null  B x n matrix of the same statistic under the B restricted
#           permutations already computed for the DMSA test
# so switching procedure costs nothing at run time. That is the point: the
# permutations are the expensive part and we already pay for them.
# ============================================================================

## ---------------------------------------------------------------------------
## 1. WESTFALL-YOUNG SINGLE-STEP maxT.  FWER over the n units.
##    p_adj(j) = P( max_k |Z_k| >= |z_j| ) under the permutation null.
## ---------------------------------------------------------------------------
wy_single <- function(z_obs, Z_null) {
  mx <- apply(abs(Z_null), 1, max, na.rm = TRUE)
  mx <- mx[is.finite(mx)]
  vapply(abs(z_obs), function(a)
    if (!is.finite(a)) NA_real_ else (1 + sum(mx >= a)) / (length(mx) + 1),
    numeric(1))
}

## ---------------------------------------------------------------------------
## 2. WESTFALL-YOUNG / ROMANO-WOLF STEP-DOWN maxT.  FWER, uniformly at least as
##    powerful as single-step: after the most extreme unit is rejected it is
##    removed from the maximum, so later units face a lower bar. Monotonicity is
##    enforced so adjusted p-values never decrease down the ordering.
## ---------------------------------------------------------------------------
wy_stepdown <- function(z_obs, Z_null) {
  a <- abs(z_obs); n <- length(a)
  ok <- which(is.finite(a))
  padj <- rep(NA_real_, n)
  if (!length(ok)) return(padj)
  o <- ok[order(a[ok], decreasing = TRUE)]
  A <- abs(Z_null[, o, drop = FALSE])
  A[!is.finite(A)] <- -Inf
  ## running maximum from the tail: max over the units still in play
  M <- t(apply(A, 1, function(r) rev(cummax(rev(r)))))
  B <- nrow(M)
  p <- vapply(seq_along(o), function(i)
    (1 + sum(M[, i] >= a[o[i]], na.rm = TRUE)) / (B + 1), numeric(1))
  padj[o] <- cummax(p)                       # enforce monotonicity
  padj
}

## ---------------------------------------------------------------------------
## 3. BH on the marginal permutation p-values. The baseline that threw AVP away:
##    it never sees the correlation the permutations encode.
## ---------------------------------------------------------------------------
bh_marginal <- function(z_obs, Z_null) {
  B <- nrow(Z_null)
  p <- vapply(seq_along(z_obs), function(j) {
    nl <- abs(Z_null[, j]); nl <- nl[is.finite(nl)]
    if (!is.finite(z_obs[j]) || !length(nl)) return(NA_real_)
    (1 + sum(nl >= abs(z_obs[j]))) / (length(nl) + 1)
  }, numeric(1))
  q <- rep(NA_real_, length(p)); k <- which(is.finite(p))
  if (length(k)) q[k] <- stats::p.adjust(p[k], "BH")
  q
}

## ---------------------------------------------------------------------------
## 4. PERMUTATION CLOSED TESTING FOR SUM TESTS  (Vesely, Finos & Goeman 2023).
##    Simultaneous lower confidence bound on the NUMBER of true discoveries in
##    any set S, valid for every S at once, so S may be chosen post hoc.
##
##    Centred statistics  c_j^pi = t_j^obs - t_j^pi  (the observed row is 0).
##    Single-step shortcut: phi(z) = 1 iff every set overlapping S in at least z
##    units is rejected. d(S) = |S| - max{ z : phi(z) = 0 }.
##
##    Input here is a PROBE-level statistic (t) and its permutation matrix, not
##    a z; any sum-type statistic is admissible.
## ---------------------------------------------------------------------------
.tdp_setup <- function(t_obs, T_null, alpha = .05) {
  ## The identity transformation belongs in the group and contributes a row of
  ## zeros; leaving it out shifts the quantile and breaks the guarantee.
  C <- rbind(0, sweep(-T_null, 2, t_obs, "+"))   # (B+1) x m centred statistics
  C[!is.finite(C)] <- 0
  list(C = C, omega = floor(alpha * nrow(C)) + 1L, m = ncol(C))
}

.phi <- function(z, idx, S) {
  C <- S$C; m <- S$m; om <- S$omega
  ins <- idx; out <- setdiff(seq_len(m), idx)
  if (z > length(ins)) return(1L)
  sIn  <- if (length(ins)) t(apply(C[, ins, drop = FALSE], 1, sort)) else NULL
  sOut <- if (length(out)) t(apply(C[, out, drop = FALSE], 1, sort)) else NULL
  cumIn  <- if (z > 0) rowSums(sIn[, seq_len(z), drop = FALSE]) else rep(0, nrow(C))
  cOut <- if (!is.null(sOut)) cbind(0, t(apply(sOut, 1, cumsum))) else
    matrix(0, nrow(C), 1)
  for (k in 0:(ncol(cOut) - 1L)) {
    b <- cumIn + cOut[, k + 1L]
    if (sort(b, partial = om)[om] <= 0) return(0L)
  }
  1L
}

#' Lower confidence bound on the number of true discoveries in set S
tdp_bound <- function(t_obs, T_null, idx, alpha = .05) {
  S <- .tdp_setup(t_obs, T_null, alpha)
  n <- length(idx)
  ## d(S) = n - q, where q = max{ z in 0..n : phi(z) = 0 }. phi is monotone
  ## non-decreasing in z. The whole range MUST be searched: if phi(n) = 0 the
  ## answer is 0, and an earlier version that never evaluated z = n could not
  ## return 0 at all - it bottomed out at 1 and so "found" a discovery under
  ## every global null.
  if (.phi(n, idx, S) == 0L) return(0L)
  if (.phi(0L, idx, S) == 1L) return(n)
  lo <- 0L; hi <- n                     # phi(lo)=0, phi(hi)=1
  while (hi - lo > 1L) {
    mid <- (lo + hi) %/% 2L
    if (.phi(mid, idx, S) == 1L) hi <- mid else lo <- mid
  }
  n - lo
}

## ---------------------------------------------------------------------------
## 5. ADAPTIVE TRUNCATION (ARTP, Yu et al. 2009) - the dense/sparse hedge with
##    NO fixed penalty. Instead of min(1, 2*min(p_dense, p_sparse)) we take the
##    minimum over a family of truncation points and calibrate that minimum
##    against the same permutations, so the correlation between arms is priced
##    by the data rather than by Bonferroni.
## ---------------------------------------------------------------------------
artp <- function(z_obs, Z_null, Ks = c(1, 2, 3, 5, 10, Inf)) {
  a <- abs(z_obs); A <- abs(Z_null)
  a[!is.finite(a)] <- 0; A[!is.finite(A)] <- 0
  n <- length(a); B <- nrow(A)
  Ks <- unique(pmin(Ks, n)); Ks <- Ks[Ks >= 1]
  srt <- function(v) sort(v, decreasing = TRUE)
  stat <- function(v) vapply(Ks, function(K) sum(srt(v)[seq_len(K)]), numeric(1))
  s_obs <- stat(a)
  S_null <- t(apply(A, 1, stat))                       # B x length(Ks)
  ## per-K permutation p, for the observed and for every permutation
  pK_obs <- vapply(seq_along(Ks), function(k)
    (1 + sum(S_null[, k] >= s_obs[k])) / (B + 1), numeric(1))
  pK_null <- vapply(seq_along(Ks), function(k) {
    r <- rank(-S_null[, k], ties.method = "min"); r / B }, numeric(B))
  minp_obs  <- min(pK_obs)
  minp_null <- apply(pK_null, 1, min)
  list(p = (1 + sum(minp_null <= minp_obs)) / (B + 1),
       best_K = Ks[which.min(pK_obs)], per_K = stats::setNames(pK_obs, Ks))
}

## ---------------------------------------------------------------------------
## 6. The gate as currently shipped, for reference.
## ---------------------------------------------------------------------------
gate_both <- function(p_dense, p_sparse) min(1, 2 * min(p_dense, p_sparse))
