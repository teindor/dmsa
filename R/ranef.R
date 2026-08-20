## ---------------------------------------------------------------------------
## Chip as a RANDOM INTERCEPT, per the Project Alpha covariate contract
## (`1 parent T1, chip_T1, random, (1 | chip_T1)`).
##
## DMSA entered chip as a FIXED factor through 1.17.0. On the parent build that
## is 61 levels and 60 df, and it fits 7 singleton chips perfectly. The contract
## asks for a random intercept, which costs a median of ~16 effective df on
## these data and shrinks chip away entirely on the ~20% of probes where the
## chip variance is estimated as zero.
##
## No mixed-model package is required. For a one-way random intercept the GLS
## solution is exactly ordinary least squares on the quasi-demeaned data
## (Balestra-Nerlove):
##
##     theta_c = 1 - 1/sqrt(1 + n_c * gamma),   gamma = sigma^2_u / sigma^2_e
##     v~_i    = v_i - theta_c * mean_c(v)
##
## applied to the response AND every design column. gamma is estimated by REML,
## profiled to a 1-D optimisation. Validated against lme4::lmer(REML = TRUE) on
## the Alpha parent build: coefficients and t-statistics agree to 5e-13, and the
## variance components to 4e-06, including the sigma^2_u = 0 boundary.
##
## Because the transform is a per-probe LINEAR map applied once, before any
## permutation, Freedman-Lane block permutation is unaffected.
## ---------------------------------------------------------------------------

## quasi-demeaning transform for one gamma
.ri_transform <- function(V, grp, nc, gi, gamma) {
  if (!is.finite(gamma) || gamma <= 0) return(V)
  th <- 1 - 1/sqrt(1 + nc[gi] * gamma)
  gm <- apply(V, 2, function(v) tapply(v, grp, mean)[gi])
  V - th * gm
}

## REML for gamma, profiled. Returns gamma = 0 when the chip variance is nil.
.ri_reml <- function(y, X, grp) {
  grp <- as.factor(grp)
  nc <- as.numeric(table(grp)); names(nc) <- levels(grp)
  gi <- as.character(grp); n <- length(y); p <- ncol(X)
  if (n <= p + 1L) return(0)
  nll <- function(lg) {
    gam <- exp(lg)
    Xt <- .ri_transform(X, grp, nc, gi, gam)
    yt <- as.numeric(.ri_transform(cbind(y), grp, nc, gi, gam))
    A <- crossprod(Xt)
    ch <- tryCatch(chol(A), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    bh <- backsolve(ch, backsolve(ch, crossprod(Xt, yt), transpose = TRUE))
    rss <- sum((yt - Xt %*% bh)^2)
    if (!is.finite(rss) || rss <= 0) return(1e10)
    (n - p) * log(rss/(n - p)) + sum(log(1 + nc * gam)) + 2 * sum(log(diag(ch)))
  }
  o <- tryCatch(stats::optimize(nll, c(log(1e-8), log(1e4))),
                error = function(e) NULL)
  if (is.null(o)) return(0)
  ## boundary: if gamma -> 0 is at least as good, there is no chip variance
  if (nll(log(1e-8)) <= o$objective + 1e-8) return(0)
  exp(o$minimum)
}

## Per-probe gamma for a matrix of responses, given one design.
.ri_gammas <- function(Y, X, grp) {
  vapply(seq_len(ncol(Y)), function(j) {
    y <- Y[, j]
    if (!all(is.finite(y))) return(0)
    tryCatch(.ri_reml(y, X, grp), error = function(e) 0)
  }, numeric(1))
}
