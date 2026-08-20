# ============================================================================
# LOCKED GENERATOR  -  v1.0, 13 Aug 2026
#
# One data-generating process for all three panels, so nothing differs between
# them except the axis under study.
#
# THE CAUSAL STORY, stated once. The exposure x raises the set's GENE-EXPRESSION
# tone. For a probe whose methylation predicts expression in direction d_j
# (d_j = +1 means more methylation -> more expression), raising expression means
# methylation moves by +d_j. So
#
#     M_j = correlated noise + h_j * d_j * f(x)
#
# and the BIOLOGICAL truth is always "+1: expression tone goes up", whatever the
# mix of d_j. The mean METHYLATION effect, which is all an unaligned method can
# see, is h * (2 * p_plus - 1) - it shrinks to nothing at p_plus = .5 and
# reverses sign across it. That gap is the whole point.
#
#   shape    linear     f(x) = z(x)
#            quadratic  f(x) = -(z(x)^2 - 1), an optimum in mid-range: too low
#                       and too high both cost. Its linear component is ~0.
#            threshold  f(x) = 1(x > q80) - .2, nothing until a level is
#                       crossed, then a step (allostatic load).
#   conc     how many of the G genes carry signal. h is scaled by sqrt(G/conc)
#            so total signal sum(h_j^2) is held constant - dense and sparse
#            differ in WHERE the signal is, not how much there is.
# ============================================================================

N_DEF <- 396L; PER_GENE_DEF <- 5L; RHO_DEF <- .5; NBG_DEF <- 800L

make_gen <- function(N = N_DEF, K = 60L, per_gene = PER_GENE_DEF,
                     rho = RHO_DEF, nbg = NBG_DEF, h = .10, p_plus = .30,
                     shape = c("linear", "quadratic", "threshold"),
                     conc = NULL) {
  shape <- match.arg(shape)
  G  <- as.integer(K %/% per_gene)
  GB <- as.integer(nbg %/% per_gene)
  gid  <- rep(seq_len(G),  each = per_gene)[1:K]
  gidB <- rep(seq_len(GB), each = per_gene)[1:nbg]
  conc <- if (is.null(conc)) G else as.integer(min(conc, G))
  ## matched total signal: concentrating into fewer genes raises each h
  h_gene <- numeric(G)
  h_gene[seq_len(conc)] <- h * sqrt(G / conc)
  hj <- h_gene[gid]

  fx <- switch(shape,
    linear    = function(x) as.numeric(scale(x)),
    quadratic = function(x) { z <- as.numeric(scale(x)); -(z^2 - 1) },
    threshold = function(x) as.numeric(x > stats::quantile(x, .8)) - .2)

  function() {
    x <- stats::rnorm(N)
    d <- sample(c(1L, -1L), K, TRUE, prob = c(p_plus, 1 - p_plus))
    u <- matrix(stats::rnorm(N * G), N, G)[, gid, drop = FALSE]
    M <- sqrt(rho) * u + sqrt(1 - rho) * matrix(stats::rnorm(N * K), N, K) +
      outer(fx(x), hj * d)
    uB <- matrix(stats::rnorm(N * GB), N, GB)[, gidB, drop = FALSE]
    B <- sqrt(rho) * uB + sqrt(1 - rho) * matrix(stats::rnorm(N * nbg), N, nbg)
    list(M = M, B = B, x = x, d = d, gid = gid, gidB = gidB,
         signal_gene = seq_len(conc), signal_probe = which(hj > 0),
         K = K, G = G, conc = conc, h = h, p_plus = p_plus, shape = shape)
  }
}

## the null version of any configuration: same correlation structure, h = 0
make_gen0 <- function(...) {
  a <- list(...); a$h <- 0
  do.call(make_gen, a)
}
