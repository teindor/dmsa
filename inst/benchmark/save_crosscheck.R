# Saves a stratified sample of benchmark datasets, together with the
# re-implemented engines' own p-values on them, so the Mac can re-run the same
# data through the real packages. Kept small enough to email.
source("/home/claude/bench/engines.R")
source("/home/claude/bench/generator.R")
set.seed(4242)

CFG <- list(list(h = 0, p_plus = .30), list(h = .03, p_plus = .30),
            list(h = .03, p_plus = .41), list(h = .05, p_plus = .30),
            list(h = .05, p_plus = .74), list(h = .08, p_plus = .50))
PER <- 8L; K <- 60L; ACC <- .85

## p-values, not the thresholded hits, so the cross-check can compare
## continuously as well as at .05
gt_null <- NULL
cat("building a null reference for the p-value versions ...\n")
g0 <- make_gen0(K = K)
NULLS <- t(vapply(seq_len(600), function(i) {
  g <- g0(); f <- allfit(cbind(g$M, g$B), g$x)
  c(gt = gt_stat(g$M, g$x),
    es = abs(gsea_es(f$t, c(rep(TRUE, K), rep(FALSE, ncol(g$B))))))
}, numeric(2)))

sets <- list()
for (cf in CFG) {
  gen <- make_gen(K = K, h = cf$h, p_plus = cf$p_plus, shape = "linear")
  for (j in seq_len(PER)) {
    g <- gen()
    Y <- cbind(g$M, g$B); f <- allfit(Y, g$x)
    in_set <- c(rep(TRUE, K), rep(FALSE, ncol(g$B)))
    gt <- gt_stat(g$M, g$x)
    es <- abs(gsea_es(f$t, in_set))
    sets[[length(sets) + 1L]] <- list(
      M = round(g$M, 4), x = round(g$x, 5), t = round(f$t, 5), p = f$p,
      in_set = in_set,
      gene_all = c(g$gid, max(g$gid) + g$gidB),
      h = cf$h, p_plus = cf$p_plus,
      mine = c(globaltest_p = (1 + sum(NULLS[, "gt"] >= gt)) / (nrow(NULLS) + 1),
               gsea_p = (1 + sum(NULLS[, "es"] >= es)) / (nrow(NULLS) + 1),
               rra_p = rra_ora(f$p[1:K], g$gid, f$p[-(1:K)], g$gidB)))
  }
}
out <- list(sets = sets,
            note = paste0("n=", N_DEF, ", K=", K, " set probes in genes of 5, ",
                          "800 background probes, rho=.5 within gene, ",
                          "direction accuracy ", ACC, ". 6 configurations x ",
                          PER, " draws."))
saveRDS(out, "/home/claude/outputs/bench_crosscheck_data.rds", compress = "xz")
cat("saved", length(sets), "datasets ->",
    round(file.size("/home/claude/outputs/bench_crosscheck_data.rds") / 1e6, 1),
    "MB\n")
