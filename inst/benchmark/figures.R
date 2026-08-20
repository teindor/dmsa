suppressPackageStartupMessages({library(data.table); library(viridisLite)})
A <- rbindlist(lapply(list.files("/home/claude/outputs", "^bench_panelA_.*csv$",
                                 full.names = TRUE), fread)); setorder(A, h, p_plus)
B <- fread("/home/claude/outputs/bench_panelB.csv"); setorder(B, -conc, h)
C <- fread("/home/claude/outputs/bench_panelC.csv")

VC  <- viridis(9, end = .95)
COL <- c(dmsa_expected = VC[1], ndmsa = VC[2], dmsa_signblind = "grey55",
         globaltest = VC[4], roast_mixed = VC[5], ora = VC[6],
         camera = VC[7], fry = VC[8], roast_dir = VC[8], gsea_mcsea = VC[9],
         methylRRA = "grey75", ewas_bh = "grey35")
LTY <- c(dmsa_expected = 1, ndmsa = 1, dmsa_signblind = 3, globaltest = 2,
         roast_mixed = 2, ora = 2, camera = 1, fry = 5, roast_dir = 4,
         gsea_mcsea = 1, methylRRA = 3, ewas_bh = 3)

## ===========================================================================
## FIGURE 7 - Panel A: the direction ratio
## ===========================================================================
png("/home/claude/outputs/fig7_bench_ratio.png", 2300, 900, res = 190)
par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3.4, .8), mgp = c(2.5, .7, 0))

sh <- c("dmsa_expected", "ndmsa", "globaltest", "roast_mixed", "ora",
        "camera", "fry", "gsea_mcsea", "ewas_bh")
nm <- c("DMSA", "nDMSA", "globaltest", "roast (mixed)", "ORA", "camera",
        "fry", "GSEA/mCSEA", "EWAS+BH")
d <- A[h == .03]
plot(NA, xlim = c(.18, .76), ylim = c(0, 1.28), xlab = expression(P(d[j] == +1)),
     ylab = "power", main = "A. Power at h = .03")
abline(v = .5, col = "grey85", lwd = 6); abline(h = .05, lty = 3, col = "grey50")
for (i in seq_along(sh)) lines(d$p_plus, d[[paste0("pow_", sh[i])]],
                               col = COL[sh[i]], lty = LTY[sh[i]], lwd = 2.4,
                               type = "b", pch = 16, cex = .6)
text(.5, .02, "50/50 tissue", cex = .62, col = "grey30", srt = 90, adj = 0)
legend("top", bty = "n", cex = .58, lwd = 2.2, col = COL[sh], lty = LTY[sh], ncol = 3,
       legend = nm, seg.len = 2.4)

## B: biological direction
plot(NA, xlim = c(.18, .76), ylim = c(0, 1.02), xlab = expression(P(d[j] == +1)),
     ylab = "P(direction reported is biologically correct)",
     main = "B. Which way does the exposure push?")
abline(v = .5, col = "grey85", lwd = 6); abline(h = .5, lty = 3, col = "grey60")
lines(d$p_plus, d$bio_dmsa_expected, col = COL["dmsa_expected"], lwd = 3.2,
      type = "b", pch = 16)
for (e in c("camera", "fry", "roast_dir", "gsea_mcsea"))
  lines(d$p_plus, d[[paste0("bio_", e)]], col = COL[e], lty = LTY[e], lwd = 2.2,
        type = "b", pch = 1, cex = .6)
text(.24, .93, "DMSA", col = COL["dmsa_expected"], cex = .78, font = 2, adj = 0)
text(.24, .08, "camera / fry / roast / GSEA", col = COL["camera"], cex = .7,
     adj = 0)
text(.62, .34, "same biology,\nopposite conclusion", cex = .64, col = "grey20",
     font = 3)

## C: DMSA advantage
plot(NA, xlim = c(.18, .76), ylim = c(0, 12), xlab = expression(P(d[j] == +1)),
     ylab = "DMSA power / competitor power",
     main = "C. How much the ratio costs")
abline(v = .5, col = "grey85", lwd = 6); abline(h = 1, lty = 3, col = "grey50")
for (e in c("camera", "fry", "roast_dir", "gsea_mcsea"))
  lines(d$p_plus, d$pow_dmsa_expected / d[[paste0("pow_", e)]], col = COL[e],
        lty = LTY[e], lwd = 2.2, type = "b", pch = 1, cex = .6)
for (e in c("globaltest", "roast_mixed", "ora", "dmsa_signblind"))
  lines(d$p_plus, d$pow_dmsa_expected / d[[paste0("pow_", e)]], col = COL[e],
        lty = LTY[e], lwd = 2.2, type = "b", pch = 16, cex = .5)
legend("topright", bty = "n", cex = .62,
       legend = c("directional competitors (open)",
                  "directionless competitors (filled)"), pch = c(1, 16))
dev.off()

## ===========================================================================
## FIGURE 8 - Panel B: where the signal sits
## ===========================================================================
png("/home/claude/outputs/fig8_bench_cascade.png", 2300, 900, res = 190)
par(mfrow = c(1, 3), mar = c(5.2, 4.3, 3.4, .8), mgp = c(2.5, .7, 0))

d5 <- B[h == .05]; d5 <- d5[order(-conc)]
m <- rbind(d5$pow_dmsa_expected, d5$casc_sys, d5$pow_globaltest,
           d5$pow_camera, d5$pow_ewas_bh)
bp <- barplot(m, beside = TRUE, ylim = c(0, 1.30), border = NA,
              col = VC[c(1, 2, 4, 7, 9)], names.arg = rep("", 3),
              ylab = "P(set detected)", main = "A. Detection, matched signal")
box(bty = "l")
text(colMeans(bp), -.12, c("dense\n12 genes", "4 genes", "sparse\n1 gene"),
     xpd = NA, cex = .78)
legend("top", bty = "n", cex = .6, ncol = 3, fill = VC[c(1, 2, 4, 7, 9)],
       border = NA,
       legend = c("DMSA (pooled)", "cascade-DMSA", "globaltest", "camera",
                  "EWAS+BH"))
mtext("the ordering reverses", 3, .1, cex = .62, col = "grey30")

m2 <- rbind(d5$casc_gene_sens, d5$bhgw_gene_sens, d5$bhset_gene_sens)
bp2 <- barplot(m2, beside = TRUE, ylim = c(0, 1.05), border = NA,
               col = VC[c(2, 5, 8)], names.arg = rep("", 3),
               ylab = "gene-level sensitivity",
               main = "B. Localization: which gene?")
box(bty = "l")
text(colMeans(bp2), -.10, c("dense\n12 genes", "4 genes", "sparse\n1 gene"),
     xpd = NA, cex = .78)
legend("topleft", bty = "n", cex = .66, fill = VC[c(2, 5, 8)], border = NA,
       legend = c("cascade-DMSA", "BH genome-wide", "BH within set"))

m3 <- rbind(d5$casc_gene_fdr, d5$bhgw_gene_fdr, d5$bhset_gene_fdr)
bp3 <- barplot(m3, beside = TRUE, ylim = c(0, .12), border = NA,
               col = VC[c(2, 5, 8)], names.arg = rep("", 3),
               ylab = "gene-level FDR", main = "C. And at what error rate")
abline(h = .05, lty = 2, col = "grey25"); box(bty = "l")
text(colMeans(bp3), -.012, c("dense\n12 genes", "4 genes", "sparse\n1 gene"),
     xpd = NA, cex = .78)
text(bp3[3, 3], d5$bhset_gene_fdr[3] + .006, "BH within set\nexceeds q = .05",
     cex = .6, col = "grey15", font = 3)
dev.off()

## ===========================================================================
## FIGURE 9 - Panel C: functional form
## ===========================================================================
png("/home/claude/outputs/fig9_bench_shape.png", 2300, 900, res = 190)
par(mfrow = c(1, 3), mar = c(4.3, 4.3, 3.4, .8), mgp = c(2.5, .7, 0))
shp <- c("linear", "quadratic", "threshold")
ttl <- c("A. Linear - nDMSA pays its price",
         "B. Optimum in mid-range - linear tests are blind",
         "C. Threshold - a step, not a slope")
eg <- c("ndmsa", "dmsa_expected", "globaltest", "ora", "camera", "gsea_mcsea",
        "ewas_bh")
for (k in seq_along(shp)) {
  d <- C[shape == shp[k]][order(h)]
  plot(NA, xlim = range(d$h), ylim = c(0, 1.02), xlab = "effect size h",
       ylab = "power", main = ttl[k], cex.main = .96)
  abline(h = .05, lty = 3, col = "grey50")
  for (e in eg) lines(d$h, d[[paste0("pow_", e)]], col = COL[e], lty = LTY[e],
                      lwd = 2.4, type = "b", pch = 16, cex = .6)
  if (k == 1) legend("topleft", bty = "n", cex = .62, lwd = 2.2, col = COL[eg],
                     lty = LTY[eg], seg.len = 2.4,
                     legend = c("nDMSA", "DMSA (linear)", "globaltest", "ORA",
                                "camera", "GSEA/mCSEA", "EWAS+BH"))
  if (k == 2) text(mean(range(d$h)), .45,
                   "every linear method\nstays at its null", cex = .72,
                   col = "grey15", font = 3)
}
dev.off()
cat("wrote fig7, fig8, fig9\n")
