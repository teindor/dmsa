# ============================================================================
# DMSA REPORT FIGURES
#
# Three panels, and the second and third are produced ONCE PER FIRING SYSTEM
# rather than once per analysis, because a system that did not fire has nothing
# to drill into and a figure that mixes systems hides which one carried what.
#
#   panel 1  every system x outcome, dense arm against sparse arm. One figure.
#   panel 2  gene level inside one firing system, with the permutation
#            max-|z| thresholds drawn. One figure per firing system.
#   panel 3  probe level inside the genes that were named. One figure per
#            firing system.
#
# THE SIGN CONVENTION IN PANEL 3, which is the whole point of the panel.
#
# A probe's fitted slope b_j is an effect on METHYLATION. What the reader wants
# is the effect in a single, biologically consistent direction, and that
# requires the CpG's own methylation-to-expression direction d_j from
# cpgdirection. Probes with one value of d_j must therefore be reflected before
# they can be shown on the same axis as the rest.
#
# `invert` names WHICH d_j gets reflected, and the plot marks every reflected
# probe with an open symbol and states the convention in the subtitle. Nothing
# is silently flipped: a reader can always recover the raw slope by reflecting
# the open symbols back.
# ============================================================================

.dmsa_pal <- function(n = 5) {
  if (requireNamespace("viridisLite", quietly = TRUE))
    viridisLite::viridis(n, end = .92)
  else grDevices::hcl.colors(n, "viridis")
}

#' Panel 1: dense arm against sparse arm, every system
#'
#' Distance from the diagonal is how concentrated a system's response is: a
#' system null on the pooled arm and strong on the max-gene arm carries its
#' signal in one gene.
#'
#' @param systems data.frame with columns \code{system}, \code{outcome},
#'   \code{p_dense}, \code{p_sparse}; optionally \code{primary} (logical) and
#'   \code{top_gene}.
#' @param alpha Gate level. The gate is \code{min(1, 2*min(p_dense, p_sparse))},
#'   so the boundary drawn is at \code{alpha/2} on each arm.
#' @param file Optional png path. If \code{NULL} draws to the current device.
#' @param width,height,res png dimensions.
#' @return Invisibly, the plotted data frame.
#' @examples
#' ## Distance from the diagonal says how concentrated a system's response is:
#' ## strong on the sparse arm only means the signal sits in a single gene.
#' sys <- data.frame(
#'   system   = c("HPA axis", "Oxytocin", "Serotonin", "Inflammation"),
#'   outcome  = "PTSD severity",
#'   p_dense  = c(0.0004, 0.08, 0.41, 0.62),
#'   p_sparse = c(0.0012, 0.0003, 0.33, 0.55),
#'   primary  = c(TRUE, TRUE, FALSE, FALSE),
#'   top_gene = c("NR3C1", "OXTR", "TPH2", "IL6"))
#' dmsa_plot_systems(sys, file = tempfile(fileext = ".png"))
#' @export
dmsa_plot_systems <- function(systems, alpha = .05, file = NULL,
                              width = 1500, height = 1400, res = 190) {
  s <- as.data.frame(systems)
  for (v in c("p_dense", "p_sparse")) if (!v %in% names(s))
    stop("systems needs a '", v, "' column", call. = FALSE)
  if (is.null(s$primary)) s$primary <- FALSE
  if (is.null(s$outcome)) s$outcome <- ""
  VC <- .dmsa_pal(5)
  s$ld <- -log10(pmax(s$p_dense, 1e-6)); s$ls <- -log10(pmax(s$p_sparse, 1e-6))
  if (!is.null(file)) { grDevices::png(file, width, height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE) }
  graphics::par(mar = c(4.4, 4.4, 3.4, 1), mgp = c(2.6, .7, 0))
  lim <- c(0, max(s$ld, s$ls, -log10(alpha / 2)) * 1.14)
  gt <- -log10(alpha / 2)
  plot(NA, xlim = lim, ylim = lim,
       xlab = expression(-log[10] ~ italic(p) ~ ", pooled (dense) arm"),
       ylab = expression(-log[10] ~ italic(p) ~ ", max-gene (sparse) arm"),
       main = "Systems: how, not just whether")
  graphics::polygon(c(gt, lim[2], lim[2], gt), c(0, 0, lim[2], lim[2]),
                    col = grDevices::adjustcolor(VC[2], .07), border = NA)
  graphics::polygon(c(0, lim[2], lim[2], 0), c(gt, gt, lim[2], lim[2]),
                    col = grDevices::adjustcolor(VC[2], .07), border = NA)
  graphics::abline(v = gt, h = gt, lty = 2, col = VC[2])
  graphics::abline(0, 1, col = "grey85")
  ou <- unique(s$outcome); pchs <- c(19, 17, 15, 18)[seq_along(ou)]
  graphics::points(s$ld, s$ls, pch = pchs[match(s$outcome, ou)],
                   cex = ifelse(s$primary, 1.7, .95),
                   col = ifelse(s$primary, VC[1], "grey58"))
  fired <- which(pmin(s$p_dense, s$p_sparse) < alpha / 2)
  if (length(fired)) {
    lab <- paste0(substr(s$system[fired], 1, 22),
                  if (!is.null(s$top_gene)) paste0("  ", s$top_gene[fired]) else "")
    graphics::text(s$ld[fired], s$ls[fired], lab, pos = 4, cex = .6, font = 2,
                   col = VC[1], xpd = NA)
  }
  if (length(ou) > 1 || any(nzchar(ou)))
    graphics::legend("bottomright", bty = "n", cex = .68,
                     pch = pchs[seq_along(ou)], col = "grey45", legend = ou)
  graphics::mtext("distance from the diagonal = how concentrated the signal is",
                  3, .2, cex = .6, col = "grey35")
  invisible(s)
}

#' Panel 2: gene level inside one system
#'
#' @param genes data.frame for ONE system and outcome: \code{gene}, \code{z},
#'   \code{n_probes}, optionally \code{p_adj} (max-T adjusted) and
#'   \code{selected}.
#' @param null_max Optional numeric vector of the permutation null of
#'   \code{max|z|} over the same genes; its .95 and .99 quantiles are drawn as
#'   the multiplicity-adjusted thresholds.
#' @param title System label.
#' @inheritParams dmsa_plot_systems
#' @return Called for its side effect of drawing a plot on the active
#'   graphics device. Returns \code{NULL} invisibly.
#' @examples
#' set.seed(1)
#' genes <- data.frame(gene = paste0("G", 1:8), n_probes = sample(3:9, 8, TRUE),
#'                     z = c(rnorm(7), 3.6))
#' genes$p_adj <- 2 * pnorm(-abs(genes$z))
#' ## the permutation null of max|z| over this same family of genes
#' null_max <- apply(matrix(rnorm(200 * 8), 200, 8), 1, function(r) max(abs(r)))
#' dmsa_plot_genes(genes, null_max = null_max, title = "HPA axis",
#'                 file = tempfile(fileext = ".png"))
#' @export
dmsa_plot_genes <- function(genes, null_max = NULL, title = "", alpha = .05,
                            file = NULL, width = 1500, height = 1500,
                            res = 190) {
  g <- as.data.frame(genes)
  n_all <- nrow(g)
  g <- g[is.finite(g$z), , drop = FALSE]
  g <- g[order(g$z), , drop = FALSE]
  VC <- .dmsa_pal(5)
  sel <- if (!is.null(g$selected)) as.logical(g$selected) %in% TRUE else
    if (!is.null(g$p_adj)) is.finite(g$p_adj) & g$p_adj < alpha else
      rep(FALSE, nrow(g))
  if (!is.null(file)) { grDevices::png(file, width, height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE) }
  graphics::par(mar = c(4.4, 7.4, 3.6, 1.2), mgp = c(2.6, .6, 0))
  xr <- range(c(g$z, if (!is.null(null_max)) c(-1, 1) *
                  stats::quantile(null_max, .99, na.rm = TRUE)), na.rm = TRUE)
  xr <- xr + c(-.6, .6)
  plot(NA, xlim = xr, ylim = c(.5, nrow(g) + .5), yaxt = "n",
       xlab = "gene-level aligned z", ylab = "",
       main = paste0("Genes: ", substr(title, 1, 42)))
  if (n_all > nrow(g))
    graphics::mtext(sprintf("%d of %d genes testable; the rest have no probe with a direction call",
                            nrow(g), n_all), side = 1, line = 3.0, cex = .58,
                    col = "grey35")
  if (!is.null(null_max)) {
    q95 <- stats::quantile(null_max, .95, na.rm = TRUE)
    q99 <- stats::quantile(null_max, .99, na.rm = TRUE)
    graphics::polygon(c(-q95, q95, q95, -q95),
                      c(.4, .4, nrow(g) + .6, nrow(g) + .6),
                      col = grDevices::adjustcolor("grey60", .09), border = NA)
    graphics::abline(v = c(-q95, q95), lty = 2, col = "grey45")
    graphics::abline(v = c(-q99, q99), lty = 3, col = VC[2])
    graphics::mtext("dashed = 95% of the permutation null of max|z|; dotted = 99%",
                    3, .2, cex = .56, col = "grey35")
  }
  graphics::abline(v = 0, col = "grey40")
  graphics::axis(2, at = seq_len(nrow(g)), labels = g$gene, las = 1,
                 cex.axis = .62, tick = FALSE)
  graphics::segments(0, seq_len(nrow(g)), g$z, seq_len(nrow(g)),
                     col = "grey86", lwd = 1.6)
  graphics::points(g$z, seq_len(nrow(g)), pch = 19,
                   cex = ifelse(sel, 1.4, .55),
                   col = ifelse(sel, VC[1], "grey62"))
  if (any(sel)) {
    lb <- paste0(g$gene[sel], if (!is.null(g$p_adj))
      sprintf("  p = %.4f", g$p_adj[sel]) else "")
    ## label toward whichever side has room, so it never lands on the axis
    side <- ifelse(g$z[sel] < mean(xr), 4, 2)
    graphics::text(g$z[sel], which(sel), paste0(" ", lb, " "), pos = side,
                   cex = .64, font = 2, col = VC[1], xpd = NA)
  }
  invisible(g)
}

#' Panel 3: probe level, with the CpG-to-expression direction applied
#'
#' The panel that declares a direction. Each probe's slope is reflected
#' according to its \code{cpgdirection} call so that every probe is read on one
#' axis; reflected probes are drawn as open symbols and the convention is
#' printed under the title.
#'
#' @param probes data.frame for ONE system and outcome, restricted to the genes
#'   worth showing: \code{gene}, \code{probe}, \code{b}, \code{se}, \code{d}
#'   (the cpgdirection call, +1 or -1); optionally \code{p} and \code{selected}.
#' @param invert Which \code{d} is reflected onto the common axis: \code{"+1"}
#'   (default) reflects the probes whose methylation predicts HIGHER expression,
#'   \code{"-1"} reflects the others, \code{"none"} shows raw slopes.
#'
#'   With \code{"+1"} the axis reads: \strong{positive = the exposure moves
#'   methylation in the direction that LOWERS this gene's expression}. Every
#'   d = -1 probe is shown as fitted (more methylation, less expression) and
#'   every d = +1 probe is reflected so that it means the same thing.
#' @param gene_summary Optional data.frame with \code{gene} and \code{z}, and
#'   optionally \code{p_adj}: the pooled gene-level result, drawn as a diamond
#'   on its own row. This is usually the point of the figure - individual
#'   probes are rarely significant on their own, and the evidence lives in
#'   their agreement.
#' @param axis_label Overrides the x-axis label.
#' @param title System label.
#' @inheritParams dmsa_plot_systems
#' @return Called for its side effect of drawing a plot on the active
#'   graphics device. Returns \code{NULL} invisibly.
#' @examples
#' ## Every probe is read on one axis. A CpG whose methylation predicts HIGHER
#' ## expression (d = +1) is reflected, so positive always means the exposure
#' ## moved methylation the way that LOWERS this gene's expression.
#' pr <- data.frame(
#'   gene  = "NR3C1",
#'   probe = paste0("cg", 10001:10006),
#'   b     = c(-0.021, -0.014, 0.018, -0.009, 0.026, -0.017),
#'   se    = c(0.008, 0.007, 0.009, 0.006, 0.010, 0.008),
#'   d     = c(-1, -1, 1, -1, 1, -1))
#' dmsa_plot_probes(pr, gene_summary = data.frame(gene = "NR3C1", z = -2.9),
#'                  title = "HPA axis", file = tempfile(fileext = ".png"))
#' @export
dmsa_plot_probes <- function(probes, invert = c("+1", "-1", "none"),
                             gene_summary = NULL,
                             title = "", axis_label = NULL, alpha = .05,
                             file = NULL, width = 1700, height = 1400,
                             res = 190) {
  invert <- match.arg(invert)
  p <- as.data.frame(probes)
  for (v in c("gene", "probe", "b", "se", "d")) if (!v %in% names(p))
    stop("probes needs a '", v, "' column", call. = FALSE)
  p <- p[is.finite(p$b) & is.finite(p$se), , drop = FALSE]
  p$d <- sign(as.numeric(p$d))
  flip <- switch(invert, "+1" = p$d > 0, "-1" = p$d < 0,
                 "none" = rep(FALSE, nrow(p)))
  flip[is.na(flip)] <- FALSE
  p$shown <- ifelse(flip, -p$b, p$b)
  p$flipped <- flip
  p <- p[order(p$gene, p$shown), , drop = FALSE]
  VC <- .dmsa_pal(5)
  ## Probes are drawn PROMINENTLY whether or not each one clears its own
  ## multiplicity-adjusted threshold. Almost none ever does; the gene is
  ## significant because its probes agree, and greying them out would tell the
  ## reader the opposite of the truth. Individually significant probes are
  ## marked separately.
  sel <- if (!is.null(p$p)) is.finite(p$p) & p$p < alpha else rep(FALSE, nrow(p))
  if (!is.null(file)) { grDevices::png(file, width, height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE) }
  extra <- if (!is.null(gene_summary)) 1.6 else 0
  graphics::par(mar = c(6.4, 8.6, 4.2, 2.6), mgp = c(2.6, .6, 0))
  lo <- p$shown - 1.96 * p$se; hi <- p$shown + 1.96 * p$se
  xr <- range(c(lo, hi, 0), na.rm = TRUE); xr <- xr + c(-.08, .08) * diff(xr)
  ttl <- if (invert == "none") "Probes: raw slopes" else
    "Probes, on one direction"
  plot(NA, xlim = xr, ylim = c(.5 - extra, nrow(p) + .5), yaxt = "n",
       xlab = if (!is.null(axis_label)) axis_label else
         if (invert == "none") "slope on methylation" else
           "effect  <-  raises expression tone   |   lowers expression tone  ->",
       ylab = "", cex.main = 1.0,
       main = paste0(ttl, if (nzchar(title)) paste0("  -  ", substr(title, 1, 40))))
  graphics::abline(v = 0, col = "grey40")
  ## gene separators and gene labels
  gb <- cumsum(rle(as.character(p$gene))$lengths)
  if (length(gb) > 1)
    graphics::abline(h = utils::head(gb, -1) + .5, col = "grey90")
  gl <- rle(as.character(p$gene))
  mids <- gb - gl$lengths / 2 + .5
  graphics::mtext(gl$values, side = 2, at = mids, line = 5.6, las = 1,
                  cex = .68, font = 2, col = "grey20")
  graphics::axis(2, at = seq_len(nrow(p)), labels = p$probe, las = 1,
                 cex.axis = .55, tick = FALSE)
  graphics::segments(lo, seq_len(nrow(p)), hi, seq_len(nrow(p)),
                     col = VC[1], lwd = 1.9)
  ## FILLED = shown as fitted; OPEN = reflected by its cpgdirection call
  graphics::points(p$shown, seq_len(nrow(p)),
                   pch = ifelse(p$flipped, 21, 19), cex = 1.15,
                   bg = "white", lwd = 1.8, col = VC[1])
  if (any(sel)) graphics::points(p$shown[sel], which(sel), pch = 8, cex = .8,
                                 col = VC[2])
  graphics::axis(4, at = seq_len(nrow(p)),
                 labels = ifelse(p$d > 0, "+1", "-1"), las = 1, cex.axis = .55,
                 tick = FALSE, col.axis = "grey40")
  graphics::mtext("d", side = 4, line = 1.5, at = nrow(p) + .5, cex = .62,
                  col = "grey40", las = 1)

  ## pooled gene-level result, which is where the evidence actually is
  if (!is.null(gene_summary)) {
    gs <- as.data.frame(gene_summary)
    gs <- gs[gs$gene %in% unique(p$gene), , drop = FALSE]
    if (nrow(gs)) {
      agree <- tapply(sign(p$shown), p$gene, function(v) max(table(v)) / length(v))
      lab <- paste0(gs$gene, ": pooled z = ", sprintf("%+.2f", gs$z),
                    if (!is.null(gs$p_adj)) sprintf(",  maxT p = %.4f", gs$p_adj) else "",
                    sprintf(",  %.0f%% of probes agree in sign",
                            100 * agree[gs$gene]))
      graphics::abline(h = .5 - extra * .15, col = "grey88")
      ## one diamond PER summarised gene - a multi-gene panel used to draw
      ## only the first gene's pooled result and silently drop the rest
      for (gi in seq_len(nrow(gs))) {
        yy <- .5 - extra / 2 - (gi - 1L) * max(.9, extra * .35)
        dx <- mean(p$shown[p$gene == gs$gene[gi]])
        graphics::points(dx, yy, pch = 18, cex = 2.2, col = VC[2])
        ## label away from the diamond so the two never overlap
        graphics::text(dx, yy, paste0(" ", lab[gi], " "),
                       pos = if (dx > mean(xr)) 2 else 4, cex = .66, font = 2,
                       col = VC[2], xpd = NA)
      }
    }
  }

  nflip <- sum(p$flipped)
  sub <- if (invert == "none") "raw fitted slopes; no reflection applied" else
    if (nflip == 0)
      sprintf("no probe needed reflecting: every probe here has d = %s",
              if (invert == "+1") "-1" else "+1")
    else sprintf("open circles (%d of %d): cpgdirection d = %s, slope REFLECTED; filled: shown as fitted",
                 nflip, nrow(p), invert)
  graphics::mtext(sub, 3, .25, cex = .58, col = "grey25")
  graphics::mtext(sprintf(
    "bars 95%% CI   %d probes, %d gene(s)%s",
    nrow(p), length(unique(p$gene)),
    if (any(sel)) sprintf("   * = probe clears maxT on its own (%d)", sum(sel))
    else "   no single probe clears maxT alone - the evidence is the agreement"),
    side = 1, line = 4.6, cex = .6, col = "grey30")
  invisible(p)
}

#' Build the whole report: panel 1 once, panels 2 and 3 per firing system
#'
#' @param systems,genes,probes Data frames as described in the panel functions.
#'   \code{genes} and \code{probes} must carry \code{system_id} and
#'   \code{outcome} so they can be split.
#' @param nulls Optional named list of permutation \code{max|z|} null vectors,
#'   named \code{paste(system_id, outcome, sep = "|")}.
#' @param outdir Directory for the png files.
#' @param genes_shown Which genes go into panel 3: \code{"selected"} (default)
#'   or \code{"top"} (the single strongest) or \code{"all"}.
#' @inheritParams dmsa_plot_systems
#' @inheritParams dmsa_plot_probes
#' @return Invisibly, a data frame of the files written.
#' @examples
#' ## one system fires, so it gets a gene panel and a probe panel of its own
#' systems <- data.frame(system_id = 1:2, system = c("HPA axis", "Oxytocin"),
#'                       outcome = "anx", p_dense = c(0.004, 0.40),
#'                       p_sparse = c(0.02, 0.55), top_gene = c("NR3C1", "OXTR"))
#' genes <- data.frame(system_id = 1L, outcome = "anx",
#'                     gene = c("NR3C1", "FKBP5", "CRH"), z = c(3.4, -1.1, 0.6),
#'                     n_probes = c(6L, 4L, 3L), p_adj = c(0.008, 0.51, 0.88),
#'                     selected = c(TRUE, FALSE, FALSE))
#' probes <- data.frame(system_id = 1L, outcome = "anx", gene = "NR3C1",
#'                      probe = paste0("cg", 1:6), se = rep(.08, 6),
#'                      b = c(.21, .18, .09, -.14, -.20, .16),
#'                      d = c(-1, -1, -1, 1, 1, -1))
#'
#' out <- dmsa_report_panels(systems, genes, probes, outdir = tempfile())
#' out[, c("panel", "system_id", "outcome")]
#' @export
dmsa_report_panels <- function(systems, genes, probes, nulls = NULL, outdir = ".",
                        alpha = .05, invert = "+1",
                        genes_shown = c("selected", "top", "all")) {
  genes_shown <- match.arg(genes_shown)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  s <- as.data.frame(systems); g <- as.data.frame(genes)
  p <- as.data.frame(probes)
  out <- list()
  f1 <- file.path(outdir, "panel1_systems.png")
  dmsa_plot_systems(s, alpha = alpha, file = f1)
  out[[1]] <- data.frame(panel = 1, system_id = NA, outcome = NA, file = f1)
  fired <- s[pmin(s$p_dense, s$p_sparse) < alpha / 2, , drop = FALSE]
  if (!nrow(fired)) {
    message("no system fired at alpha = ", alpha, "; panels 2 and 3 skipped")
    return(invisible(do.call(rbind, out)))
  }
  for (i in seq_len(nrow(fired))) {
    sid <- fired$system_id[i]; oc <- fired$outcome[i]
    key <- paste(sid, oc, sep = "|")
    tag <- paste0("sys", sid, "_", gsub("[^A-Za-z0-9]", "", oc))
    gi <- g[g$system_id == sid & g$outcome == oc, , drop = FALSE]
    if (nrow(gi)) {
      f2 <- file.path(outdir, paste0("panel2_", tag, "_genes.png"))
      dmsa_plot_genes(gi, null_max = if (!is.null(nulls)) nulls[[key]] else NULL,
                      title = paste(fired$system[i], "x", oc), alpha = alpha,
                      file = f2)
      out[[length(out) + 1L]] <- data.frame(panel = 2, system_id = sid,
                                            outcome = oc, file = f2)
    }
    pi <- p[p$system_id == sid & p$outcome == oc, , drop = FALSE]
    keep <- switch(genes_shown,
      selected = if (!is.null(gi$selected)) gi$gene[as.logical(gi$selected)] else
        if (!is.null(gi$p_adj)) gi$gene[gi$p_adj < alpha] else character(0),
      top = gi$gene[which.max(abs(gi$z))],
      all = gi$gene)
    if (!length(keep)) keep <- gi$gene[which.max(abs(gi$z))]   # always show one
    pi <- pi[pi$gene %in% keep, , drop = FALSE]
    if (nrow(pi)) {
      f3 <- file.path(outdir, paste0("panel3_", tag, "_probes.png"))
      dmsa_plot_probes(pi, invert = invert,
                       gene_summary = gi[gi$gene %in% keep, , drop = FALSE],
                       title = paste(fired$system[i], "x", oc), alpha = alpha,
                       file = f3)
      out[[length(out) + 1L]] <- data.frame(panel = 3, system_id = sid,
                                            outcome = oc, file = f3)
    }
  }
  res <- do.call(rbind, out)
  message("wrote ", nrow(res), " figure(s) to ", normalizePath(outdir))
  invisible(res)
}
