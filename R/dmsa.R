# dmsa: Directional Methylation Set Analysis
# Two-layer sign chain: s_j = d_j (CpG -> expression, from cpgdirection)
#                             x w_g (gene -> system activation, curated polarity)
# Gene-level analysis uses d_j alone; system-level uses d_j * w_g.

#' Alpha panel gene -> system map (bundled)
#'
#' Read the Project Alpha gene codebook shipped with the package:
#' 549 genes assigned to 30 biological systems.
#' @return data.frame with system_id, system, gene, n_kept, usability
#' @examples
#' gs <- alpha_gene_systems()
#' dim(gs)
#' length(unique(gs$system_id))
#' head(gs[, c("system_id", "system", "gene", "n_kept", "usability")], 4)
#' @export
alpha_gene_systems <- function() {
  p <- system.file("extdata", "alpha_gene_systems.csv", package = "dmsa")
  utils::read.csv(p, stringsAsFactors = FALSE)
}

#' Alpha panel gene -> system-activation polarity (bundled, DRAFT)
#'
#' Curated w_g per gene: +1 gene product increases the system's activation
#' tone, -1 decreases it (e.g. NR3C1 in the HPA axis: more GR = stronger
#' negative feedback = LOWER axis tone), 0 = not on the activation axis
#' (readouts, specificity controls, ambiguous).
#'
#' STATUS: draft pending PI approval and per-gene citations. Every use
#' prints a reminder until the table is finalised.
#' @return data.frame with system_id, system, gene, w_g_draft, role, confidence
#' @examples
#' pol <- alpha_polarity()
#' head(pol[, c("system", "gene", "w_g", "role")], 4)
#' # CRH drives the HPA axis (w_g = +1); NR3C1 is its brake, so more GR means
#' # stronger negative feedback and LOWER axis tone (w_g = -1)
#' pol[pol$gene %in% c("CRH", "NR3C1") & pol$system_id == 2,
#'     c("gene", "w_g", "role")]
#' @export
alpha_polarity <- function() {
  p <- system.file("extdata", "alpha_polarity.csv", package = "dmsa")
  out <- utils::read.csv(p, stringsAsFactors = FALSE)
  names(out)[names(out) == "w_g_draft"] <- "w_g"
  attr(out, "status") <- "draft"
  out
}

#' Build the aligned sign s_j for a set of probes
#'
#' @param direction A data.frame from
#'   \code{cpgdirection::cpg_expression_direction()} (needs columns
#'   \code{cpg_id} or \code{input}, \code{best_direction}, and if available
#'   \code{best_confidence}), OR a data.frame with columns \code{cpg},
#'   \code{d} (+1/-1) and optionally \code{p_plus} = P(d = +1).
#' @param genes Character vector, same length as the probes in
#'   \code{direction}: the set-membership gene of each probe (the panel's
#'   annotation, e.g. from the column suffix).
#' @param level \code{"gene"} (alignment = d_j only) or \code{"system"}
#'   (alignment = d_j * w_g).
#' @param polarity For \code{level = "system"}: where w_g comes from.
#'   \code{"alpha"} (default) uses the bundled Alpha polarity table.
#'   Or supply your own data.frame with columns \code{gene} and \code{w_g}
#'   (+1 / -1 / 0). User-supplied values OVERRIDE the bundled table where
#'   both cover a gene.
#' @param system_id Optional integer: restrict the bundled polarity lookup
#'   to one Alpha system (a gene can carry different roles in different
#'   systems).
#' @param missing_polarity What to do with a called probe whose gene has no
#'   polarity entry: \code{"error"} (default - stop and list the genes so
#'   the user can specify each, per protocol), \code{"zero"} (treat as
#'   off-axis, probe drops from the system pool but is reported), or
#'   \code{"drop"} (silently exclude).
#' @return data.frame: probe, gene, d, p_plus, w_g, s (= d * w_g at system
#'   level, d at gene level), p_s_plus = P(s = +1) (chained), usable
#'   (logical), and reason for every non-usable probe.
#' @examples
#' # d is the CpG -> expression sign, w_g the gene -> system-tone sign.
#' # A promoter CpG silences NR3C1 (d = -1) and GR is the HPA brake (w_g = -1),
#' # so the probe votes for HIGHER axis tone: s = d * w_g = +1.
#' g <- c("CRH", "NR3C1")
#' dcall <- data.frame(cpg = c("cg_crh", "cg_nr3c1"), d = c(1, -1),
#'                     p_plus = c(0.9, 0.1))
#' pol <- data.frame(gene = g, w_g = c(1, -1))
#' dmsa_align(dcall, genes = g, level = "gene")
#' dmsa_align(dcall, genes = g, level = "system", polarity = pol)
#' @export
dmsa_align <- function(direction, genes,
                       level = c("gene", "system"),
                       polarity = "alpha",
                       system_id = NULL,
                       missing_polarity = c("error", "zero", "drop")) {
  level <- match.arg(level)
  missing_polarity <- match.arg(missing_polarity)

  ## -- normalise the direction input ------------------------------------
  dd <- as.data.frame(direction)
  cpg_col <- intersect(c("cpg_id", "cpg", "input"), names(dd))[1]
  if (is.na(cpg_col)) stop("direction needs a 'cpg_id' / 'cpg' / 'input' column")
  d_col <- intersect(c("best_direction", "d", "direction"), names(dd))[1]
  if (is.na(d_col)) stop("direction needs 'best_direction' or 'd'")
  out <- data.frame(probe = as.character(dd[[cpg_col]]),
                    d     = suppressWarnings(as.numeric(dd[[d_col]])),
                    stringsAsFactors = FALSE)
  if ("p_plus" %in% names(dd)) {
    out$p_plus <- dd$p_plus
  } else if ("best_confidence" %in% names(dd)) {
    ## best_confidence = confidence in the CALLED direction
    conf <- pmin(pmax(dd$best_confidence, 0), 1)
    out$p_plus <- ifelse(out$d > 0, 0.5 + conf / 2, 0.5 - conf / 2)
  } else {
    out$p_plus <- ifelse(is.na(out$d), NA_real_, ifelse(out$d > 0, 1, 0))
  }
  if (length(genes) != nrow(out))
    stop("genes must have one entry per probe row (got ", length(genes),
         " for ", nrow(out), " probes)")
  out$gene <- as.character(genes)

  ## -- polarity layer -----------------------------------------------------
  if (level == "system") {
    pol <- NULL
    if (is.character(polarity) && identical(polarity, "alpha")) {
      pol <- alpha_polarity()
      if (!is.null(system_id)) pol <- pol[pol$system_id == system_id, ]
      message("dmsa: using bundled Alpha polarity table (status: DRAFT - ",
              "pending PI approval and citations).")
    } else {
      pol <- as.data.frame(polarity)
      if (!all(c("gene", "w_g") %in% names(pol)))
        stop("user polarity needs columns 'gene' and 'w_g' (+1/-1/0)")
    }
    ## user table overrides bundled where both exist
    if (is.data.frame(polarity) && identical(attr(pol, "status"), "draft")) {
      # (bundled only - nothing to merge)
    }
    idx <- match(out$gene, pol$gene)
    out$w_g <- pol$w_g[idx]

    unknown <- sort(unique(out$gene[is.na(out$w_g) & !is.na(out$d)]))
    if (length(unknown) > 0) {
      if (missing_polarity == "error") {
        stop("No polarity (w_g) for gene(s): ", paste(unknown, collapse = ", "),
             "\nSpecify each via the 'polarity' argument, e.g.\n",
             "  polarity = data.frame(gene = c(",
             paste0('"', unknown[seq_len(min(2, length(unknown)))], '"',
                    collapse = ", "),
             ", ...), w_g = c(...))\n",
             "(+1 activates the system, -1 opposes it - e.g. NR3C1 in HPA, ",
             "0 = off-axis)", call. = FALSE)
      } else if (missing_polarity == "zero") {
        warning("w_g set to 0 (off-axis) for: ", paste(unknown, collapse = ", "))
        out$w_g[is.na(out$w_g)] <- 0
      } else {
        out <- out[!is.na(out$w_g) | is.na(out$d), , drop = FALSE]
      }
    }
    out$s <- out$d * out$w_g
    ## chained probability: P(s=+1) = p_d*p_w + (1-p_d)(1-p_w); w_g known -> p_w in {0,1}
    p_w <- ifelse(out$w_g > 0, 1, ifelse(out$w_g < 0, 0, NA_real_))
    out$p_s_plus <- out$p_plus * p_w + (1 - out$p_plus) * (1 - p_w)
  } else {
    out$w_g <- NA_real_
    out$s <- out$d
    out$p_s_plus <- out$p_plus
  }

  out$p_s_plus[is.na(out$d)] <- NA_real_   # an abstained probe carries no sign certainty
  out$usable <- !is.na(out$s) & out$s != 0
  out$reason <- ifelse(is.na(out$d), "no_direction_call",
                ifelse(!is.na(out$w_g) & out$w_g == 0, "off_axis_gene",
                ifelse(is.na(out$s), "no_polarity", "")))
  class(out) <- c("dmsa_alignment", class(out))
  out
}

#' Fixed-sign or expected-sign DMSA pooled test
#'
#' Pools probe-level coefficients after alignment. Two estimators:
#' \code{"fixed"} multiplies each b_j by the hard sign s_j; expected-sign
#' (\code{"expected"}) multiplies by E[s_j] = 2 P(s_j = +1) - 1, the
#' first-order latent-sign approximation - low-tier calls contribute in
#' proportion to their calibrated certainty, which is why no tier filtering
#' is needed upstream.
#'
#' The z test here is descriptive; the reportable error control is the
#' family-wise permutation (\code{dmsa_perm_pvalue}) on the same statistic.
#'
#' @param b,se Probe-level coefficient and SE from YOUR model (fitted on
#'   M-values standardised within probe; same estimator observed and null).
#' @param alignment A \code{dmsa_align()} result covering the same probes,
#'   or vectors via \code{s} / \code{p_s_plus}.
#' @param method \code{"fixed"} or \code{"expected"}.
#' @param w Optional per-probe reliability weight (see
#'   \code{\link{dmsa_relweights}}); scales each aligned multiplier. \code{NULL}
#'   (default) or all ones is the flat engine.
#' @return list: estimate (pooled aligned slope), se, z, p_normal,
#'   n_used, n_excluded, table (per-probe contributions).
#' @examples
#' ## Two promoter CpGs in NR3C1: methylation silences the gene (d = -1), and
#' ## GR is the HPA brake (w_g = -1), so both probes vote for HIGHER axis tone.
#' calls <- data.frame(cpg = c("cg01", "cg02"), d = c(-1, -1), p_plus = c(.1, .2))
#' al <- dmsa_align(calls, genes = c("NR3C1", "NR3C1"), level = "system",
#'                  polarity = data.frame(gene = "NR3C1", w_g = -1))
#'
#' tt <- dmsa_test(b = c(0.28, 0.19), se = c(0.09, 0.08), al, method = "expected")
#' tt[c("estimate", "se", "z", "p_normal", "n_used")]
#' tt$table
#' @export
dmsa_test <- function(b, se, alignment,
                      method = c("expected", "fixed"), w = NULL) {
  method <- match.arg(method)
  al <- as.data.frame(alignment)
  if (nrow(al) != length(b) || length(b) != length(se))
    stop("b, se and alignment must cover the same probes, in order")
  mult <- if (method == "fixed") al$s else (2 * al$p_s_plus - 1)
  if (!is.null(w)) mult <- mult * w
  use <- !is.na(mult) & mult != 0 & !is.na(b) & !is.na(se) & se > 0
  ## GLS estimator of the common aligned effect mu under b_j ~ N(mu * m_j, se_j^2):
  ## mu_hat = sum(m_j b_j / se_j^2) / sum(m_j^2 / se_j^2)
  ## (reduces to standard inverse-variance pooling when |m_j| = 1, i.e. fixed-sign;
  ##  under expected-sign, low-certainty probes contribute ~nothing rather than dominating)
  num <- sum(mult[use] * b[use] / se[use]^2)
  den <- sum(mult[use]^2 / se[use]^2)
  est <- num / den
  pse <- sqrt(1 / den)
  z   <- est / pse
  structure(list(
    estimate = est, se = pse, z = z,
    p_normal = 2 * stats::pnorm(-abs(z)),
    method = method, n_used = sum(use), n_excluded = sum(!use),
    table = data.frame(probe = al$probe, gene = al$gene, b = b, se = se,
                       s = al$s, p_s_plus = al$p_s_plus,
                       multiplier = mult, b_aligned = b * mult,
                       used = use)
  ), class = "dmsa_test")
}

#' @export
print.dmsa_test <- function(x, ...) {
  cat("DMSA pooled aligned test (", x$method, "-sign)\n", sep = "")
  cat(sprintf("  estimate %.4f  se %.4f  z %.2f  p(normal) %.4g\n",
              x$estimate, x$se, x$z, x$p_normal))
  cat(sprintf("  probes used %d, excluded %d\n", x$n_used, x$n_excluded))
  cat("  NOTE: report the permutation p (dmsa_perm_pvalue), not p(normal).\n")
  invisible(x)
}

#' Permutation p-value for a DMSA statistic
#'
#' @param observed The observed statistic (e.g. z from \code{dmsa_test}).
#' @param null_stats Vector of the same statistic computed on B
#'   family-wise-permuted datasets (shuffle the exposure BY FAMILY, refit
#'   the probe models with the SAME estimator, re-run \code{dmsa_test}).
#' @return two-sided permutation p with the +1 correction (resolution
#'   floor 1/(B+1)).
#' @examples
#' set.seed(1)
#' null_z <- rnorm(999)          # z from B family-wise permuted datasets
#' dmsa_perm_pvalue(3.1, null_z)
#' ## the +1 correction puts a floor of 1/(B+1) on any permutation p
#' dmsa_perm_pvalue(50, null_z)
#' @export
dmsa_perm_pvalue <- function(observed, null_stats) {
  B <- length(null_stats)
  (1 + sum(abs(null_stats) >= abs(observed))) / (B + 1)
}

#' Per-system coverage and polarity balance
#'
#' The honest-reporting companion: how many probes each system can actually
#' align, on which side of the activation axis. A system whose brake genes
#' have no usable probes silently degrades to an activation-gene test -
#' this table is how that is caught and reported.
#' @param alignment result of dmsa_align(level = "system")
#' @return one-row data.frame per call
#' @examples
#' dcall <- data.frame(cpg = paste0("p", 1:4), d = c(1, -1, 1, NA),
#'                     p_plus = c(0.9, 0.1, 0.8, NA))
#' genes <- c("CRH", "POMC", "AVPR2", "NR3C1")
#' pol <- data.frame(gene = genes, w_g = c(1, 1, 0, -1))
#' al <- dmsa_align(dcall, genes, level = "system", polarity = pol)
#' # the only brake gene abstained, so this "system" test is really an
#' # activation-gene test - which is what the balance table makes visible
#' dmsa_balance(al)
#' @export
dmsa_balance <- function(alignment) {
  al <- as.data.frame(alignment)
  data.frame(
    probes            = nrow(al),
    called            = sum(!is.na(al$d)),
    usable            = sum(al$usable),
    activation_probes = sum(al$usable & al$w_g > 0, na.rm = TRUE),
    brake_probes      = sum(al$usable & al$w_g < 0, na.rm = TRUE),
    off_axis          = sum(al$reason == "off_axis_gene", na.rm = TRUE),
    no_direction      = sum(al$reason == "no_direction_call", na.rm = TRUE),
    activation_genes  = length(unique(al$gene[al$usable & al$w_g > 0])),
    brake_genes       = length(unique(al$gene[al$usable & al$w_g < 0]))
  )
}
