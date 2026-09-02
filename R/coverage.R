# ============================================================================
# Spec 15/16: dmsa_coverage() - what the biology asked for vs what the data
# can test, itemised so nothing disappears in silence.
#
# Three levels:
#   system  - genes in the reference vs genes with a usable pair here, plus
#             the size of each level-local correction family (what maxT/minP
#             actually corrects within; it GROWS as coverage grows, which is
#             why adjusted p-values are not comparable across versions)
#   gene    - per-gene pair counts and the reasons pairs were not used
#   pairs   - the full CpG x gene pair ledger, unfiltered
# ============================================================================

#' Coverage report: the reference vs the testable data
#'
#' What the biological reference asked for, what the data at hand can test,
#' and why every discovered CpG x gene pair was or was not used. Only
#' available for frames built on the pair path
#' (\code{direction_source = "cpgdirection"}); a frame on the bundled
#' snapshot has no pair ledger and gets an error saying so.
#'
#' @param frame A \code{dmsa_frame}.
#' @param level \code{"system"} (default): per-system counts - reference
#'   genes, testable genes, untestable genes, measurements, pairs, and the
#'   gene-level correction family size. \code{"gene"}: per-gene pair counts
#'   with the dominant reason unused pairs were abstained. \code{"pairs"}:
#'   the full pair ledger, one row per discovered CpG x gene pair.
#' @return A data.frame (invisible its print). For \code{"system"} the
#'   columns are system_id, system, n_genes_reference, n_genes_testable,
#'   n_genes_untestable, n_measurements, n_pairs_used, n_pairs_discovered,
#'   family_size_gene (= n_genes_testable, the size of the maxT/minP family
#'   at the gene level).
#' @examples
#' \donttest{
#' ## needs cpgdirection; see the vignette for a full run
#' }
#' @export
dmsa_coverage <- function(frame, level = c("system", "gene", "pairs")) {
  stopifnot(inherits(frame, "dmsa_frame"))
  level <- match.arg(level)
  led <- frame$pair_ledger
  if (is.null(led))
    stop("this frame was built from the bundled per-column snapshot ",
         "(direction_source = \"", frame$direction_source %||% "bundled",
         "\"), which carries no pair ledger.\ndmsa_coverage() reports the ",
         "CpG x gene pair path; rebuild the frame with cpgdirection ",
         "installed (direction_source = \"auto\" finds it).", call. = FALSE)

  if (level == "pairs") return(led)

  ## restrict to the systems the frame analyses
  sys_ids <- as.character(unique(frame$map$system_id))
  led_in <- led[as.character(led$system_id) %in% sys_ids |
                  is.na(led$system_id), , drop = FALSE]

  if (level == "gene") {
    sp <- split(led_in, paste(led_in$system_id, led_in$target_gene, sep = "\r"))
    out <- do.call(rbind, lapply(sp, function(g) {
      why <- g$reason_not_used[!g$used]
      data.frame(system_id = as.character(g$system_id[1]),
                 system = g$system[1],
                 gene = g$target_gene[1],
                 n_pairs = nrow(g),
                 n_used = sum(g$used),
                 n_measurements = length(unique(g$measurement_id)),
                 testable = any(g$used),
                 main_reason_not_used = if (all(g$used)) NA_character_
                   else names(sort(table(why), decreasing = TRUE))[1],
                 stringsAsFactors = FALSE)
    }))
    out <- out[order(out$system_id, out$gene), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }

  ## level == "system": reference counts come from the frame's untestable
  ## table plus the tested map - together they ARE the reference restricted
  ## to the analysed systems, by construction (spec 15)
  tested <- unique(frame$map[, c("system_id", "system", "gene")])
  untest <- frame$untestable
  if (is.null(untest)) untest <- tested[0, ]
  allg <- unique(rbind(tested[, c("system_id", "system", "gene")],
                       untest[, c("system_id", "system", "gene")]))
  sp <- split(allg, allg$system_id)
  out <- do.call(rbind, lapply(sp, function(g) {
    sid <- as.character(g$system_id[1])
    t_g <- tested$gene[as.character(tested$system_id) == sid]
    l_g <- led_in[as.character(led_in$system_id) == sid, , drop = FALSE]
    m_g <- frame$map[as.character(frame$map$system_id) == sid, , drop = FALSE]
    data.frame(system_id = sid, system = g$system[1],
               n_genes_reference = length(unique(g$gene)),
               n_genes_testable = length(unique(t_g)),
               n_genes_untestable = length(unique(g$gene)) -
                 length(unique(t_g)),
               n_measurements = length(unique(m_g$column)),
               n_pairs_used = nrow(m_g),
               n_pairs_discovered = nrow(l_g),
               family_size_gene = length(unique(t_g)),
               stringsAsFactors = FALSE)
  }))
  out <- out[order(out$system_id), , drop = FALSE]
  rownames(out) <- NULL
  ## spec 40: polarity coverage per system - how many of the testable genes
  ## carry a signed activation polarity, and how many are unresolved
  .pa <- tryCatch(.rp_polarity_audit(frame), error = function(e) NULL)
  if (!is.null(.pa) && isTRUE(attr(.pa, "has_polarity"))) {
    i <- match(as.character(out$system_id), .pa$system_id)
    out$n_genes_polarity_signed <- .pa$n_polarity_signed[i]
    out$n_genes_polarity_unresolved <- .pa$n_polarity_unresolved[i]
    out$genes_polarity_unresolved <- .pa$genes_unresolved[i]
  }
  for (at in c("qc_excluded_cpgs", "unmapped_cpgs", "cpgs_outside_reference",
               "n_cpgs_submitted", "n_cpgs_mapped", "tissue"))
    attr(out, at) <- attr(led, at)
  class(out) <- c("dmsa_coverage", class(out))
  out
}

#' @export
print.dmsa_coverage <- function(x, ...) {
  cat("DMSA coverage - the reference vs the testable data\n")
  ## reconcile the CpG accounting up front: submitted = mapped + unmapped +
  ## QC-excluded, so a "missing" CpG always has a stated fate
  .ns <- attr(x, "n_cpgs_submitted"); .nm <- attr(x, "n_cpgs_mapped")
  .qc <- length(attr(x, "qc_excluded_cpgs"))
  .un <- length(attr(x, "unmapped_cpgs"))
  .or <- length(attr(x, "cpgs_outside_reference"))
  if (!is.null(.ns))
    cat(sprintf(paste0("  CpGs: %d submitted -> %d mapped to reference ",
                       "genes, %d mapped only outside the reference, ",
                       "%d unmapped, %d QC-excluded\n"),
                .ns, .nm %||% NA_integer_, .or, .un, .qc))
  cat(sprintf("  %d system(s); reference genes %d, testable here %d\n",
              nrow(x), sum(x$n_genes_reference), sum(x$n_genes_testable)))
  if ("n_genes_polarity_unresolved" %in% names(x) &&
      any(x$n_genes_polarity_unresolved > 0, na.rm = TRUE))
    cat(sprintf(paste0("  polarity: %d testable gene(s) have NO resolved ",
                       "activation sign and are weighted 0 in the system ",
                       "score (see genes_polarity_unresolved)\n"),
                sum(x$n_genes_polarity_unresolved, na.rm = TRUE)))
  cat("  family_size_gene is the size of the gene-level maxT/minP family:\n")
  cat("  it grows with coverage, so adjusted p-values are not comparable\n")
  cat("  across versions or datasets with different coverage.\n\n")
  print.data.frame(x, row.names = FALSE)
  invisible(x)
}
