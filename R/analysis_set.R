# The analysis set, saved for publication and validation ---------------------
#
# A reviewer or a replicating lab needs ONE table that answers "which CpGs did
# this analysis actually test, against which genes, under which direction
# calls" - not the full pair ledger (which also carries every pair that was
# NOT used) and not the frame object (which needs R and the package to open).
# This is that table, written as a plain .csv (PI request, 2026-08-29).

#' Save the probes a frame actually analyses, annotated, as a CSV
#'
#' One row per CpG x gene pair that entered the analysis: exactly the probes
#' \code{dmsa_frame()} selected after mapping, probe QC, replicate collapse
#' and any \code{cpgs_include}/\code{cpgs_exclude} filter. Each row carries
#' the annotation the frame resolved for it - gene, system, module, the
#' direction call with its evidence tier - plus genomic coordinates when the
#' bundled cascade covers the probe. Intended for the supplement of a paper
#' and for user validation; the same information, plus every pair that was
#' \emph{not} used and why, is in the pair ledger
#' (\code{frame$pair_ledger}, written by the report as
#' \code{tables/cpg_gene_pair_ledger.csv}).
#'
#' A report run writes this file automatically as
#' \code{tables/analysis_set.csv}; call this yourself to write it anywhere
#' else, or with \code{file = NULL} to get the table without writing.
#'
#' @param frame A \code{\link{dmsa_frame}}.
#' @param file Path of the CSV to write, or \code{NULL} to only return the
#'   table.
#' @return The analysis-set data frame, invisibly when \code{file} is given.
#' @examples
#' set.seed(1)
#' map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
#'                   probe = c("cg01", "cg02"), column = c("cg01", "cg02"),
#'                   best_direction = c(-1, 1), p_plus = c(.1, .9))
#' dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40),
#'                   cg01 = plogis(rnorm(40)), cg02 = plogis(rnorm(40)))
#' fr <- dmsa_frame(dat, map = map, outcomes = "anx", covariates = "cov1",
#'                  random_effects = NULL, B = 49, seed = 1, plots = FALSE,
#'                  tables = FALSE, summary = FALSE, progress = FALSE,
#'                  beep = FALSE, outdir = tempfile())
#' dmsa_save_analysis_set(fr, file = NULL)
#' @export
dmsa_save_analysis_set <- function(frame, file = "analysis_set.csv") {
  if (!inherits(frame, "dmsa_frame"))
    stop("`frame` must be a dmsa_frame", call. = FALSE)
  mp <- frame$map
  ## the columns a publication needs, in reading order, kept only when the
  ## frame's path produced them (the bundled path has no evidence tiers)
  want <- c(probe = "probe", gene = "gene", system = "system",
            module = "module", direction = "best_direction",
            p_plus1 = "p_plus", direction_tier = "direction_tier",
            evidence = "best_evidence", co_effect = "is_coeffect_selected",
            genes_for_this_cpg = "n_targets_selected")
  have <- want[want %in% names(mp)]
  out <- mp[, have, drop = FALSE]
  names(out) <- names(have)
  ## measurement column(s) behind the probe: identical to `probe` unless
  ## replicates were collapsed, in which case every replicate is named
  if (!is.null(frame$measurements) &&
      all(c("cpg_id", "column") %in% names(frame$measurements))) {
    reps <- vapply(out$probe, function(p) paste(
      frame$measurements$column[frame$measurements$cpg_id == p],
      collapse = ";"), character(1))
    if (!all(reps == out$probe | reps == "")) out$data_columns <- unname(reps)
  } else if (!identical(mp$column, mp$probe)) out$data_columns <- mp$column
  ## coordinates from the bundled cascade (hg38) when it covers the probe -
  ## the same source the locus figure uses, for the same reason: it is the
  ## table the gene models are aligned to
  cas <- tryCatch(dmsa_sets(if (is.null(frame$sets_source)) "alpha"
                            else frame$sets_source)$cascade,
                  error = function(e) NULL)
  if (!is.null(cas) && all(c("cpg", "chr", "pos_hg38") %in% names(cas))) {
    i <- match(out$probe, cas$cpg)
    if (any(!is.na(i))) {
      out$chr <- cas$chr[i]
      out$pos_hg38 <- suppressWarnings(as.numeric(cas$pos_hg38[i]))
    }
  }
  ## provenance, constant per frame, so the CSV stands alone
  out$direction_source <- frame$direction_source %||% "bundled"
  if (!is.null(frame$tissue)) out$tissue <- frame$tissue
  if (!is.null(frame$reference_name)) out$reference <- frame$reference_name
  rownames(out) <- NULL
  if (is.null(file)) return(out)
  utils::write.csv(out, file, row.names = FALSE)
  invisible(out)
}
