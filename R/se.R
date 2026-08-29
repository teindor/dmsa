## ===========================================================================
## SummarizedExperiment INTEROPERABILITY
##
## Bioconductor's core container for assay + phenotype data is
## SummarizedExperiment. dmsa's own functions take a samples x probes matrix
## and a separate data.frame of phenotypes, because that is what the model
## fitting needs; an SE holds the transpose (features x samples) plus colData.
## This file is the one place that conversion happens, so the transpose lives
## in exactly one function rather than in every user's script - which is where
## a silent samples/probes mix-up would otherwise come from.
## ===========================================================================

#' Split a SummarizedExperiment into the pieces DMSA tests
#'
#' DMSA models probes as the multivariate response, so it wants a
#' \code{samples x probes} matrix and a \code{data.frame} of phenotypes with
#' rows in the same order. A \code{SummarizedExperiment} stores the assay the
#' other way round - \code{features x samples} - with the phenotypes in
#' \code{colData}. This function performs that transposition once and returns
#' the two objects DMSA takes, so the orientation is not something every
#' analysis script has to get right on its own.
#'
#' Probe names are taken from the assay's rownames and become the column names
#' of \code{M}, which is what \code{dmsa_align()} matches its direction calls
#' against. An assay with no rownames is an error rather than a warning: with
#' no probe identifiers there is nothing to align to, and a positional match
#' would be a guess.
#'
#' @param se A \code{SummarizedExperiment}, features (probes) in rows and
#'   samples in columns, as returned by resources such as
#'   \code{recountmethylation}.
#' @param assay Which assay to use: a name or a positive integer. Defaults to
#'   the first. On a methylation object this is normally the beta or M-value
#'   matrix; DMSA standardises probes before aggregating, so either scale
#'   works, but M-values are the better-behaved choice for linear models.
#' @param probes Optional character vector restricting the probes returned, in
#'   the order given. Probes absent from the assay are an error, since a
#'   silently shortened panel changes which family a unit is corrected in.
#' @return A list with two elements: \code{M}, a \code{samples x probes}
#'   numeric matrix, and \code{data}, a \code{data.frame} of the object's
#'   \code{colData} with rows in the same order as \code{M}. Pass them to
#'   \code{dmsa_triangulate()} or \code{dmsa_frame()}.
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   set.seed(1)
#'   counts <- matrix(rnorm(20 * 30), nrow = 20,
#'                    dimnames = list(sprintf("cg%04d", 1:20), NULL))
#'   se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(mval = counts),
#'     colData = data.frame(age = rnorm(30), sex = rep(c("F", "M"), 15)))
#'   parts <- dmsa_se(se)
#'   dim(parts$M)        # samples x probes, the orientation DMSA wants
#'   head(parts$data, 3)
#' }
#' @export
dmsa_se <- function(se, assay = 1L, probes = NULL) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
    stop("SummarizedExperiment is required to use dmsa_se()", call. = FALSE)
  if (!methods::is(se, "SummarizedExperiment"))
    stop("`se` must be a SummarizedExperiment, not ", class(se)[1],
         call. = FALSE)

  nass <- length(SummarizedExperiment::assays(se))
  if (is.character(assay)) {
    nms <- SummarizedExperiment::assayNames(se)
    if (!assay %in% nms)
      stop("no assay named '", assay, "'. Available: ",
           paste(nms, collapse = ", "), call. = FALSE)
  } else {
    assay <- as.integer(assay)
    if (is.na(assay) || assay < 1L || assay > nass)
      stop("`assay` must name an assay or index one of the ", nass,
           " present", call. = FALSE)
  }

  a <- SummarizedExperiment::assay(se, assay)
  if (is.null(rownames(a)))
    stop("the assay has no rownames, so its probes cannot be identified. ",
         "DMSA matches direction calls to probes by name; set rownames(se) ",
         "to the probe IDs before calling.", call. = FALSE)

  if (!is.null(probes)) {
    probes <- as.character(probes)
    miss <- setdiff(probes, rownames(a))
    if (length(miss))
      stop(length(miss), " requested probe(s) are not in the assay, e.g. ",
           paste(utils::head(miss, 3), collapse = ", "), call. = FALSE)
    a <- a[probes, , drop = FALSE]
  }

  M <- t(as.matrix(a))
  storage.mode(M) <- "double"

  dat <- as.data.frame(SummarizedExperiment::colData(se),
                       stringsAsFactors = FALSE)
  if (nrow(dat) != nrow(M))
    stop("colData has ", nrow(dat), " rows but the assay has ", nrow(M),
         " samples", call. = FALSE)
  rownames(dat) <- NULL

  list(M = M, data = dat)
}

