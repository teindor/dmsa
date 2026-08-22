## ---------------------------------------------------------------------------
## dmsa_import(): one door in from every methylation pipeline
##
## Every mainstream preprocessing pipeline ends in one of a handful of shapes:
##   minfi        GenomicRatioSet / MethylSet (SummarizedExperiment subclasses)
##   sesame       openSesame() beta matrix (probes x samples), or list of SigDF
##   meffil       meffil.normalize.samples() matrix (probes x samples)
##   ewastools    dont_normalize()/betas matrix (probes x samples)
##   ChAMP        myNorm matrix (probes x samples) + pd data.frame
##   RnBeads      RnBSet (meth() matrix, probes x samples)
##   methylprep   beta CSV export (probes in the first column or rownames)
##
## dmsa_frame() wants the OTHER orientation: a samples x probes matrix whose
## rows line up with the phenotype data.frame and whose column names are bare
## probe ids (cg...). dmsa_import() is the adapter: it takes whichever end
## object a pipeline produced, extracts betas, fixes the orientation, aligns
## the sample sheet, and returns a `dmsa_import` object that dmsa_frame()
## accepts directly as its first argument.
##
## Design rules, in order of importance:
##   1. NEVER guess silently. Orientation, value scale and sample alignment are
##      each either unambiguous or an error that says what to pass.
##   2. No hard dependencies on any pipeline. SummarizedExperiment is consulted
##      through requireNamespace(); minfi, sesame and RnBeads are never
##      imported - their classes are recognised structurally.
##   3. Betas are the canonical scale. M-values are converted here (with a
##      note), because dmsa_frame() converts betas to M-values itself and
##      double conversion would corrupt the frame.
## ---------------------------------------------------------------------------

#' Import a preprocessing pipeline's result for DMSA
#'
#' Takes the end object of any mainstream methylation preprocessing pipeline
#' (minfi, sesame, meffil, ewastools, ChAMP, RnBeads, methylprep - see
#' Details) and returns a \code{dmsa_import} object holding a samples-by-probes
#' beta matrix and an aligned phenotype data.frame. Pass the result straight to
#' \code{\link{dmsa_frame}} as its first argument.
#'
#' @param x A pipeline result: a beta/M matrix in either orientation, a
#'   \code{data.frame} export (probes in the rownames or first column), a
#'   \code{SummarizedExperiment} (covers minfi's \code{GenomicRatioSet},
#'   \code{MethylSet} and relatives without needing minfi), a list of sesame
#'   \code{SigDF} objects, or a path to a \code{.csv}/\code{.rds} file.
#' @param pheno Optional phenotype/sample sheet: a \code{data.frame}, or a path
#'   to a CSV. When \code{x} is a \code{SummarizedExperiment} its
#'   \code{colData} is used automatically unless \code{pheno} overrides it.
#'   When no phenotype is available at all, a one-column sheet holding the
#'   sample names is built so the object still works, and the print method
#'   says so.
#' @param align_by Name of the \code{pheno} column holding the array/sample
#'   identifiers that match the methylation sample names (e.g.
#'   \code{"Basename"} or \code{"Sample_Name"}). \code{NULL} (default)
#'   auto-detects: each character/factor column is tried literally, then with
#'   directory paths stripped (minfi sample sheets), then the
#'   \code{Sentrix_ID}/\code{Sentrix_Position} pair pasted with "_". Ambiguity
#'   or failure is an error naming what was tried - never a silent guess.
#' @param values \code{"auto"} (default) detects the value scale: everything
#'   inside [0, 1] is treated as betas; anything outside as M-values, which
#'   are converted to betas via \code{2^m / (2^m + 1)} with a note.
#'   \code{"beta"} and \code{"m"} override the detection.
#' @param orientation \code{"auto"} (default) reads the dimension names: the
#'   axis whose names look like probe ids (\code{cg}/\code{ch.}/\code{rs}
#'   prefixes) is the probe axis. \code{"probes_rows"} and
#'   \code{"probes_cols"} override - required when a matrix has no dimnames
#'   on the probe axis.
#' @param id,format Optional study-shape declaration carried to the analysis
#'   step: the participant-id column in \code{pheno} and the declared
#'   \code{"wide"}/\code{"long"} format, exactly as in
#'   \code{\link{dmsa_design}}. Stored and printed; pass them on to
#'   \code{dmsa_design()}/\code{alpha_design()} when modelling.
#' @param source Optional label naming the pipeline the object came from
#'   (\code{"minfi"}, \code{"sesame"}, ...). Cosmetic: shown by the print
#'   method and recorded in the object; auto-filled when detectable.
#' @param ... Passed between methods.
#' @return An object of class \code{dmsa_import}: a list with \code{data} (the
#'   aligned phenotype data.frame), \code{methylation} (samples x probes beta
#'   matrix), \code{source}, \code{values} (scale detected), \code{align}
#'   (how samples were matched), \code{id}, \code{format} and \code{notes}.
#'   \code{dmsa_frame()} accepts it directly as its first argument.
#' @examples
#' ## a probes-x-samples beta matrix, the shape most pipelines produce
#' set.seed(1)
#' B <- matrix(runif(60), nrow = 6,
#'             dimnames = list(sprintf("cg%08d", 1:6), paste0("S", 1:10)))
#' ph <- data.frame(Sample_Name = paste0("S", 1:10), age = rnorm(10))
#' imp <- dmsa_import(B, pheno = ph)
#' imp
#' @export
dmsa_import <- function(x, ...) UseMethod("dmsa_import")

## ---- shared helpers --------------------------------------------------------

.imp_probe_like <- function(nm) {
  if (is.null(nm) || !length(nm)) return(0)
  mean(grepl("^(cg|ch\\.|rs|nv-|ctl_)", nm))
}

.imp_note <- function(obj, note) { obj$notes <- c(obj$notes, note); obj }

## Decide the value scale and return betas. M-values are converted; betas
## pass through. The 0/1-inclusive check tolerates the exact 0s and 1s that
## some pipelines emit at fully (un)methylated probes.
.imp_values <- function(M, values) {
  rng <- suppressWarnings(range(M, na.rm = TRUE))
  looks_beta <- is.finite(rng[1]) && rng[1] >= 0 && rng[2] <= 1
  scale <- switch(values,
    auto = if (looks_beta) "beta" else "m",
    beta = "beta",
    m    = "m")
  if (values == "beta" && !looks_beta)
    warning("values = \"beta\" declared, but values outside [0, 1] found ",
            "(range ", signif(rng[1], 3), " to ", signif(rng[2], 3),
            "); passing through unchanged", call. = FALSE)
  converted <- FALSE
  if (scale == "m") { M <- 2^M / (2^M + 1); converted <- TRUE }
  list(M = M, scale = scale, converted = converted)
}

## Fix orientation so samples are rows. `orientation` overrides; otherwise the
## probe axis is the one whose names look like probe ids.
.imp_orient <- function(M, orientation) {
  if (orientation == "probes_rows") return(list(M = t(M), how = "declared"))
  if (orientation == "probes_cols") return(list(M = M,    how = "declared"))
  pr <- .imp_probe_like(rownames(M))
  pc <- .imp_probe_like(colnames(M))
  if (pr >= 0.5 && pc < 0.5) return(list(M = t(M), how = "probe ids in rownames"))
  if (pc >= 0.5 && pr < 0.5) return(list(M = M,    how = "probe ids in colnames"))
  stop("cannot tell which axis holds the probes: ",
       if (pr >= 0.5 && pc >= 0.5) "BOTH dimension names look like probe ids"
       else "neither dimension's names look like probe ids (cg.../ch..../rs...)",
       ". Pass orientation = \"probes_rows\" (the usual pipeline shape, ",
       "samples in columns) or \"probes_cols\".", call. = FALSE)
}

## Align the phenotype sheet to the methylation rows. Returns the sheet
## reordered to the matrix's row order plus a description of how the match
## was made. Never guesses: no usable key is an error that lists what was
## tried, with the one exception of a same-size sheet with NO overlapping
## identifiers anywhere, which errors too - positional alignment must be
## asked for explicitly via align_by = ".position".
.imp_align <- function(M, pheno, align_by) {
  sm <- rownames(M)
  if (is.null(pheno)) {
    ph <- data.frame(.sample = if (is.null(sm)) paste0("S", seq_len(nrow(M)))
                               else sm, stringsAsFactors = FALSE)
    return(list(data = ph, how = "no pheno supplied; sample-name sheet built"))
  }
  if (is.character(pheno) && length(pheno) == 1L) {
    if (!file.exists(pheno)) stop("pheno file not found: ", pheno, call. = FALSE)
    pheno <- utils::read.csv(pheno, stringsAsFactors = FALSE, check.names = FALSE)
  }
  pheno <- as.data.frame(pheno)
  if (identical(align_by, ".position")) {
    if (nrow(pheno) != nrow(M))
      stop("align_by = \".position\" needs pheno and methylation to have the ",
           "same number of rows (", nrow(pheno), " vs ", nrow(M), ")",
           call. = FALSE)
    return(list(data = pheno, how = "positional (declared via align_by)"))
  }
  if (is.null(sm)) {
    if (nrow(pheno) == nrow(M))
      stop("the methylation matrix has no sample names, so rows cannot be ",
           "matched to `pheno` by identifier. If the two are already in the ",
           "same order, say so explicitly with align_by = \".position\".",
           call. = FALSE)
    stop("the methylation matrix has no sample names and a different number ",
         "of rows (", nrow(M), ") than `pheno` (", nrow(pheno), ")",
         call. = FALSE)
  }
  key_of <- function(v) as.character(v)
  candidates <- list()
  if (!is.null(align_by)) {
    if (!align_by %in% names(pheno))
      stop("align_by column '", align_by, "' is not in `pheno`", call. = FALSE)
    candidates[[align_by]] <- key_of(pheno[[align_by]])
  } else {
    for (cl in names(pheno)) {
      v <- pheno[[cl]]
      if (is.character(v) || is.factor(v)) candidates[[cl]] <- key_of(v)
    }
    ## minfi sample sheets: Basename holds full idat paths
    for (cl in names(candidates))
      candidates[[paste0(cl, " (basename)")]] <- basename(candidates[[cl]])
    ## Illumina sheets: barcode = Sentrix_ID _ Sentrix_Position
    sid <- intersect(c("Sentrix_ID", "Slide", "sentrix_id"), names(pheno))
    spos <- intersect(c("Sentrix_Position", "Array", "sentrix_position"),
                      names(pheno))
    if (length(sid) && length(spos))
      candidates[["Sentrix_ID_Sentrix_Position"]] <-
        paste(pheno[[sid[1]]], pheno[[spos[1]]], sep = "_")
  }
  hits <- vapply(candidates, function(k)
    !anyDuplicated(k) && all(sm %in% k), logical(1))
  if (!any(hits)) {
    stop("no `pheno` column matches the methylation sample names. Tried: ",
         paste(names(candidates), collapse = ", "),
         ". First sample names: ", paste(utils::head(sm, 3), collapse = ", "),
         ". Pass align_by = to name the identifier column, or ",
         "align_by = \".position\" if the rows are already in the same order.",
         call. = FALSE)
  }
  use <- names(candidates)[which(hits)[1]]
  ord <- match(sm, candidates[[use]])
  list(data = pheno[ord, , drop = FALSE],
       how = paste0("matched on pheno column '", use, "'"))
}

## Assemble the final object from a samples x probes matrix + alignment.
.imp_build <- function(M, pheno, align_by, values, source, id, format,
                       notes = character()) {
  vv <- .imp_values(M, values)
  al <- .imp_align(vv$M, pheno, align_by)
  rownames(al$data) <- NULL
  out <- list(data = al$data, methylation = vv$M, source = source,
              values = if (vv$converted)
                "M-values (converted to betas)" else "betas",
              align = al$how, id = id, format = format, notes = notes)
  if (vv$converted)
    out <- .imp_note(out, "M-values detected and converted to betas via 2^m/(2^m+1)")
  class(out) <- "dmsa_import"
  out
}

## ---- methods ---------------------------------------------------------------

#' @rdname dmsa_import
#' @export
dmsa_import.matrix <- function(x, pheno = NULL, align_by = NULL,
                               values = c("auto", "beta", "m"),
                               orientation = c("auto", "probes_rows",
                                               "probes_cols"),
                               id = NULL, format = NULL, source = "matrix",
                               ...) {
  values <- match.arg(values); orientation <- match.arg(orientation)
  if (!is.numeric(x)) stop("the matrix must be numeric", call. = FALSE)
  ori <- .imp_orient(x, orientation)
  .imp_build(ori$M, pheno, align_by, values, source, id, format,
             notes = paste0("orientation: ", ori$how))
}

#' @rdname dmsa_import
#' @export
dmsa_import.data.frame <- function(x, pheno = NULL, align_by = NULL,
                                   values = c("auto", "beta", "m"),
                                   orientation = c("auto", "probes_rows",
                                                   "probes_cols"),
                                   id = NULL, format = NULL,
                                   source = "data.frame", ...) {
  values <- match.arg(values); orientation <- match.arg(orientation)
  notes <- character()
  ## methylprep-style export: probe ids in the first column
  if (ncol(x) >= 2 && (is.character(x[[1]]) || is.factor(x[[1]])) &&
      .imp_probe_like(as.character(x[[1]])) >= 0.5) {
    rn <- as.character(x[[1]]); x <- x[, -1, drop = FALSE]; rownames(x) <- rn
    notes <- c(notes, "probe ids taken from the first column")
  }
  num <- vapply(x, is.numeric, logical(1))
  if (any(!num)) {
    notes <- c(notes, paste0(sum(!num), " non-numeric column(s) dropped: ",
                             paste(utils::head(names(x)[!num], 5),
                                   collapse = ", ")))
    x <- x[, num, drop = FALSE]
  }
  if (!ncol(x)) stop("no numeric methylation columns remain", call. = FALSE)
  M <- as.matrix(x)
  ori <- .imp_orient(M, orientation)
  .imp_build(ori$M, pheno, align_by, values, source, id, format,
             notes = c(notes, paste0("orientation: ", ori$how)))
}

#' @rdname dmsa_import
#' @export
dmsa_import.character <- function(x, ...) {
  if (length(x) != 1L)
    stop("pass one file path (or a matrix/data.frame/pipeline object)",
         call. = FALSE)
  if (!file.exists(x)) stop("file not found: ", x, call. = FALSE)
  ext <- tolower(tools::file_ext(x))
  obj <- switch(ext,
    csv = utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = utils::read.delim(x, stringsAsFactors = FALSE, check.names = FALSE),
    rds = readRDS(x),
    stop("unsupported file type '.", ext, "': use .csv, .tsv or .rds ",
         "(for methylprep, export betas to CSV first)", call. = FALSE))
  dmsa_import(obj, ...)
}

#' @rdname dmsa_import
#' @export
dmsa_import.list <- function(x, pheno = NULL, align_by = NULL,
                             id = NULL, format = NULL, ...) {
  ## a list of sesame SigDF objects (one per sample)
  is_sigdf <- length(x) > 0 && all(vapply(x, function(e)
    inherits(e, "SigDF") ||
      (is.data.frame(e) && all(c("Probe_ID", "MG", "UG") %in% names(e))),
    logical(1)))
  if (is_sigdf)
    stop("this is a list of sesame SigDF objects - per-sample signal, not ",
         "yet betas. Run openSesame(<idat dir>) (it returns the beta matrix ",
         "directly), or getBetas() on each SigDF and cbind() the results, ",
         "then dmsa_import() that matrix.", call. = FALSE)
  stop("dmsa_import() on a list expects sesame SigDF objects; for other ",
       "pipelines pass the beta matrix or the SummarizedExperiment.",
       call. = FALSE)
}

#' @rdname dmsa_import
#' @export
dmsa_import.default <- function(x, pheno = NULL, align_by = NULL,
                                values = c("auto", "beta", "m"),
                                id = NULL, format = NULL, source = NULL,
                                ...) {
  values <- match.arg(values)
  ## SummarizedExperiment covers minfi's whole preprocessed family
  ## (GenomicRatioSet, GenomicMethylSet, MethylSet, RatioSet) structurally -
  ## no minfi import needed, the assays tell us what we hold.
  if (requireNamespace("SummarizedExperiment", quietly = TRUE) &&
      methods::is(x, "SummarizedExperiment"))
    return(.imp_se(x, pheno, align_by, values, id, format, source))
  cls <- class(x)[1]
  if (cls %in% c("RGChannelSet", "RGChannelSetExtended"))
    stop("this is a raw red/green RGChannelSet - no methylation values have ",
         "been computed yet. Preprocess first (e.g. minfi::preprocessNoob(rgSet)",
         " then minfi::ratioConvert(), or any preprocess* route) and ",
         "dmsa_import() the result.", call. = FALSE)
  if (cls %in% c("RnBSet", "RnBeadSet", "RnBeadRawSet"))
    stop("this is an RnBeads object. Extract with RnBeads itself - ",
         "M <- meth(rnb.set, row.names = TRUE); ph <- pheno(rnb.set) - ",
         "then dmsa_import(M, pheno = ph).", call. = FALSE)
  if (cls == "SigDF")
    stop("this is a single sesame SigDF (one sample). Run ",
         "openSesame(<idat dir>) to get the beta matrix for all samples, or ",
         "pass a list of SigDF objects.", call. = FALSE)
  stop("dmsa_import() does not know class '", cls, "'. Supported: a beta/M ",
       "matrix, a data.frame export, a SummarizedExperiment (minfi results), ",
       "a list of sesame SigDF objects, or a .csv/.rds path.", call. = FALSE)
}

## SummarizedExperiment route: pick the right assay, compute betas, take
## colData as the sheet.
.imp_se <- function(x, pheno, align_by, values, id, format, source) {
  an <- SummarizedExperiment::assayNames(x)
  if (is.null(an)) an <- character()
  pick <- function(nm) SummarizedExperiment::assay(x, nm)
  notes <- character(); src <- source
  if (any(c("Green", "Red") %in% an))
    stop("this SummarizedExperiment holds raw Green/Red channels - ",
         "preprocess first (e.g. minfi::preprocessNoob) and dmsa_import() ",
         "the result.", call. = FALSE)
  if ("Beta" %in% an) {
    M <- pick("Beta"); values <- "beta"
    notes <- c(notes, "assay used: Beta")
    if (is.null(src)) src <- "SummarizedExperiment (Beta assay)"
  } else if ("M" %in% an) {
    M <- pick("M"); values <- "m"
    notes <- c(notes, "assay used: M (converted to betas)")
    if (is.null(src)) src <- "SummarizedExperiment (M assay)"
  } else if (all(c("Meth", "Unmeth") %in% an)) {
    M <- pick("Meth") / (pick("Meth") + pick("Unmeth") + 100)
    values <- "beta"
    notes <- c(notes, "betas computed as Meth/(Meth+Unmeth+100)")
    if (is.null(src)) src <- "SummarizedExperiment (Meth/Unmeth assays)"
  } else if (length(an)) {
    M <- pick(an[1])
    notes <- c(notes, paste0("assay used: ", an[1], " (first of ",
                             length(an), ")"))
    if (is.null(src)) src <- paste0("SummarizedExperiment (", an[1], " assay)")
  } else stop("the SummarizedExperiment carries no assays", call. = FALSE)
  M <- as.matrix(M)
  if (is.null(rownames(M)) && !is.null(rownames(x))) rownames(M) <- rownames(x)
  if (is.null(colnames(M)) && !is.null(colnames(x))) colnames(M) <- colnames(x)
  ph <- if (!is.null(pheno)) pheno else {
    cd <- as.data.frame(SummarizedExperiment::colData(x))
    if (ncol(cd)) { cd$.sample <- rownames(cd); cd } else NULL
  }
  ## probes are rows in every SummarizedExperiment methylation container
  .imp_build(t(M), ph, align_by, values, src, id, format, notes = notes)
}

## ---- print -----------------------------------------------------------------

#' @export
print.dmsa_import <- function(x, ...) {
  cat("dmsa_import:", nrow(x$methylation), "samples x",
      ncol(x$methylation), "probes\n")
  cat("  source:   ", x$source, "\n")
  cat("  values:   ", x$values, "\n")
  cat("  alignment:", x$align, "\n")
  if (!is.null(x$id) || !is.null(x$format))
    cat("  declared: ",
        if (!is.null(x$format)) paste0("format = ", x$format) else "",
        if (!is.null(x$id)) paste0(" id = ", x$id) else "", "\n", sep = "")
  cov <- tryCatch({
    mp <- .frame_read_map("alpha")
    sum(unique(mp$probe) %in% colnames(x$methylation))
  }, error = function(e) NULL)
  if (!is.null(cov))
    cat("  alpha map:", cov, "of", ncol(x$methylation),
        "probes are Alpha-mapped\n")
  for (n in x$notes) cat("  note:     ", n, "\n")
  cat("  next:      dmsa_frame(<this object>, outcome = ..., ...)\n")
  invisible(x)
}
