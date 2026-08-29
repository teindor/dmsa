# Single-file reference input
# --------------------------------------------------------------------------
# dmsa_reference_read() wants a directory of four files, which is fine for a
# generated bundle and a poor way to ask a person for their own pathway. This
# accepts ONE flat csv in system / module / gene form, so a user whose system is
# nowhere in Reactome or the Alpha panel can still run every level of DMSA.

.norm_names <- function(x) tolower(gsub("[^a-z0-9]+", "", tolower(x)))

.pick_col <- function(df, wanted, required = FALSE, what = wanted[1]) {
  nn <- .norm_names(names(df))
  hit <- which(nn %in% .norm_names(wanted))[1]
  if (is.na(hit)) {
    if (required)
      stop("no column found for ", what, ". Accepted names: ",
           paste(wanted, collapse = ", "),
           ".\nColumns present: ", paste(names(df), collapse = ", "),
           call. = FALSE)
    return(NULL)
  }
  names(df)[hit]
}

#' Read a user-supplied system / module / gene csv as a reference bundle
#'
#' One flat file, one row per gene-within-module (or per gene-within-system if
#' there is no module layer). Column names are matched case- and
#' punctuation-insensitively, so \code{System}, \code{system_id} and
#' \code{"System ID"} are all accepted.
#'
#' \strong{Required:} a system column (\code{system}, \code{system_id},
#' \code{pathway}, \code{set}) and a gene column (\code{gene}, \code{symbol},
#' \code{gene_symbol}).
#'
#' \strong{Optional:}
#' \itemize{
#'   \item module (\code{module}, \code{subsystem}, \code{submodule},
#'     \code{module_id}) - adds the module level to the cascade.
#'   \item \code{w_g} (\code{wg}, \code{weight}, \code{polarity}, \code{sign}) -
#'     the gene's effect on system activation, \code{-1}/\code{0}/\code{+1} or
#'     continuous in \code{[-1, 1]}. Without it the bundle carries no polarity
#'     and only gene-level DMSA is valid.
#'   \item anchor (\code{anchor}, \code{is_anchor}, \code{driver}) - TRUE for the
#'     genes that define what "more activation" means. Without it, anchors are
#'     taken from \code{role == "driver"} if a role column exists, else inferred
#'     from the most positive \code{w_g} with a warning.
#'   \item role, confidence, source - carried through untouched.
#' }
#'
#' @param path Path to the csv (anything \code{utils::read.csv} accepts).
#' @param name,version,notes Provenance for the bundle; \code{name} defaults to
#'   the file name.
#' @param min_genes Systems with fewer genes than this are dropped with a
#'   message - a set of two genes is not a set test. Set to 1 to keep everything.
#' @param quiet Suppress the summary message.
#' @return A \code{dmsa_reference}.
#' @examples
#' f <- tempfile(fileext = ".csv")
#' dmsa_reference_template(f)
#' ref <- dmsa_reference_csv(f)
#' @export
dmsa_reference_csv <- function(path, name = NULL, version = NA_character_,
                              notes = NULL, min_genes = 3L, quiet = FALSE) {
  if (!file.exists(path)) stop("no such file: ", path, call. = FALSE)
  ## peek at the header so any id-like column is read as character - numeric
  ## parsing merges module ids "2.10" and "2.1" (same rule as the cascade
  ## reader). Names are matched loosely because this reader accepts aliases.
  .hd <- names(utils::read.csv(path, nrows = 1, check.names = FALSE))
  .idn <- .hd[grepl("(^|_)(system|module|set)_?id$|^systemid$|^moduleid$",
                    tolower(.hd))]
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                        colClasses =
                          if (length(.idn))
                            stats::setNames(rep("character", length(.idn)),
                                            .idn) else NA)
  if (!nrow(df)) stop("the file has no rows", call. = FALSE)

  c_sys  <- .pick_col(df, c("system", "system_id", "systemid", "pathway", "set",
                            "set_id"), TRUE, "the system")
  c_gene <- .pick_col(df, c("gene", "genes", "symbol", "gene_symbol",
                            "genesymbol", "hgnc"), TRUE, "the gene")
  c_mod  <- .pick_col(df, c("module", "module_id", "moduleid", "subsystem",
                            "submodule", "subset", "component"))
  c_w    <- .pick_col(df, c("w_g", "wg", "weight", "polarity", "sign",
                            "direction_on_system"))
  c_anc  <- .pick_col(df, c("anchor", "is_anchor", "isanchor", "driver",
                            "is_driver"))
  c_role <- .pick_col(df, c("role"))
  c_conf <- .pick_col(df, c("confidence", "conf"))
  c_src  <- .pick_col(df, c("source", "reference", "citation", "pmid"))

  gene <- trimws(as.character(df[[c_gene]]))
  sys  <- trimws(as.character(df[[c_sys]]))
  keep <- nzchar(gene) & nzchar(sys) & !is.na(gene) & !is.na(sys)
  if (any(!keep)) {
    message("dropping ", sum(!keep), " row(s) with a blank system or gene")
    df <- df[keep, , drop = FALSE]; gene <- gene[keep]; sys <- sys[keep]
  }

  systems <- data.frame(gene = gene, system_id = sys, system = sys,
                        stringsAsFactors = FALSE)
  if (!is.null(c_mod)) {
    md <- trimws(as.character(df[[c_mod]]))
    md[!nzchar(md) | is.na(md)] <- NA_character_
    ## a module id must be unique across systems, so qualify it
    systems$module_id <- ifelse(is.na(md), NA_character_, paste(sys, md, sep = "::"))
    systems$module <- md
  }
  systems <- unique(systems)

  ## drop unusably small systems
  n_by <- table(unique(systems[c("gene", "system_id")])$system_id)
  small <- names(n_by)[n_by < min_genes]
  if (length(small)) {
    if (!quiet) message("dropping ", length(small), " system(s) with fewer than ",
                        min_genes, " genes: ",
                        paste(utils::head(small, 5), collapse = ", "),
                        if (length(small) > 5) ", ..." else "")
    systems <- systems[!systems$system_id %in% small, , drop = FALSE]
  }
  if (!nrow(systems)) stop("no system survived the min_genes filter", call. = FALSE)

  ## polarity
  polarity <- NULL
  if (!is.null(c_w)) {
    w <- suppressWarnings(as.numeric(df[[c_w]]))
    bad <- which(!is.na(df[[c_w]]) & nzchar(as.character(df[[c_w]])) & is.na(w))
    if (length(bad))
      stop("w_g column '", c_w, "' has ", length(bad),
           " value(s) that are not numeric, first at row ", bad[1],
           ": '", df[[c_w]][bad[1]], "'", call. = FALSE)
    if (any(abs(w) > 1, na.rm = TRUE)) {
      rng <- max(abs(w), na.rm = TRUE)
      stop("w_g must lie in [-1, 1]; largest absolute value is ", rng,
           ". Rescale, or use -1/0/+1.", call. = FALSE)
    }
    polarity <- data.frame(gene = gene, system_id = sys, w_g = w,
                           stringsAsFactors = FALSE)
    if (!is.null(c_role)) polarity$role <- as.character(df[[c_role]])
    if (!is.null(c_conf)) polarity$confidence <- as.character(df[[c_conf]])
    if (!is.null(c_src))  polarity$source <- as.character(df[[c_src]])
    polarity <- polarity[!is.na(polarity$w_g), , drop = FALSE]
    ## duplicates with CONFLICTING signs must not resolve silently first-wins
    .k <- paste(polarity$gene, polarity$system_id)
    if (anyDuplicated(.k)) {
      .con <- vapply(split(polarity$w_g, .k),
                     function(z) length(unique(z)) > 1L, TRUE)
      if (any(.con))
        warning("conflicting w_g for the same gene within one system (",
                paste(utils::head(names(.con)[.con], 3), collapse = "; "),
                "); keeping the first value of each. Resolve the duplicates ",
                "in the CSV.", call. = FALSE)
    }
    polarity <- polarity[!duplicated(polarity[c("gene", "system_id")]), , drop = FALSE]
    polarity <- polarity[paste(polarity$gene, polarity$system_id) %in%
                           paste(systems$gene, systems$system_id), , drop = FALSE]
    if (!nrow(polarity)) polarity <- NULL
  }

  ## anchors
  anchors <- NULL; amethod <- "none"
  if (!is.null(polarity)) {
    if (!is.null(c_anc)) {
      a <- df[[c_anc]]
      flag <- if (is.logical(a)) a else
        tolower(trimws(as.character(a))) %in% c("1", "true", "t", "yes", "y")
      anchors <- unique(data.frame(system_id = sys[flag], gene = gene[flag],
                                   stringsAsFactors = FALSE))
      amethod <- "user"
    } else if (!is.null(c_role) &&
               any(tolower(trimws(as.character(df[[c_role]]))) == "driver",
                   na.rm = TRUE)) {
      ## EXACT role match. The role vocabulary contains "driver-adjacent" and
      ## "brake-of-driver"; a substring match promoted a BRAKE to the anchor
      ## that defines what activation means for its system.
      dr <- !is.na(df[[c_role]]) &
        tolower(trimws(as.character(df[[c_role]]))) == "driver"
      anchors <- unique(data.frame(system_id = sys[dr], gene = gene[dr],
                                   stringsAsFactors = FALSE))
      amethod <- "user"   # derived from the user's role column, not curation
    } else {
      top <- do.call(rbind, lapply(split(polarity, polarity$system_id), function(p) {
        if (all(p$w_g <= 0)) return(NULL)
        p[p$w_g == max(p$w_g), c("system_id", "gene"), drop = FALSE]
      }))
      if (!is.null(top) && nrow(top)) {
        anchors <- unique(top); amethod <- "user"
        warning("no anchor or role column: anchors inferred from the most ",
                "positive w_g in each system. Supply an 'anchor' column to be ",
                "explicit about what defines activation.", call. = FALSE)
      }
    }
    ## anchors were built from the pre-filter row vectors: keep only anchors
    ## whose system survived min_genes (polarity above got the same treatment)
    if (!is.null(anchors))
      anchors <- anchors[anchors$system_id %in% systems$system_id, , drop = FALSE]
    if (!is.null(anchors) && !nrow(anchors)) { anchors <- NULL; amethod <- "none" }
  }
  if (is.null(polarity) && !quiet)
    message("no w_g column: the bundle carries no polarity, so system-level ",
            "signs are undefined and only GENE-level DMSA is valid.")

  ref <- dmsa_reference(systems, polarity, anchors, anchor_method = amethod,
                        name = if (is.null(name)) basename(path) else name,
                        version = version, notes = notes)
  if (!quiet) print(ref)
  ref
}

#' Write a template csv for a user-supplied reference
#'
#' Produces a small, filled-in example with every accepted column, so the shape
#' is unambiguous. Overwrite the rows with your own.
#'
#' @param path Where to write it.
#' @param with_modules Include the module column.
#' @return \code{path}, invisibly.
#' @examples
#' f <- tempfile(fileext = ".csv")
#' dmsa_reference_template(f)
#' head(utils::read.csv(f), 3)
#' ## the filled-in rows show every accepted column; replace them with your own
#' ## panel, then read the result back as a reference bundle
#' dmsa_reference_csv(f, quiet = TRUE)
#' @export
dmsa_reference_template <- function(path, with_modules = TRUE) {
  d <- data.frame(
    system = c(rep("HPA axis", 5), rep("Oxytocin", 3)),
    module = c("hypothalamic", "hypothalamic", "pituitary", "adrenal",
               "negative feedback", "ligand", "receptor", "receptor"),
    gene   = c("CRH", "CRHBP", "POMC", "MC2R", "NR3C1",
               "OXT", "OXTR", "AVPR1A"),
    w_g    = c(1, -1, 1, 1, -1, 1, 1, 0.5),
    anchor = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE),
    role   = c("driver", "brake", "driver", "driver", "feedback",
               "driver", "transducer", "transducer"),
    confidence = c("high", "high", "high", "medium", "high",
                   "high", "high", "low"),
    source = c("PMID:...", "PMID:...", "PMID:...", "", "PMID:...",
               "PMID:...", "PMID:...", ""),
    stringsAsFactors = FALSE)
  if (!with_modules) d$module <- NULL
  utils::write.csv(d, path, row.names = FALSE)
  message("template written to ", path,
          "\n  system + gene are required; module, w_g, anchor, role, ",
          "confidence, source are optional.",
          "\n  w_g: +1 raises system activation, -1 opposes it, 0 off-axis; ",
          "continuous values in [-1, 1] are allowed.",
          "\n  anchor: TRUE for the genes that DEFINE activation for that system.")
  invisible(path)
}

#' Write a reference bundle to a directory
#'
#' The inverse of \code{dmsa_reference_read()}: useful for saving a bundle built
#' from a one-file csv, or a generated bundle you have edited.
#'
#' @param reference A \code{dmsa_reference}.
#' @param dir Directory to create and write into.
#' @return \code{dir}, invisibly.
#' @examples
#' f <- tempfile(fileext = ".csv")
#' dmsa_reference_template(f)
#' ref <- dmsa_reference_csv(f, quiet = TRUE)
#'
#' d <- tempfile()
#' dmsa_reference_write(ref, d)
#' list.files(d)
#' nrow(dmsa_reference_read(d)$systems)
#' @export
dmsa_reference_write <- function(reference, dir) {
  stopifnot(inherits(reference, "dmsa_reference"))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(reference$systems, file.path(dir, "systems.csv"),
                   row.names = FALSE)
  if (!is.null(reference$polarity))
    utils::write.csv(reference$polarity, file.path(dir, "polarity.csv"),
                     row.names = FALSE)
  if (!is.null(reference$anchors))
    utils::write.csv(reference$anchors, file.path(dir, "anchors.csv"),
                     row.names = FALSE)
  writeLines(c(paste0("name: ", reference$name),
               paste0("version: ", reference$version),
               paste0("anchor_method: ", reference$anchor_method),
               paste0("notes: ", if (is.null(reference$notes)) ""
                      else reference$notes)),
             file.path(dir, "manifest.dcf"))
  invisible(dir)
}
