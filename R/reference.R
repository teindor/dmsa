# dmsa reference bundles
# --------------------------------------------------------------------------
# Everything above this file is Alpha-specific: 549 genes, 30 curated systems,
# a hand-made polarity table. For DMSA to run on any panel and any organism-wide
# gene set, three things have to come from outside:
#
#   systems   gene -> system, and gene -> module inside that system
#   polarity  gene -> w_g within each system (the activation/brake sign)
#   anchors   which genes DEFINE "more activation" for each system
#
# A reference bundle carries all three. `alpha_reference()` wraps the bundled
# Alpha tables in the same shape, so Alpha becomes one bundle among others
# rather than a special case.
#
# The anchor set is the part that cannot be automated away, and it is worth
# being explicit about why. "System activation tone" is only defined once you say
# what the system's output is. Alpha has that because the PI curated drivers.
# For an arbitrary Reactome pathway or GO term it has to be supplied, derived
# from graph topology, or the system has to be treated as unsigned. A bundle
# records which of those was done, in `anchor_method`, so a reader can tell.

#' Construct a DMSA reference bundle
#'
#' @param systems data.frame with columns \code{gene}, \code{system_id},
#'   \code{system}, and optionally \code{module_id} and \code{module}. One row
#'   per gene-within-system (a gene may appear in several systems).
#' @param polarity data.frame with columns \code{gene}, \code{system_id},
#'   \code{w_g}. \code{w_g} may be \code{-1/0/+1} or continuous in
#'   \code{[-1, 1]}. Optional columns \code{evidence}, \code{role},
#'   \code{confidence}, \code{source} are carried through.
#' @param anchors Optional data.frame with \code{system_id}, \code{gene} - the
#'   genes that define the direction of "more activation" for each system.
#' @param anchor_method Character, how the anchors were obtained. One of
#'   \code{"curated"}, \code{"graph_sink"}, \code{"user"}, or \code{"none"}
#'   (systems then carry no sign and only gene-level analysis is valid).
#' @param name,version,notes Provenance, printed and stored.
#' @return An object of class \code{dmsa_reference}.
#' @export
dmsa_reference <- function(systems, polarity = NULL, anchors = NULL,
                           anchor_method = c("curated", "graph_sink", "user", "none"),
                           name = "unnamed", version = NA_character_,
                           notes = NULL) {
  anchor_method <- match.arg(anchor_method)
  systems <- as.data.frame(systems, stringsAsFactors = FALSE)
  req <- c("gene", "system_id", "system")
  miss <- setdiff(req, names(systems))
  if (length(miss)) stop("systems is missing: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  systems$gene <- as.character(systems$gene)
  systems$system_id <- as.character(systems$system_id)
  if (!"module_id" %in% names(systems)) {
    systems$module_id <- NA_character_; systems$module <- NA_character_
  } else systems$module_id <- as.character(systems$module_id)
  if (anyDuplicated(systems[c("gene", "system_id", "module_id")]))
    systems <- unique(systems)

  if (!is.null(polarity)) {
    polarity <- as.data.frame(polarity, stringsAsFactors = FALSE)
    miss <- setdiff(c("gene", "system_id", "w_g"), names(polarity))
    if (length(miss)) stop("polarity is missing: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    polarity$gene <- as.character(polarity$gene)
    polarity$system_id <- as.character(polarity$system_id)
    polarity$w_g <- as.numeric(polarity$w_g)
    if (any(abs(polarity$w_g) > 1 + 1e-9, na.rm = TRUE))
      stop("w_g must lie in [-1, 1]", call. = FALSE)
    if (anyDuplicated(polarity[c("gene", "system_id")]))
      stop("polarity has more than one w_g for the same gene and system",
           call. = FALSE)
    unknown <- setdiff(paste(polarity$gene, polarity$system_id),
                       paste(systems$gene, systems$system_id))
    if (length(unknown))
      warning(length(unknown), " polarity row(s) refer to gene-system pairs ",
              "absent from `systems`; they will never be used")
  }
  if (!is.null(anchors)) {
    anchors <- as.data.frame(anchors, stringsAsFactors = FALSE)
    miss <- setdiff(c("gene", "system_id"), names(anchors))
    if (length(miss)) stop("anchors is missing: ", paste(miss, collapse = ", "),
                           call. = FALSE)
    anchors$system_id <- as.character(anchors$system_id)
  }
  if (anchor_method == "none" && !is.null(polarity) &&
      any(polarity$w_g != 0, na.rm = TRUE))
    warning("anchor_method = 'none' but polarity carries non-zero w_g; ",
            "system-level signs are then undefined in provenance terms")

  structure(list(name = name, version = version, notes = notes,
                 anchor_method = anchor_method, systems = systems,
                 polarity = polarity, anchors = anchors),
            class = "dmsa_reference")
}

#' @export
print.dmsa_reference <- function(x, ...) {
  cat("dmsa reference: ", x$name,
      if (!is.na(x$version)) paste0(" (", x$version, ")"), "\n", sep = "")
  s <- x$systems
  cat(sprintf("  %d genes across %d systems", length(unique(s$gene)),
              length(unique(s$system_id))))
  nm <- sum(!is.na(s$module_id))
  cat(if (nm) sprintf(", %d modules\n", length(unique(stats::na.omit(s$module_id))))
      else "  (no module layer)\n")
  if (is.null(x$polarity)) cat("  polarity: none - gene-level analysis only\n")
  else {
    p <- x$polarity
    cont <- any(!p$w_g %in% c(-1, 0, 1), na.rm = TRUE)
    cat(sprintf("  polarity: %d gene-system pairs, %s (%d activating, %d braking, %d off-axis)\n",
                nrow(p), if (cont) "continuous" else "signed +-1",
                sum(p$w_g > 0, na.rm = TRUE), sum(p$w_g < 0, na.rm = TRUE),
                sum(p$w_g == 0, na.rm = TRUE)))
  }
  cat("  anchors: ", x$anchor_method,
      if (!is.null(x$anchors)) sprintf(" (%d genes)", nrow(x$anchors)) else "",
      "\n", sep = "")
  if (!is.null(x$notes)) cat("  note: ", x$notes, "\n", sep = "")
  invisible(x)
}

#' The bundled Project Alpha tables as a reference bundle
#'
#' @param modules Optional data.frame with \code{system_id}, \code{gene},
#'   \code{module_id}, \code{module} to add a module layer (e.g. the HPA
#'   H / P / A / negative-feedback split).
#' @return A \code{dmsa_reference}.
#' @export
alpha_reference <- function(modules = NULL) {
  sysd <- alpha_gene_systems()
  pol <- alpha_polarity()
  s <- data.frame(gene = sysd$gene, system_id = as.character(sysd$system_id),
                  system = sysd$system, stringsAsFactors = FALSE)
  if (!is.null(modules)) {
    modules <- as.data.frame(modules, stringsAsFactors = FALSE)
    modules$system_id <- as.character(modules$system_id)
    s <- merge(s, modules[c("gene", "system_id", "module_id", "module")],
               by = c("gene", "system_id"), all.x = TRUE)
  }
  p <- data.frame(gene = pol$gene, system_id = as.character(pol$system_id),
                  w_g = pol$w_g, role = pol$role, confidence = pol$confidence,
                  stringsAsFactors = FALSE)
  a <- p[p$w_g == 1 & !is.na(p$role) & p$role == "driver", c("system_id", "gene")]
  dmsa_reference(s, p, a, anchor_method = "curated", name = "Project Alpha 2026",
                 version = "draft polarity",
                 notes = "polarity is a DRAFT pending per-gene citations")
}

#' Build a cascade tree for a set of probes from a reference bundle
#'
#' Returns the \code{tree} argument \code{dmsa_cascade()} expects: one row per
#' probe, columns ordered outermost first.
#'
#' @param reference A \code{dmsa_reference}.
#' @param genes Character vector, the set-membership gene of each probe.
#' @param system_id Optional: restrict to one system (a gene can sit in several,
#'   which would otherwise duplicate probes).
#' @param levels Which levels to include, in order. Defaults to
#'   \code{c("system", "module")} when a module layer exists, else
#'   \code{"system"}. \code{"gene"} is always appended.
#' @return data.frame with one row per probe. Probes whose gene is not in the
#'   reference get \code{NA} and should be routed to the flat unannotated arm.
#' @export
dmsa_tree <- function(reference, genes, system_id = NULL, levels = NULL) {
  stopifnot(inherits(reference, "dmsa_reference"))
  s <- reference$systems
  if (!is.null(system_id)) s <- s[s$system_id == as.character(system_id), , drop = FALSE]
  has_mod <- any(!is.na(s$module_id))
  if (is.null(levels)) levels <- if (has_mod) c("system", "module") else "system"
  genes <- as.character(genes)
  key <- s[!duplicated(s$gene), , drop = FALSE]
  if (is.null(system_id) && anyDuplicated(s$gene))
    warning("some genes belong to several systems; supply system_id to avoid ",
            "silently keeping only the first")
  i <- match(genes, key$gene)
  out <- data.frame(row.names = seq_along(genes))
  if ("system" %in% levels) out$system <- key$system_id[i]
  if ("module" %in% levels) out$module <- key$module_id[i]
  out$gene <- genes
  out
}

#' Polarity table for one system, in the form dmsa_align() wants
#' @param reference A \code{dmsa_reference}.
#' @param system_id The system to extract.
#' @return data.frame with \code{gene} and \code{w_g}, or NULL if the bundle
#'   carries no polarity.
#' @export
dmsa_polarity_for <- function(reference, system_id) {
  stopifnot(inherits(reference, "dmsa_reference"))
  if (is.null(reference$polarity)) return(NULL)
  p <- reference$polarity
  p <- p[p$system_id == as.character(system_id), c("gene", "w_g"), drop = FALSE]
  if (!nrow(p)) return(NULL)
  p
}

#' Read a reference bundle from CSV files
#'
#' The format written by the bundle builder: \code{systems.csv},
#' \code{polarity.csv} and optionally \code{anchors.csv} in one directory,
#' plus an optional \code{manifest.dcf} carrying name, version and notes.
#'
#' @param dir Directory holding the files.
#' @param anchor_method Passed to \code{dmsa_reference()} when no manifest is
#'   present.
#' @return A \code{dmsa_reference}.
#' @export
dmsa_reference_read <- function(dir, anchor_method = "user") {
  f <- function(x) file.path(dir, x)
  if (!file.exists(f("systems.csv")))
    stop("no systems.csv in ", dir, call. = FALSE)
  s <- utils::read.csv(f("systems.csv"), stringsAsFactors = FALSE)
  p <- if (file.exists(f("polarity.csv")))
    utils::read.csv(f("polarity.csv"), stringsAsFactors = FALSE) else NULL
  a <- if (file.exists(f("anchors.csv")))
    utils::read.csv(f("anchors.csv"), stringsAsFactors = FALSE) else NULL
  meta <- list(name = basename(dir), version = NA_character_, notes = NULL,
               anchor_method = anchor_method)
  if (file.exists(f("manifest.dcf"))) {
    d <- read.dcf(f("manifest.dcf"))
    for (k in intersect(colnames(d), names(meta))) meta[[k]] <- as.character(d[1, k])
  }
  dmsa_reference(s, p, a, anchor_method = meta$anchor_method,
                 name = meta$name, version = meta$version, notes = meta$notes)
}
