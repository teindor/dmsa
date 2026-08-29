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
#' @examples
#' ## One row per gene-within-system, a sign per gene, and the anchor genes that
#' ## fix which way "more activation" points for each system.
#' sys <- data.frame(gene = c("CRH", "POMC", "NR3C1", "OXT", "OXTR"),
#'                   system_id = c("1", "1", "1", "2", "2"),
#'                   system = c(rep("HPA axis", 3), rep("Oxytocin", 2)))
#' pol <- data.frame(gene = sys$gene, system_id = sys$system_id,
#'                   w_g = c(1, 1, -1, 1, 0.6))
#' anc <- data.frame(system_id = c("1", "2"), gene = c("CRH", "OXT"))
#' dmsa_reference(sys, pol, anc, anchor_method = "curated", name = "toy",
#'                version = "0.1")
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

#' The bundled Project Alpha biological reference
#'
#' Reads the audited 2026c gene codebook that ships with the package and
#' returns it as a \code{\link{dmsa_reference}}: one row per gene within its
#' system and module, with the curated gene-to-system polarity attached.
#'
#' The reference is BIOLOGICAL - it defines \code{system > module > gene} and
#' stops there. Which CpGs exist is a property of the methylation data a user
#' supplies, not of the biology, so a gene with no probe in any particular
#' cohort still belongs to its system here. Genes with no curated polarity are
#' left UNRESOLVED; they are never silently treated as \code{+1}.
#'
#' @param modules Optional data.frame with \code{system_id}, \code{gene},
#'   \code{module_id}, \code{module} to add a module layer (e.g. the HPA
#'   H / P / A / negative-feedback split).
#' @return A \code{dmsa_reference}.
#' @examples
#' ref <- alpha_reference()
#' ref
#' # the polarity dmsa_align() consumes for one system (2 = HPA axis)
#' head(dmsa_polarity_for(ref, "2"), 4)
#' @export
alpha_reference <- function(modules = NULL) {
  ## spec 2 + 5: the reference is BIOLOGICAL - system > module > gene, and
  ## nothing below gene. It is read from the audited 2026c gene codebook, NOT
  ## from alpha_gene_systems(), which is derived from retained probe coverage
  ## (n_probes / n_kept / usability) and therefore drops any gene this cohort
  ## happened to have no usable probe for. Building biology from one array's
  ## coverage silently redefines the biology; the codebook keeps all genes,
  ## including those with zero retained Alpha probes, so another lab supplying
  ## a CpG for such a gene can analyse it.
  p <- system.file("extdata", "alpha_reference_2026c.csv.gz", package = "dmsa")
  if (nzchar(p)) {
    s <- utils::read.csv(gzfile(p), stringsAsFactors = FALSE)
    for (v in c("system_id", "module_id", "gene")) s[[v]] <- as.character(s[[v]])
    s <- s[, c("gene", "system_id", "system", "module_id", "module")]
  } else {
    ## fallback for an install without the codebook resource
    sysd <- alpha_gene_systems()
    s <- data.frame(gene = sysd$gene, system_id = as.character(sysd$system_id),
                    system = sysd$system, stringsAsFactors = FALSE)
  }
  pol <- alpha_polarity()
  if (!is.null(modules)) {
    ## a user-supplied module layer OVERRIDES the codebook's for the genes it
    ## names, so an alternative partition (e.g. an H/P/A/negative-feedback
    ## split) can still be declared
    modules <- as.data.frame(modules, stringsAsFactors = FALSE)
    modules$system_id <- as.character(modules$system_id)
    modules$module_id <- as.character(modules$module_id)
    ## override PER GENE: genes the user table does not name keep the
    ## codebook's module. Deleting the whole layer first silently stripped
    ## ~1,250 unnamed genes of their modules whenever a user re-partitioned
    ## one system.
    if (anyDuplicated(modules[c("gene", "system_id")]))
      stop("`modules` lists the same gene twice within one system; a gene ",
           "belongs to exactly one module of a system.", call. = FALSE)
    i <- match(paste(s$gene, s$system_id),
               paste(modules$gene, modules$system_id))
    hit <- !is.na(i)
    s$module_id[hit] <- modules$module_id[i[hit]]
    s$module[hit] <- modules$module[i[hit]]
  }
  p <- data.frame(gene = pol$gene, system_id = as.character(pol$system_id),
                  w_g = pol$w_g, role = pol$role, confidence = pol$confidence,
                  stringsAsFactors = FALSE)
  a <- p[!is.na(p$w_g) & p$w_g == 1 & !is.na(p$role) & p$role == "driver",
         c("system_id", "gene")]
  dmsa_reference(s, p, a, anchor_method = "curated",
                 name = "Project Alpha (2026c, biological)", version = "2026c",
                 notes = paste("biological reference: system > module > gene.",
                               "Genes with no curated polarity are UNRESOLVED",
                               "(w_g absent), never treated as +1.",
                               "Polarity remains a DRAFT pending per-gene",
                               "citations."))
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
#' @examples
#' systems <- data.frame(
#'   gene = c("CRH", "NR3C1", "FKBP5", "OXT"),
#'   system_id = c("HPA", "HPA", "HPA", "OXT"),
#'   system = c("HPA axis", "HPA axis", "HPA axis", "Oxytocin"),
#'   module_id = c("HPA.a", "HPA.b", "HPA.b", "OXT.a"),
#'   module = c("Drive", "Receptors", "Receptors", "Ligand"))
#' ref <- dmsa_reference(systems, anchor_method = "none")
#'
#' ## a probe whose gene is not in the bundle gets NA and goes to the flat arm
#' dmsa_tree(ref, genes = c("CRH", "NR3C1", "FKBP5", "NOTINSET"),
#'           system_id = "HPA")
#' @export
dmsa_tree <- function(reference, genes, system_id = NULL, levels = NULL) {
  stopifnot(inherits(reference, "dmsa_reference"))
  s <- reference$systems
  if (!is.null(system_id)) s <- s[s$system_id == as.character(system_id), , drop = FALSE]
  has_mod <- any(!is.na(s$module_id))
  if (is.null(levels)) levels <- if (has_mod) c("system", "module") else "system"
  genes <- as.character(genes)
  key <- s[!duplicated(s$gene), , drop = FALSE]
  if (anyDuplicated(s$gene))
    warning(if (is.null(system_id))
              "some genes belong to several systems; supply system_id to avoid silently keeping only the first"
            else
              "some genes appear in several modules of this system; keeping each gene's first module")
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
#' @examples
#' f <- tempfile(fileext = ".csv")
#' dmsa_reference_template(f)
#' ref <- dmsa_reference_csv(f, quiet = TRUE)
#' ## w_g is the gene's sign toward its system's activation tone: +1 raises it,
#' ## -1 opposes it. dmsa_align() multiplies w_g by each CpG's
#' ## methylation-to-expression direction to get the probe's aligned sign.
#' dmsa_polarity_for(ref, "HPA axis")
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
#' @examples
#' ## A bundle on disk is systems.csv, polarity.csv, anchors.csv and an optional
#' ## manifest.dcf carrying the provenance, all in one directory.
#' sys <- data.frame(gene = c("CRH", "POMC", "NR3C1"), system_id = "1",
#'                   system = "HPA axis")
#' pol <- data.frame(gene = sys$gene, system_id = "1", w_g = c(1, 1, -1))
#' d <- file.path(tempdir(), "hpa_bundle")
#' dmsa_reference_write(dmsa_reference(sys, pol, anchor_method = "curated",
#'                                     name = "hpa-toy", version = "1"), d)
#' list.files(d)
#' dmsa_reference_read(d)
#' @export
dmsa_reference_read <- function(dir, anchor_method = "user") {
  f <- function(x) file.path(dir, x)
  if (!file.exists(f("systems.csv")))
    stop("no systems.csv in ", dir, call. = FALSE)
  ## id columns are read as character ALWAYS: parsed as numeric, module ids
  ## "2.10" and "2.1" become the same module (the cascade reader documents the
  ## same rule; these readers must obey it too).
  .rd <- function(fp) {
    hd <- names(utils::read.csv(fp, nrows = 1))
    idc <- intersect(c("system_id", "module_id"), hd)
    utils::read.csv(fp, stringsAsFactors = FALSE,
                    colClasses = if (length(idc))
                      stats::setNames(rep("character", length(idc)), idc)
                    else NA)
  }
  s <- .rd(f("systems.csv"))
  p <- if (file.exists(f("polarity.csv"))) .rd(f("polarity.csv")) else NULL
  a <- if (file.exists(f("anchors.csv"))) .rd(f("anchors.csv")) else NULL
  meta <- list(name = basename(dir), version = NA_character_, notes = NULL,
               anchor_method = anchor_method)
  if (file.exists(f("manifest.dcf"))) {
    d <- read.dcf(f("manifest.dcf"))
    for (k in intersect(colnames(d), names(meta))) meta[[k]] <- as.character(d[1, k])
  }
  dmsa_reference(s, p, a, anchor_method = meta$anchor_method,
                 name = meta$name, version = meta$version, notes = meta$notes)
}

## ---------------------------------------------------------------------------
## THE BUNDLED ALPHA BIOLOGICAL REFERENCE  (spec sections 2, 5)
##
## system > module > gene, and NOTHING below gene. Which CpGs exist is a
## property of the user's methylation data, not of the biology, so the
## biological reference carries no probe column and no probe-derived filter.
##
## This matters because the previous Alpha resources were built from RETAINED
## PROBE COVERAGE: a gene Project Alpha happened to have no usable probe for
## simply vanished from its system, which silently redefined the biology to
## match one cohort's array. The reference below is regenerated from the
## audited 2026c gene codebook instead, so all 1,282 genes remain selectable -
## including the 48 with zero retained Alpha probes. If another lab supplies a
## CpG mapping to one of those genes, DMSA can analyse it.
##
## Integrity is checked by DERIVING the counts from the shipped file rather
## than asserting literals, so a future codebook revision cannot silently
## disagree with a hard-coded number (2026c holds 187 modules, not the 188 an
## earlier draft assumed).
## ---------------------------------------------------------------------------

#' Integrity of the bundled Alpha biological reference
#'
#' Derives the system, module and gene counts from the shipped file and checks
#' the reference is well formed: every gene in exactly one system and module,
#' no probe-derived filtering, and polarity either curated or explicitly
#' unresolved. Counts are DERIVED, never asserted against literals.
#'
#' @param verbose Print the report? Default TRUE.
#' @return Invisibly, a list with the derived counts and the check results.
#' @examples
#' alpha_reference_check(verbose = FALSE)$counts
#' @export
alpha_reference_check <- function(verbose = TRUE) {
  ref <- alpha_reference()
  s <- ref$systems
  counts <- c(systems = length(unique(s$system_id)),
              modules = length(unique(stats::na.omit(s$module_id))),
              genes   = length(unique(s$gene)),
              rows    = nrow(s))
  chk <- list()
  chk[["one row per gene"]] <- counts[["rows"]] == counts[["genes"]]
  chk[["every gene has a system"]] <- !any(is.na(s$system_id))
  chk[["every gene has a module"]] <- !any(is.na(s$module_id))
  chk[["gene in exactly one system"]] <-
    max(tapply(s$system_id, s$gene, function(z) length(unique(z)))) == 1L
  chk[["gene in exactly one module"]] <-
    max(tapply(s$module_id, s$gene, function(z) length(unique(z)))) == 1L
  chk[["no probe/CpG column"]] <-
    !any(c("cpg", "probe", "probe_id", "column") %in% names(s))
  npol <- if (is.null(ref$polarity)) 0L else
    sum(ref$systems$gene %in% ref$polarity$gene)
  if (verbose) {
    cat(sprintf("Alpha biological reference (%s)\n", ref$version))
    cat(sprintf("  %d systems | %d modules | %d genes | %d rows\n",
                counts[["systems"]], counts[["modules"]], counts[["genes"]],
                counts[["rows"]]))
    for (nm in names(chk))
      cat(sprintf("  [%s] %s\n", if (chk[[nm]]) "ok" else "FAIL", nm))
    cat(sprintf("  polarity: %d of %d genes curated, %d unresolved (NA, never +1)\n",
                npol, counts[["genes"]], counts[["genes"]] - npol))
  }
  invisible(list(counts = counts, checks = chk, n_polarity = npol,
                 ok = all(unlist(chk))))
}
