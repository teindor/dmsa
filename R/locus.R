# ============================================================================
# LOCUS PANEL: where the CpGs are, and what each one did
#
# Up to three stacked strips for one gene:
#   A  the chromosome, with the gene's position marked   [needs coordinates]
#   B  the gene span, every measured CpG drawn as a lollipop - grey when it
#      carries no signal, and coloured DIFFERENTLY by its cpgdirection call d
#   C  the effect for each CpG with a 95% interval and the number printed,
#      in the same order as strip B
#
# DRAW THE GENE, NOT A BAR SHAPED LIKE ONE. Strip B takes `gene_model =`, a
# table of real exons, UTRs and strand from dmsa_gene_model(). Coding exons are
# full-height boxes, UTRs half-height, introns a line with chevrons running
# 5' to 3', and the TSS is marked - the conventions every genome browser uses,
# so a reader already knows how to read it. Earlier versions drew the gene as a
# featureless grey rectangle, which told the reader nothing while looking like
# it told them something. A CpG in the first exon and a CpG 40 kb into an
# intron are different claims about regulation, and the panel now shows which
# one it is.
#
# COORDINATES ARE STILL OPTIONAL, BUT NO LONGER FREE.
#   * With `pos`, strip B is drawn to genomic scale, strip A appears, and a
#     gene model can be drawn under the probes.
#   * With a `gene_region` column but no `pos`, the probes are ordered
#     TSS1500 -> TSS200 -> 5'UTR -> 1stExon -> Body -> 3'UTR and spaced evenly.
#   * With neither, the probes keep the order you supplied, spaced evenly.
# The last two are fallbacks, and the panel now says NOT TO SCALE in bold under
# the axis and messages once on the console when it takes one. Evenly-spaced
# probes look like a map and are not one: two probes drawn adjacent may be 40 kb
# apart, which is exactly the distinction that decides a probe's direction call.
# The bundled Alpha cascade carries an hg38 position for every probe, so the
# default path is always to scale.
#
# WHY BOTH ENGINES. Everything here is base graphics on a plain data frame, so
# the figure runs anywhere dmsa runs, with no Bioconductor install. When you
# want what Gviz gives and this does not - cytoband ideograms, transcript
# collapsing, and any other GdObject track stacked in the same frame -
# dmsa_plot_locus_gviz() takes the same probes and the same model and hands
# them to Gviz. Gviz is a Suggests; neither engine is required by the other.
#
# GETTING COORDINATES, if you want them: dmsa_probe_coords(probes) pulls them
# from cpgdirection, which is already installed because it is where the `d`
# calls come from - nothing new to install. It also reads any plain-text
# manifest with base R alone. dmsa_probe_annotation_template() writes the older
# minfi/Bioconductor recipe, which is about a gigabyte of annotation packages
# for two numbers per probe.
# ============================================================================

## hg19 chromosome lengths, used only to draw strip A when no length is given.
## Verify against your own annotation's seqinfo before publishing - pass
## `chrom_length` explicitly and this table is never touched.
.HG19_LEN <- c(
  "1" = 249250621, "2" = 243199373, "3" = 198022430, "4" = 191154276,
  "5" = 180915260, "6" = 171115067, "7" = 159138663, "8" = 146364022,
  "9" = 141213431, "10" = 135534747, "11" = 135006516, "12" = 133851895,
  "13" = 115169878, "14" = 107349540, "15" = 102531392, "16" = 90354753,
  "17" = 81195210, "18" = 78077248, "19" = 59128983, "20" = 63025520,
  "21" = 48129895, "22" = 51304566, "X" = 155270560, "Y" = 59373566)

## Canonical 5' -> 3' order of the UCSC RefGene groups. Used to lay probes out
## along the gene when no genomic coordinate is available.
.REGION_ORDER <- c("TSS1500", "TSS200", "5'UTR", "5UTR", "1stExon", "Body",
                   "ExonBnd", "3'UTR", "3UTR")

.region_rank <- function(x) {
  x <- sub(";.*$", "", trimws(as.character(x)))
  r <- match(x, .REGION_ORDER)
  r[is.na(r)] <- length(.REGION_ORDER) + 1L      # unknown regions sort last
  r
}

## ---------------------------------------------------------------------------
## The gene track: exons, introns, UTRs and strand, drawn to genomic scale.
##
## Conventions, which are the ones every genome browser uses and therefore the
## ones a reader already knows how to read:
##   - coding exon      full-height filled box
##   - UTR              half-height box, same colour
##   - intron           thin line with chevrons pointing 5' -> 3'
##   - TSS              a tick and an arrow at the transcript's 5' end
##   - unresolved exon  full-height box, hatched, because "we do not know which
##                      part codes" must not look like "all of it codes"
##
## Drawn into the CURRENT plot between y0 and y1, so the caller keeps control
## of the layout and the probes above it.
## ---------------------------------------------------------------------------
.locus_gene_track <- function(gm, y0, y1, col = "grey35", label = TRUE,
                              cex = .5) {
  if (is.null(gm) || !nrow(gm)) return(invisible(NULL))
  txs <- unique(gm$transcript)
  n <- length(txs)
  lane_h <- (y1 - y0) / n
  for (i in seq_along(txs)) {
    g <- gm[gm$transcript == txs[i], , drop = FALSE]
    ## lanes run top-down so the canonical transcript, sorted first, sits top
    top <- y1 - (i - 1) * lane_h
    bot <- top - lane_h
    mid <- (top + bot) / 2
    hf <- (top - bot) * .38          # half-height of a coding exon
    st <- g$strand[1]
    tx_s <- min(g$start); tx_e <- max(g$end)

    ## intron line across the whole transcript
    graphics::segments(tx_s, mid, tx_e, mid, col = col, lwd = 1.1)
    ## chevrons: direction is a property of the gene, not decoration. Spaced by
    ## device width so a 2 kb gene and a 200 kb gene get the same visual rate.
    usr <- graphics::par("usr")
    step <- (usr[2] - usr[1]) / 34
    if (is.finite(step) && step > 0 && (tx_e - tx_s) > step) {
      at <- seq(tx_s + step / 2, tx_e - step / 4, by = step)
      ## chevrons belong in introns. Drawn across an exon they read as texture
      ## on the box and hide its edges, which are the thing being shown.
      if (length(at)) {
        in_ex <- vapply(at, function(z) any(z >= g$start & z <= g$end), TRUE)
        at <- at[!in_ex]
      }
      if (length(at)) {
        dx <- step * .16 * (if (identical(st, "-")) -1 else 1)
        graphics::segments(at - dx / 2, mid - hf * .30, at + dx / 2, mid,
                           col = col, lwd = .9)
        graphics::segments(at + dx / 2, mid, at - dx / 2, mid + hf * .30,
                           col = col, lwd = .9)
      }
    }
    ## exons, UTRs last so a UTR never hides a coding box it abuts
    ord <- order(match(g$feature, c("cds", "exon", "utr5", "utr3")))
    for (k in ord) {
      f <- g$feature[k]
      h <- if (f %in% c("utr5", "utr3")) hf * .52 else hf
      dens <- if (f == "exon") 22 else NULL      # hatched = coding unresolved
      graphics::rect(g$start[k], mid - h, g$end[k], mid + h,
                     col = if (is.null(dens)) col else col, border = col,
                     density = dens, angle = 45, lwd = .7)
    }
    ## TSS: at the 5' end, which is the HIGH coordinate on the minus strand
    tss <- if (identical(st, "-")) tx_e else tx_s
    graphics::segments(tss, mid - hf, tss, top - lane_h * .04, col = col,
                       lwd = 1.3)
    arw <- (usr[2] - usr[1]) / 46 * (if (identical(st, "-")) -1 else 1)
    graphics::arrows(tss, top - lane_h * .04, tss + arw, top - lane_h * .04,
                     length = .035, col = col, lwd = 1.3)
    if (label && n > 1)
      graphics::text(usr[1], mid, txs[i], adj = c(0, .5), cex = cex * .82,
                     col = "grey45", xpd = NA)
  }
  invisible(NULL)
}

## What the panel says under strip B about where the gene model came from. A
## figure that draws exons must say whose exons they are.
.locus_model_note <- function(gm, transcripts) {
  if (is.null(gm) || !nrow(gm)) return(NULL)
  nt <- length(unique(gm$transcript))
  cn <- unique(gm$transcript[gm$canonical])
  unres <- any(gm$feature == "exon")
  paste0(gm$source[1], " ", gm$genome[1], " gene model",
         if (transcripts == "canonical" && length(cn))
           paste0(", canonical transcript ", cn[1]) else
         if (nt > 1) paste0(", all ", nt, " transcripts") else "",
         "; boxes exons, thin boxes UTRs, chevrons 5' to 3'",
         if (unres) "; hatched exons have unresolved coding status" else "")
}

#' Write a template script that fetches probe coordinates from Bioconductor
#'
#' The heavyweight route, kept for people who already have the annotation
#' packages installed. \code{\link{dmsa_probe_coords}} does the same job from a
#' plain-text manifest with base R only, and
#' \code{\link{dmsa_plot_locus}} needs no coordinates at all.
#'
#' Run the written script once; it produces a csv with \code{probe, chr, pos,
#' gene_region, island} that \code{dmsa_plot_locus()} consumes.
#'
#' @param path Where to write the script.
#' @param array \code{"EPIC"} or \code{"450K"}.
#' @return \code{path}, invisibly.
#' @seealso \code{\link{dmsa_probe_coords}}
#' @examples
#' f <- tempfile(fileext = ".R")
#' dmsa_probe_annotation_template(f, array = "EPIC")
#' ## the script is not run here: it installs and queries the minfi annotation
#' ## packages. dmsa_probe_coords() does the same job from a plain manifest.
#' head(readLines(f), 5)
#' @export
dmsa_probe_annotation_template <- function(path = "get_probe_coords.R",
                                           array = c("EPIC", "450K")) {
  array <- match.arg(array)
  pkg <- if (array == "EPIC")
    "IlluminaHumanMethylationEPICanno.ilm10b4.hg19" else
    "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  writeLines(c(
    "# Probe coordinates for dmsa_plot_locus(). Run once.",
    "# Lighter alternatives: dmsa_probe_coords(), or no coordinates at all.",
    "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')",
    sprintf("BiocManager::install(c('minfi', '%s'))", pkg),
    sprintf("library(%s); library(minfi)", pkg),
    sprintf("ann <- as.data.frame(minfi::getAnnotation(%s))", pkg),
    "out <- data.frame(probe = rownames(ann),",
    "                  chr = sub('^chr', '', ann$chr),",
    "                  pos = ann$pos,",
    "                  gene_region = ann$UCSC_RefGene_Group,",
    "                  island = ann$Relation_to_Island,",
    "                  stringsAsFactors = FALSE)",
    "write.csv(out, 'probe_coords.csv', row.names = FALSE)",
    "cat('wrote probe_coords.csv with', nrow(out), 'probes\\n')",
    "",
    "# gene spans, for the gene_start / gene_end arguments:",
    "# BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene'); ",
    "# BiocManager::install('org.Hs.eg.db')",
    "# library(TxDb.Hsapiens.UCSC.hg19.knownGene); library(org.Hs.eg.db)",
    "# g <- GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene)",
    "# sym <- AnnotationDbi::mapIds(org.Hs.eg.db, names(g), 'SYMBOL', 'ENTREZID')",
    "# genes_df <- data.frame(gene = sym, chr = sub('^chr','',as.character(GenomicRanges::seqnames(g))),",
    "#                        start = GenomicRanges::start(g), end = GenomicRanges::end(g))",
    "# write.csv(genes_df, 'gene_spans.csv', row.names = FALSE)"),
    path)
  message("wrote ", path, " - run it once, then join probe_coords.csv to your ",
          "probes. Lighter: dmsa_probe_coords(), or skip coordinates entirely.")
  invisible(path)
}

## The locus title, drawn as mtext rather than main= so a second line can sit
## under it.
##
## Strip C plots a beta per CpG on an axis labelled "raises / lowers expression
## tone". That says what the SIGN means. It does not say what the effect is an
## effect of - and "AVP   chr20" does not either. Reading the figure without
## the outcome is guessing, and a figure that goes into a paper cannot ask its
## reader to guess.
.locus_title <- function(gene, chrom, context = "") {
  ttl <- paste0(gene, if (!is.na(chrom) && nzchar(as.character(chrom)))
                        paste0("   chr", chrom) else "")
  two <- length(context) == 1L && !is.na(context) && nzchar(context)
  graphics::mtext(ttl, side = 3, line = if (two) 1.45 else 0.75,
                  cex = 1.0, font = 2)
  if (!two) return(invisible(NULL))
  ## Shrink rather than clip: an outcome label is the user's own column name
  ## or whatever they passed to labels =, and can be any length.
  cx <- 0.68
  while (cx > 0.40 &&
         graphics::strwidth(context, units = "inches", cex = cx) >
           graphics::par("din")[1] - 0.4)
    cx <- round(cx - 0.04, 2)
  graphics::mtext(context, side = 3, line = 0.45, cex = cx, col = "grey30")
  invisible(NULL)
}

#' Locus figure: chromosome, gene with CpG lollipops, and per-CpG effects
#'
#' Genomic coordinates are optional. Supply \code{pos} and the CpGs are placed
#' to scale and the chromosome strip is drawn; supply \code{gene_region} instead
#' and they are ordered along the gene (TSS1500 to 3'UTR) and spaced evenly;
#' supply neither and they keep the order you gave. The axis states which of
#' the three you got. Nothing else in the figure changes.
#'
#' @param probes data.frame for ONE gene. Required: \code{probe}, \code{b},
#'   \code{se}, \code{d} (cpgdirection, +1/-1). Optional: \code{pos} (genomic
#'   coordinate), \code{chr}, \code{gene_region}, \code{p} (per-probe adjusted
#'   p), \code{signal} (logical, overrides \code{p}).
#' @param gene Gene symbol, for labelling.
#' @param chrom Chromosome, e.g. \code{"20"}. Taken from a \code{chr} column
#'   when present and not given here.
#' @param gene_start,gene_end Gene span in the same coordinate system as
#'   \code{pos}. Defaults to the probe range padded by 10\%. Ignored when there
#'   are no coordinates.
#' @param chrom_length Length of the chromosome. Defaults to a built-in hg19
#'   table; pass your annotation's value to be sure. Strip A is omitted when
#'   neither is available.
#' @param signal_p Probes with \code{p < signal_p} are drawn as carrying signal.
#'   If \code{probes$p} is absent, all probes count as signal.
#' @param invert Which \code{d} is reflected onto the common effect axis in
#'   strip C, exactly as in \code{dmsa_plot_probes()}.
#' @param order_by \code{"auto"} uses \code{pos} if present, else
#'   \code{gene_region}, else the supplied row order. Force one with
#'   \code{"pos"}, \code{"region"} or \code{"given"}.
#' @param gene_model Real exon, UTR and strand structure to draw under the
#'   probes. A \code{\link{dmsa_gene_model}} table, a path to a GTF/GFF3,
#'   \code{TRUE} to fetch the model from Ensembl, or \code{NULL} for none.
#'   Requires genomic positions for the probes; a model on a different
#'   chromosome, or one that could not be obtained, is dropped with a message
#'   rather than drawn.
#' @param transcripts \code{"canonical"} draws the annotation's own canonical
#'   transcript and names it under the axis; \code{"all"} stacks every
#'   transcript in the model, one lane each.
#' @param context One line of context drawn under the title - in a DMSA report
#'   this is the outcome the effects belong to. Without it the panel names a
#'   gene and a chromosome but never says what strip C's effects are effects
#'   OF, which is not recoverable from the figure. Empty string draws nothing.
#' @param file Optional png path.
#' @param width,height,res png dimensions. \code{height = NULL} scales to the
#'   number of strips and probes.
#' @return Invisibly, the plotted data frame, with the columns the figure
#'   actually used (\code{x}, \code{sig}, \code{shown}, \code{flipped}).
#' @seealso \code{\link{dmsa_gene_model}} to build the model,
#'   \code{\link{dmsa_plot_locus_gviz}} for the Gviz engine,
#'   \code{\link{dmsa_probe_coords}} for probe positions.
#' @examples
#' probes <- data.frame(probe = paste0("cg", 1:5), chr = "chr20",
#'                      pos = c(3082600, 3082750, 3083050, 3084600, 3084700),
#'                      b = c(-.04, .02, .03, -.05, -.06), se = rep(.015, 5),
#'                      d = c(-1, 1, 1, -1, -1), p = c(.01, .3, .04, .002, .001))
#' ## invert = "+1" (the default) reflects the d = +1 probes, so a positive
#' ## effect always means "moved methylation the way that LOWERS expression"
#' out <- dmsa_plot_locus(probes, gene = "AVP", context = "attachment anxiety",
#'                        file = tempfile(fileext = ".png"))
#' out[, c("probe", "d", "shown", "flipped")]
#' @export
dmsa_plot_locus <- function(probes, gene = "", chrom = NA,
                            gene_start = NA, gene_end = NA,
                            chrom_length = NA, signal_p = .05,
                            invert = c("+1", "-1", "none"),
                            order_by = c("auto", "pos", "region", "given"),
                            gene_model = NULL,
                            transcripts = c("canonical", "all"),
                            context = "",
                            file = NULL, width = 1800, height = NULL,
                            res = 190) {
  invert <- match.arg(invert); order_by <- match.arg(order_by)
  transcripts <- match.arg(transcripts)

  ## ---- the gene model ------------------------------------------------------
  ## TRUE means "fetch it"; a data frame means "use this one"; NULL means the
  ## caller does not want one. A failed fetch returns zero rows and the panel
  ## degrades to a bare coordinate axis - it never falls back to drawing a
  ## featureless bar as if that were the gene.
  if (isTRUE(gene_model)) {
    gene_model <- if (nzchar(gene))
      tryCatch(dmsa_gene_model(gene, quiet = TRUE),
               error = function(e) NULL) else NULL
  } else if (is.character(gene_model) && length(gene_model) == 1L) {
    gene_model <- tryCatch(
      dmsa_gene_model(if (nzchar(gene)) gene else gene_model,
                      source = "gff", file = gene_model, quiet = TRUE),
      error = function(e) NULL)
  }
  if (!is.null(gene_model)) {
    gene_model <- as.data.frame(gene_model, stringsAsFactors = FALSE)
    if (!nrow(gene_model) || !all(c("start", "end", "feature", "transcript",
                                    "strand") %in% names(gene_model))) {
      gene_model <- NULL
    } else {
      if (transcripts == "canonical" && any(gene_model$canonical))
        gene_model <- gene_model[gene_model$canonical, , drop = FALSE]
      ## canonical first, so it takes the top lane
      gene_model <- gene_model[order(!gene_model$canonical,
                                     gene_model$transcript,
                                     gene_model$start), , drop = FALSE]
    }
  }
  p <- as.data.frame(probes, stringsAsFactors = FALSE)
  for (v in c("probe", "b", "se", "d")) if (!v %in% names(p))
    stop("probes needs a '", v, "' column", call. = FALSE)
  p <- p[is.finite(as.numeric(p$b)), , drop = FALSE]
  if (!nrow(p)) stop("no probe has a finite effect", call. = FALSE)

  ## ---- what layout information do we actually have? ------------------------
  ## A genomic axis is only used when EVERY probe can be placed on it. With
  ## partial coordinates the alternative is dropping measured CpGs from the
  ## figure, which is the one thing this panel must never do quietly.
  n_pos <- if ("pos" %in% names(p))
    sum(is.finite(suppressWarnings(as.numeric(p$pos)))) else 0L
  has_pos <- n_pos == nrow(p) && n_pos > 0L
  has_reg <- "gene_region" %in% names(p) &&
    any(nzchar(trimws(as.character(p$gene_region))))
  if (order_by == "auto" && n_pos > 0L && !has_pos)
    message(nrow(p) - n_pos, " of ", nrow(p), " probes have no position; ",
            "spacing all of them evenly instead of dropping them. ",
            "order_by = 'pos' forces the genomic axis and drops them.")
  ## An evenly-spaced layout is a fallback, never a silent default. It looks
  ## like a map and is not one: two probes drawn adjacent may be 40 kb apart,
  ## which is precisely the distinction that decides a probe's direction call.
  ## Say so once, loudly, whenever the panel is about to do it.
  if (order_by == "auto" && n_pos == 0L)
    message("no probe carries a genomic position, so this panel cannot be ",
            "drawn to scale: spacing is arbitrary and the x-axis is not ",
            "coordinates. Supply a `pos` column - dmsa_probe_coords() gets ",
            "one - or pass gene_model = TRUE for the exon structure too.")
  mode <- switch(order_by,
                 auto   = if (has_pos) "pos" else if (has_reg) "region" else "given",
                 pos    = if (n_pos > 0L) "pos" else
                   stop("order_by='pos' but there is no usable 'pos' column",
                        call. = FALSE),
                 region = if (has_reg) "region" else
                   stop("order_by='region' but there is no usable 'gene_region'",
                        " column", call. = FALSE),
                 given  = "given")

  if (mode == "pos") {
    p$pos <- suppressWarnings(as.numeric(p$pos))
    p <- p[is.finite(p$pos), , drop = FALSE]
    if (!nrow(p)) stop("no probe has both a finite position and effect",
                       call. = FALSE)
    p <- p[order(p$pos), , drop = FALSE]
    p$x <- p$pos
  } else {
    if (mode == "region")
      p <- p[order(.region_rank(p$gene_region)), , drop = FALSE]
    p$x <- seq_len(nrow(p))                     # evenly spaced, in that order
  }
  to_scale <- mode == "pos"

  p$d <- sign(suppressWarnings(as.numeric(p$d)))
  ## "signal" is the analyst's call, not a fixed rule: a logical `signal`
  ## column wins, then a p column, then everything counts. At probe level a
  ## maxT-adjusted p almost never fires, so defaulting to it would grey out
  ## every probe in a gene that was selected precisely because they agree.
  p$sig <- if (!is.null(p$signal)) as.logical(p$signal) else
    if (!is.null(p$p)) is.finite(p$p) & p$p < signal_p else rep(TRUE, nrow(p))
  p$sig[is.na(p$sig)] <- FALSE
  flip <- switch(invert, "+1" = p$d > 0, "-1" = p$d < 0,
                 "none" = rep(FALSE, nrow(p)))
  flip[is.na(flip)] <- FALSE
  p$shown <- ifelse(flip, -as.numeric(p$b), as.numeric(p$b)); p$flipped <- flip
  p$se <- as.numeric(p$se)

  VC <- .dmsa_pal(6)
  ## grey = no signal; two distinct colours for the two direction calls
  COL_UP <- VC[2]; COL_DN <- VC[5]; COL_NS <- "grey72"
  p$col <- ifelse(!p$sig, COL_NS, ifelse(p$d > 0, COL_UP, COL_DN))

  ## a model that does not sit on the same chromosome as the probes is the
  ## wrong gene: drop it rather than draw it
  if (!is.null(gene_model) && to_scale && !is.null(p$chr)) {
    pc <- unique(sub("^chr", "", as.character(p$chr)))
    mc <- unique(sub("^chr", "", as.character(gene_model$chr)))
    if (length(pc) == 1L && length(mc) == 1L && !identical(pc, mc)) {
      message("the gene model is on chr", mc, " but the probes are on chr", pc,
              " - not drawing it")
      gene_model <- NULL
    }
  }
  if (!is.null(gene_model) && !to_scale) {
    message("a gene model needs genomic positions for the probes; ",
            "drawing the probes without it")
    gene_model <- NULL
  }

  ## ---- the genome-build guard ---------------------------------------------
  ## Probe coordinates and the gene model can each be right and still disagree,
  ## because they came from different assemblies. An EPIC manifest is usually
  ## hg19; Ensembl and the bundled cascade are hg38. AVP moves 19,353 bp between
  ## them, which silently draws its promoter CpGs into the middle of the
  ## neighbouring intergene. Neither layer is wrong, so nothing errors - which
  ## is exactly why this has to be checked explicitly.
  ##
  ## The rule is distance, not overlap: real promoter probes sit OUTSIDE the
  ## gene span (all 8 of AVP's do) and must not be flagged. A build mismatch
  ## puts them further away than the gene is long.
  build_gap <- NULL
  if (!is.null(gene_model) && to_scale) {
    g1 <- min(gene_model$start); g2 <- max(gene_model$end)
    gap <- max(0, max(g1 - max(p$x), min(p$x) - g2))
    tol <- max(5000, g2 - g1)
    if (gap > tol) {
      build_gap <- gap
      gb <- unique(stats::na.omit(as.character(p$genome)))
      message("GENOME BUILD MISMATCH? every probe is ", format(round(gap),
              big.mark = ","), " bp from the ", gene_model$genome[1],
              " gene model, and the gene is only ",
              format(g2 - g1, big.mark = ","), " bp long.\n",
              "  Probe coordinates ", if (length(gb) == 1L)
                paste0("say they are ", gb) else "carry no genome label",
              ". An EPIC manifest is typically hg19 and the bundled cascade ",
              "is hg38.\n  The panel is drawn, and says so, but do not ",
              "publish it until the builds agree.")
    }
  }

  ## ---- the x window --------------------------------------------------------
  if (to_scale) {
    if (is.na(gene_start) || is.na(gene_end)) {
      ## With a model, the window is the gene: a panel that cropped to the
      ## probes would cut off the exons that explain where those probes sit.
      if (!is.null(gene_model)) {
        rg <- range(c(p$x, gene_model$start, gene_model$end))
      } else rg <- range(p$x)
      ## pad proportionally. A fixed floor of a few hundred bases would squash
      ## a tight CpG island - AVP's 8 probes span 241 bp - into the middle
      ## fifth of the strip.
      pad <- max(diff(rg) * .06, 30)
      gene_start <- rg[1] - pad; gene_end <- rg[2] + pad
    }
  } else {
    gene_start <- 0.4; gene_end <- nrow(p) + 0.6
  }
  if (is.na(chrom) && !is.null(p$chr)) {
    cc <- unique(sub("^chr", "", as.character(p$chr)))
    cc <- cc[!is.na(cc) & nzchar(cc)]
    if (length(cc) == 1L) chrom <- cc
  }
  if (is.na(chrom_length) && !is.na(chrom) &&
      as.character(chrom) %in% names(.HG19_LEN))
    chrom_length <- .HG19_LEN[[as.character(chrom)]]
  draw_A <- to_scale && !is.na(chrom) && is.finite(chrom_length)

  ## ---- device --------------------------------------------------------------
  if (is.null(height))
    height <- round((if (draw_A) 620 else 380) +
                    max(760, 62 * nrow(p)) + 640)
  if (!is.null(file)) { grDevices::png(file, width, height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE) }
  ## a gene model needs its own depth, and a stacked one needs more per lane
  hB <- if (!is.null(gene_model)) 2.15 + 0.55 * (length(unique(gene_model$transcript)) - 1L) else
    if (has_reg) 2.35 else 1.75
  if (draw_A) graphics::layout(matrix(1:3, ncol = 1),
                               heights = c(0.95, hB, 3.9)) else
                graphics::layout(matrix(1:2, ncol = 1),
                                 heights = c(hB + .3, 3.9))

  ## ---- A: the chromosome ---------------------------------------------------
  if (draw_A) {
    graphics::par(mar = c(1.4, 7.0, 3.0, 2.4))
    plot(NA, xlim = c(0, chrom_length), ylim = c(0, 1), axes = FALSE,
         xlab = "", ylab = "", main = "")
    .locus_title(gene, chrom, context)
    graphics::rect(0, .3, chrom_length, .7, col = "grey93", border = "grey55",
                   lwd = 1.2)
    mid <- (gene_start + gene_end) / 2
    graphics::rect(max(0, mid - chrom_length * .004), .18,
                   min(chrom_length, mid + chrom_length * .004), .82,
                   col = VC[1], border = NA)
    graphics::text(mid, .95, gene, cex = .66, font = 2, col = VC[1],
                   pos = if (mid > chrom_length / 2) 2 else 4, xpd = NA)
    graphics::axis(1, at = pretty(c(0, chrom_length), 6),
                   labels = paste0(round(pretty(c(0, chrom_length), 6) / 1e6),
                                   " Mb"),
                   cex.axis = .58, col = "grey60", col.axis = "grey40",
                   line = -1.2, tick = FALSE)
    graphics::mtext("A", side = 2, line = 5.4, las = 1, font = 2, cex = 1.0)
  }

  ## ---- B: the gene, with CpG lollipops ------------------------------------
  graphics::par(mar = c(if (to_scale) 3.6 else 3.2, 7.0,
                        if (draw_A) 1.6 else 3.0, 2.4))
  hgt <- ifelse(p$sig, 1, .55)
  ## a stacked model needs depth for its lanes; one transcript keeps the
  ## original proportions so nothing else in the figure shifts
  n_lane <- if (!is.null(gene_model)) length(unique(gene_model$transcript)) else 1L
  ## Exons in a large gene are a few hundred bases inside a hundred kilobases,
  ## so they are thin by construction - that is the truth of the locus, not a
  ## drawing fault. What the panel owes them is HEIGHT, so a 150 bp exon is
  ## still a legible box rather than a hairline.
  y_bot <- if (!is.null(gene_model)) -.30 - .30 * (n_lane - 1) else
    if (has_reg) -.42 else -.24
  plot(NA, xlim = c(gene_start, gene_end),
       ylim = c(min(y_bot, if (has_reg) -.42 else -.24), 1.30), axes = FALSE,
       xlab = "", ylab = "", main = "")
  if (!draw_A) .locus_title(gene, chrom, context)
  if (!is.null(gene_model)) {
    .locus_gene_track(gene_model, y0 = y_bot, y1 = .02, col = "grey38")
  } else {
    ## no model: a plain baseline, NOT a bar shaped like a gene
    graphics::segments(gene_start, -.08, gene_end, -.08, col = "grey65",
                       lwd = 1.1)
  }
  graphics::segments(p$x, .02, p$x, hgt, col = p$col, lwd = 1.6)
  graphics::points(p$x, hgt, pch = ifelse(p$flipped, 21, 19),
                   cex = ifelse(p$sig, 1.5, .9), col = p$col, bg = "white",
                   lwd = 1.8)
  if (to_scale) {
    graphics::axis(1, at = pretty(c(gene_start, gene_end), 5),
                   labels = format(pretty(c(gene_start, gene_end), 5),
                                   big.mark = ","),
                   cex.axis = .58, col = "grey60", col.axis = "grey40")
    if (!is.null(build_gap))
      graphics::mtext(sprintf(
        "CHECK GENOME BUILD - probes lie %s bp from this %s gene model",
        format(round(build_gap), big.mark = ","), gene_model$genome[1]),
        side = 3, line = -1.4, cex = .62, col = "firebrick", font = 2)
    ## A panel drawn without a model is correct but half the point: the caption
    ## says so, and captions get skimmed. Say it once on the console too, with
    ## the exact argument, so the fix does not require re-reading the figure.
    if (is.null(gene_model))
      message("no gene model drawn for ", if (nzchar(gene)) gene else "this unit",
              ": probes are on true coordinates but the exons are not shown. ",
              "Pass gene_models = TRUE to dmsa_frame() (fetches from Ensembl), ",
              "or gene_model = <table> from dmsa_gene_model() to stay offline.")
    note <- .locus_model_note(gene_model, transcripts)
    graphics::mtext(if (!is.null(note)) note else
      paste0("true ", if (!is.na(chrom)) paste0("chr", chrom, " ") else "",
             "coordinates; no gene model supplied ",
             "(gene_model = TRUE draws the exons)"),
      side = 1, line = 2.35, cex = .55, col = "grey35")
  } else {
    ## an arrow asserts 5' to 3'. Only the region mode has earned that claim.
    if (mode == "region") {
      yA <- if (has_reg) -.30 else -.21
      graphics::arrows(gene_start + .06, yA, gene_end - .06, yA, length = .06,
                       col = "grey55", lwd = 1.1, xpd = NA)
      graphics::text(c(gene_start, gene_end), yA, c("5'", "3'"), cex = .55,
                     col = "grey45", adj = c(.5, .5), xpd = NA)
    }
    ## Not an axis. Say that where an axis would be, in the figure itself, so
    ## the caveat survives being cropped into a slide.
    graphics::mtext(if (mode == "region")
      "NOT TO SCALE - ordered by RefGene region, spacing arbitrary" else
      "NOT TO SCALE - probes in the order supplied, spacing arbitrary",
      side = 1, line = 1.3, cex = .6, col = "grey20", font = 2)
  }
  graphics::mtext("B", side = 2, line = 5.4, las = 1, font = 2, cex = 1.0)
  graphics::mtext(sprintf("%d CpGs measured, %d carrying signal", nrow(p),
                          sum(p$sig)), side = 3, line = -.2, cex = .6,
                  col = "grey40")
  if (!is.null(p$gene_region)) {
    gr <- ifelse(nzchar(as.character(p$gene_region)),
                 sub(";.*$", "", as.character(p$gene_region)), "")
    ## stagger so neighbouring region labels never overprint
    graphics::text(p$x, ifelse(seq_len(nrow(p)) %% 2 == 0, -.24, -.36), gr,
                   srt = 90, adj = 1, cex = .48, col = "grey45", xpd = NA)
  }
  graphics::legend(gene_start, 1.32, bty = "n", cex = .62, horiz = TRUE,
    pch = 19, col = c(COL_NS, COL_UP, COL_DN), xjust = 0, yjust = 1,
    legend = c("no signal", "signal, d = +1", "signal, d = -1"))

  ## ---- C: the effects, in the same order, joined to their lollipops --------
  graphics::par(mar = c(6.2, 7.0, 1.0, 2.4))
  lo <- p$shown - 1.96 * p$se; hi <- p$shown + 1.96 * p$se
  xr <- range(c(lo, hi, 0), na.rm = TRUE)
  xr <- xr + c(-.10, .34) * diff(xr)          # room for the printed numbers
  n <- nrow(p)
  ## row 1 at the TOP, so reading C downwards matches reading B left to right
  yy <- rev(seq_len(n))
  plot(NA, xlim = xr, ylim = c(.4, n + .6), yaxt = "n", xlab =
    if (invert == "none") "slope on methylation" else
      "effect   <- raises expression tone  |  lowers expression tone ->",
    ylab = "")
  graphics::abline(v = 0, col = "grey40")
  graphics::axis(2, at = yy, labels = p$probe, las = 1, cex.axis = .55,
                 tick = FALSE)
  graphics::segments(lo, yy, hi, yy, col = p$col, lwd = 2)
  graphics::points(p$shown, yy, pch = ifelse(p$flipped, 21, 19),
                   cex = ifelse(p$sig, 1.25, .85), col = p$col, bg = "white",
                   lwd = 1.8)
  graphics::text(hi, yy, sprintf("  %+.3f", p$shown), adj = 0,
                 cex = .58, col = p$col, xpd = NA)
  graphics::axis(4, at = yy, labels = ifelse(p$d > 0, "+1", "-1"),
                 las = 1, cex.axis = .52, tick = FALSE, col.axis = "grey45")
  graphics::mtext("C", side = 2, line = 5.4, las = 1, font = 2, cex = 1.0)
  nf <- sum(p$flipped)
  graphics::mtext(sprintf(
    "bars 95%% CI; rows run top to bottom in the same order as strip B (%s). %s",
    switch(mode, pos = "genomic", region = "5' to 3' by region",
           given = "as supplied"),
    if (invert == "none") "raw slopes." else if (nf == 0)
      sprintf("no probe needed reflecting (all d = %s).",
              if (invert == "+1") "-1" else "+1") else
      sprintf("open symbols (%d) are reflected by their d call.", nf)),
    side = 1, line = 4.5, cex = .58, col = "grey35")
  invisible(p)
}
