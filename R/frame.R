# ============================================================================
# THE TWO-CALL INTERFACE, CALL 1: dmsa_frame()
#
# Declare everything once - data, levels, outcomes, covariates, blocks,
# moderation, map - and get back a validated "dmsa_frame" that dmsa_report()
# runs without further questions. The frame is where every silent data
# problem becomes a visible, documented correction: dmsa_frame() TEST-DRIVES
# the declaration (coercion audit, rank audit, block audit, outcome audit,
# B = 49 pilot) and records every action in frame$corrections.
#
# Moderation layout (Tsachi, 14 Aug 2026):
#   moderation = FALSE by default; mod = ""; mod2 = "".
#   moderation = TRUE + mod            -> tested term is frame x mod
#   moderation = TRUE + mod + mod2     -> frame x mod x mod2 (moderators
#     mean-centered, all lower-order terms in the model, the highest-order
#     term is the test). Interactions ride the COMPOSITE lens - the one lens
#     that carries to products - so moderation runs are composite-only and say
#     so in the report.
#
# cpg_map: "confidence" (default) keeps probes whose direction call meets the
# high-confidence bar (tier A or SMR S1); "full" keeps every probe with any
# usable direction call. BOTH alignments are built regardless; any gene or
# system whose pooled direction flips between the two maps is analysed under
# the chosen map but warned about, listed in frame$map_conflicts, and repeated
# in the report's Design notes.
# ============================================================================

## The four Project Alpha builds each document a DIFFERENT covariate contract in
## their own Covariates sheet. A single hard-coded set silently degraded three of
## the four to sex alone - and on the before/after builds that drops `time`,
## which is the entire design. Detect the build from the columns present and
## return its documented contract; refuse to guess when nothing matches.
.frame_contracts <- list(
  parent_T1 = list(
    need = c("age_at_array_T1", "Epi_T1", "Fib_T1", "ctrlSV3_T1", "ctrlSV5_T1"),
    fixed = c("sex_c", "age_at_array_T1", "Epi_T1", "Fib_T1",
              "ctrlSV3_T1", "ctrlSV5_T1"),
    block = c("cID"), chip = "chip_T1",
    label = "build 1: parent T1"),
  child_T4 = list(
    need = c("birth_week", "birth_weight", "Epi_T4"),
    fixed = c("sex_c", "birth_week", "birth_weight", "Epi_T4"),
    block = c("cID"), chip = "chip_T4",
    label = "build 4: child T4 (leanest set - n is small)"),
  before_after = list(
    need = c("time", "Epi_array", "Fib_array", "ctrlSV5_array"),
    fixed = c("time", "sex_c", "Epi_array", "Fib_array", "ctrlSV5_array"),
    block = c("ID", "cID"), chip = NULL,
    label = "build 2/3: before-after (repeated measures; chip is NEVER included)")
)

.frame_which_build <- function(nms) {
  hit <- vapply(.frame_contracts, function(z) all(z$need %in% nms), logical(1))
  if (!any(hit)) return(NULL)
  names(.frame_contracts)[which(hit)[1L]]
}

.frame_contract <- function(nms = NULL) {
  if (is.null(nms)) return(.frame_contracts$parent_T1$fixed)
  b <- .frame_which_build(nms)
  if (is.null(b))
    stop("covariates = \"contract\" could not identify which Project Alpha ",
         "build this data frame is: none of the documented contracts had its ",
         "columns present. Expected one of - parent T1 (age_at_array_T1, ",
         "Epi_T1, Fib_T1, ctrlSV3_T1, ctrlSV5_T1), child T4 (birth_week, ",
         "birth_weight, Epi_T4), or before-after (time, Epi_array, Fib_array, ",
         "ctrlSV5_array). Pass the covariates explicitly instead of ",
         "\"contract\".", call. = FALSE)
  .frame_contracts[[b]]$fixed
}

## BUILD-AWARE COLUMN NAMES. The same CpG is named differently in every
## Project Alpha workbook: cg..._TC21_T1_AVP in the parent build,
## cg..._TC21_AVP in the two before-after builds, and cg..._TC21_T4_AVP plus
## cg..._TC21_maternal_T1_AVP / _paternal_T1_AVP in the child build. The
## bundled direction map carries the parent-T1 name only, so before 1.16.0 a
## child or before-after frame matched ZERO probes and stopped with "no mapped
## methylation columns found" - three of the four builds could not be analysed
## at all. The selection cascade lists every variant explicitly; use it.
.frame_col_variants <- function(cs) {
  if (is.null(cs) || is.null(cs$cascade)) return(NULL)
  cd <- cs$cascade
  vs <- grep("^col_", names(cd), value = TRUE)
  if (!length(vs) || !"col_parent_T1" %in% vs) return(NULL)
  cd[, vs, drop = FALSE]
}

## Returns the column vector to use plus a note when it had to be translated.
## `build` is the contract-detected build; it breaks ties, because the child
## workbook carries THREE matching sets (the child's own T4 array and both
## parents' T1 arrays) and picking by match count alone would be arbitrary.
.frame_resolve_columns <- function(mp, avail, cs, build = NULL) {
  keep0 <- list(column = mp$column, note = NULL)
  tab <- .frame_col_variants(cs)
  if (is.null(tab) || !length(avail)) return(keep0)
  hit0 <- sum(mp$column %in% avail)
  i   <- match(mp$column, tab$col_parent_T1)
  if (all(is.na(i))) return(keep0)
  vs  <- names(tab)
  cnt <- vapply(vs, function(v) {
    x <- tab[[v]][i]; sum(!is.na(x) & x %in% avail)
  }, integer(1))
  if (!any(cnt > 0L)) return(keep0)
  best <- vs[which.max(cnt)]
  pref <- c(parent_T1 = "col_parent_T1", child_T4 = "col_child_T4",
            before_after = "col_long")
  if (!is.null(build) && build %in% names(pref) && pref[[build]] %in% vs &&
      cnt[[pref[[build]]]] > 0L) best <- pref[[build]]
  if (cnt[[best]] <= hit0) return(keep0)
  new  <- tab[[best]][i]
  bad  <- is.na(new)
  new[bad] <- mp$column[bad]
  other <- vs[vs != best & cnt[vs] > 0L]
  list(column = new, note = sprintf(
    paste0("the direction map names probes in the parent-T1 form, this data ",
           "frame uses the '%s' form -> translated, %d probe(s) matched%s"),
    sub("^col_", "", best), cnt[[best]],
    if (length(other))
      paste0(" (this workbook also carries: ",
             paste(sub("^col_", "", other), collapse = ", "),
             " - pass `methylation` as the column names you want to override ",
             "the default)")
    else ""))
}

## SHARED PROBES. 229 probes in the Alpha panel sit inside two genes and
## appear as a column under each name (IGF2+INS, GABRA5+GABRB3,
## HLA-DPA1+HLA-DPB1, PEG3+ZIM2 and others). Two such genes in the same
## level-local family are NOT two independent units - they are the same
## measurement twice - and if both survive they read as two findings. Say so
## at frame time, while the design can still be changed.
.frame_shared_probes <- function(map, cs = NULL) {
  need <- c("probe", "gene", "system_id")
  if (!all(need %in% names(map))) return(NULL)
  cd <- if (is.null(cs)) NULL else cs$cascade
  ## The cascade states sharing EXPLICITLY in `shared_with`. Inferring it from
  ## duplicated probe names in the analysed map does not work: an Alpha column
  ## name embeds its gene, so one probe under two genes is two differently
  ## named columns, and after the map is de-duplicated on `column` the
  ## duplication is invisible.
  if (is.null(cd) || !"shared_with" %in% names(cd)) return(NULL)
  sw <- ifelse(is.na(cd$shared_with), "", as.character(cd$shared_with))
  ok <- nzchar(sw)
  if (!any(ok)) return(NULL)
  cd <- cd[ok, , drop = FALSE]; sw <- sw[ok]
  pk <- if ("cpg" %in% names(cd)) cd$cpg else cd$probe_id
  if (is.null(pk)) return(NULL)

  ## "IGF2; INS" is two partners, not one gene called "IGF2; INS"
  parts <- strsplit(sw, ";", fixed = TRUE)
  n <- lengths(parts)
  long <- data.frame(sid   = rep(as.character(cd$system_id), n),
                     probe = rep(as.character(pk), n),
                     g1    = rep(as.character(cd$gene), n),
                     g2    = trimws(unlist(parts, use.names = FALSE)),
                     stringsAsFactors = FALSE)
  long <- long[nzchar(long$g2) & long$g1 != long$g2, , drop = FALSE]
  if (!nrow(long)) return(NULL)

  ## Only sharing that actually reached the analysis counts: the probe must have
  ## entered under BOTH gene names, in a system that was selected.
  ana  <- unique(data.frame(sid   = as.character(map$system_id),
                            gene  = as.character(map$gene),
                            probe = as.character(map$probe),
                            stringsAsFactors = FALSE))
  keyp <- paste(ana$sid, ana$gene, ana$probe)
  long <- long[paste(long$sid, long$g1, long$probe) %in% keyp &
               paste(long$sid, long$g2, long$probe) %in% keyp, , drop = FALSE]
  if (!nrow(long)) return(NULL)

  syl  <- if ("system" %in% names(map))
            stats::setNames(as.character(map$system),
                            as.character(map$system_id)) else NULL
  pair <- paste0(ifelse(is.null(syl), long$sid, syl[long$sid]), "|",
                 pmin(long$g1, long$g2), " + ", pmax(long$g1, long$g2))
  u <- unique(data.frame(pair = pair, probe = long$probe,
                         stringsAsFactors = FALSE))
  tb <- table(u$pair)
  out <- data.frame(pair = names(tb), n_shared = as.integer(tb),
                    stringsAsFactors = FALSE)
  out[order(-out$n_shared, out$pair), , drop = FALSE]
}



## OUTCOME KIND. The Results prose used one template - "Higher X was associated
## with ..." - which assumes a CONTINUOUS outcome. It produced "Higher sex was
## associated with ..." and "Higher time was associated with ...", both
## meaningless. Three cases have to be told apart, and by the DATA, not by the
## column name:
##   continuous - the template is correct
##   binary     - a between-group contrast (sex)
##   wave       - a WITHIN-PERSON change (time, in the before/after builds),
##                which is what makes those builds worth having
## `wave` requires two things: a before/after contract, and an outcome that
## actually varies inside the finest permutation block. Testing only "binary +
## has blocks" would misclassify sex in the parent build, where a family block
## (cID) legitimately contains both parents.
.frame_outcome_kind <- function(base, oc, blockv, build) {
  v <- base[[oc]]
  u <- sort(unique(v[is.finite(v)]))
  if (length(u) != 2L) return(list(kind = "continuous", lo = NA, hi = NA))
  within <- FALSE
  if (identical(build, "before_after") && length(blockv)) {
    id <- blockv[[1]]
    if (id %in% names(base)) {
      k <- tapply(v, base[[id]], function(z) length(unique(z[is.finite(z)])))
      within <- isTRUE(mean(k > 1, na.rm = TRUE) > 0.5)
    }
  }
  list(kind = if (within) "wave" else "binary", lo = u[1], hi = u[2])
}

## non-empty and not NA - `chip` can legitimately be NA (the before/after
## contract forbids chip entirely), and nzchar(NA) is NA, not FALSE.
.nz <- function(x) !is.null(x) && length(x) >= 1L && !is.na(x[1]) && nzchar(x[1])

.frame_note <- function(cor, field, issue, action)
  rbind(cor, data.frame(field = field, issue = issue, action = action,
                        stringsAsFactors = FALSE))

## ---- the bundled Alpha map -------------------------------------------------
.frame_read_map <- function(map = "alpha") {
  if (is.data.frame(map)) return(as.data.frame(map))
  f <- system.file("extdata", "coverage_v4_full.csv", package = "dmsa")
  if (!nzchar(f) || !file.exists(f))
    stop("bundled Alpha map not found; pass `map` as a data.frame with ",
         "columns probe, column, gene, system_id, system, best_direction, ",
         "p_plus (and optionally best_tier, smr_tier, module)", call. = FALSE)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

.frame_maps <- function(map = "alpha") {
  mp <- .frame_read_map(map)
  need <- c("probe", "column", "gene", "system_id", "system",
            "best_direction", "p_plus")
  miss <- setdiff(need, names(mp))
  if (length(miss))
    stop("map is missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  mp <- mp[!duplicated(mp$column), , drop = FALSE]
  usable <- is.finite(mp$best_direction)
  hi <- usable
  if (!is.null(mp$best_tier) || !is.null(mp$smr_tier)) {
    ta <- if (!is.null(mp$best_tier)) mp$best_tier %in% "A" else FALSE
    s1 <- if (!is.null(mp$smr_tier))  mp$smr_tier  %in% "S1" else FALSE
    hi <- usable & (ta | s1)
  }
  list(full = mp[usable, , drop = FALSE],
       confidence = mp[hi, , drop = FALSE])
}

.frame_modules <- function(sets = "alpha") {
  ## The selection cascade is the authority on module membership when one is
  ## available: it is the same file the user selects systems from, so the family
  ## a unit is corrected within cannot drift from the family it was declared in.
  ## modules_alpha.csv is kept only as a fallback for installs without a cascade.
  cs <- .frame_sets(sets)
  if (!is.null(cs)) {
    mm <- unique(cs$cascade[, c("system_id", "module_id", "module", "gene")])
    return(mm[order(.cas_num(mm$system_id), .cas_modnum(mm$module_id),
                    mm$gene), , drop = FALSE])
  }
  f <- system.file("extdata", "modules_alpha.csv", package = "dmsa")
  if (!nzchar(f) || !file.exists(f)) return(NULL)
  utils::read.csv(f, stringsAsFactors = FALSE)
}

## Load a cascade for use inside dmsa_frame(), tolerantly: a frame must still be
## constructible when no cascade is installed or a user passes cascade = NULL.
.frame_sets <- function(sets = "alpha") {
  if (is.null(sets) || isFALSE(sets)) return(NULL)
  if (inherits(sets, "dmsa_sets")) return(sets)
  if (inherits(sets, "dmsa_selection"))
    return(structure(list(cascade = sets$cascade,
                          systems = sets$systems,
                          modules = sets$modules,
                          source = sets$source, name = sets$name,
                          columns = grep("^col_", names(sets$cascade),
                                         value = TRUE)),
                     class = "dmsa_sets"))
  out <- try(dmsa_sets(sets), silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

## `systems =` accepts the cascade's short names ("hpa"), ids, or full names.
## Resolution runs through the cascade first so that the handle a user typed is
## the handle the cascade documents; the old exact-match behaviour is the
## fallback, so existing scripts keep working.
.frame_pick_systems <- function(systems, sys_tab, cs) {
  if (is.null(systems)) return(sys_tab)
  if (!is.null(cs)) {
    ## Only route through the cascade when the direction map is describing the
    ## same panel: system NAMES must overlap. Otherwise a numeric `systems = 2`
    ## against an unrelated map would silently pick up the cascade's system 2.
    same_panel <- any(sys_tab$system %in% cs$systems$system)
    got <- if (same_panel) try(.cas_resolve_systems(cs, systems), silent = TRUE)
           else structure("skip", class = "try-error")
    if (!inherits(got, "try-error") && nrow(got)) {
      keep <- sys_tab$system %in% got$system
      if (any(keep)) {
        out <- sys_tab[keep, , drop = FALSE]
        attr(out, "short") <- got$system_short[match(out$system, got$system)]
        attr(out, "absent") <- got$system_short[!got$system %in% out$system]
        return(out)
      }
      stop("`systems` resolved to ", paste(got$system_short, collapse = ", "),
           " but no mapped probe of those system(s) is present in the data",
           call. = FALSE)
    }
  }
  keep <- if (is.numeric(systems)) sys_tab$system_id %in% systems else
    sys_tab$system %in% systems |
    as.character(sys_tab$system_id) %in% as.character(systems)
  if (!any(keep))
    stop("`systems` matched nothing; available: ",
         paste(sprintf("%s '%s'", sys_tab$system_id, sys_tab$system),
               collapse = "; "),
         if (!is.null(cs)) paste0("\nshort names: ",
           paste(cs$systems$system_short, collapse = ", ")) else "",
         call. = FALSE)
  sys_tab[keep, , drop = FALSE]
}

## pooled direction of a unit under one map: sign of sum d * (2 p_plus - 1)
.frame_pooled_sign <- function(mp, by) {
  m <- mp$best_direction * (2 * pmin(pmax(mp$p_plus, 0), 1) - 1)
  m[!is.finite(m)] <- 0
  s <- tapply(m, mp[[by]], function(v) sign(sum(v)))
  s[is.na(s)] <- 0
  s
}

.frame_conflicts <- function(maps, systems_keep) {
  out <- data.frame(level = character(), unit = character(),
                    sign_confidence = integer(), sign_full = integer(),
                    stringsAsFactors = FALSE)
  for (lv in c("gene", "system")) {
    a <- .frame_pooled_sign(maps$confidence[maps$confidence$system_id %in%
                                              systems_keep, ], lv)
    b <- .frame_pooled_sign(maps$full[maps$full$system_id %in%
                                        systems_keep, ], lv)
    u <- intersect(names(a), names(b))
    bad <- u[a[u] != 0 & b[u] != 0 & a[u] != b[u]]
    if (length(bad))
      out <- rbind(out, data.frame(level = lv, unit = bad,
                                   sign_confidence = as.integer(a[bad]),
                                   sign_full = as.integer(b[bad]),
                                   stringsAsFactors = FALSE))
  }
  out
}

## ---- beta -> M when needed -------------------------------------------------
.frame_mvalues <- function(BE, eps = 1e-4) {
  rng <- range(BE, na.rm = TRUE)
  if (rng[1] >= 0 && rng[2] <= 1) {
    BE <- pmin(pmax(BE, eps), 1 - eps)
    list(M = log2(BE / (1 - BE)), converted = TRUE)
  } else list(M = BE, converted = FALSE)
}

#' Declare a DMSA analysis frame
#'
#' Call 1 of the two-call interface. Declares the data, the levels
#' (system / module / gene / probe), the outcomes, the covariates and blocks,
#' optional moderation, and the direction map - validates all of it, runs the
#' test drive (with \code{autofix = TRUE} correcting what is safely
#' correctable, LOUDLY), and returns a \code{dmsa_frame} for
#' \code{\link{dmsa_report}}.
#'
#' @param data data.frame holding outcomes, covariates, ids and (by default)
#'   the Alpha methylation columns.
#' @param methylation NULL (auto-detect the bundled map's columns in
#'   \code{data}), a character vector of column names, or a numeric matrix
#'   whose colnames are probe/column ids.
#' @param map "alpha" for the bundled Project Alpha map, or a data.frame.
#' @param system,module,gene,probe logical: which levels to analyse.
#' @param systems Optional subset of systems to analyse. Accepts the selection
#'   cascade's short names - \code{systems = c("hpa", "oxytocin")} - and also
#'   system ids or full names. Matching is case-insensitive and a unique prefix
#'   of a short name is enough. \code{NULL} (default) analyses every covered
#'   system. Everything below the system - modules, genes, probes - is taken in
#'   full; use \code{dmsa_select()} and pass the result as \code{sets} to
#'   narrow a lower level. Run \code{dmsa_systems()} to print the short names.
#' @param sets The selection cascade declaring system > module > gene > probe.
#'   \code{"alpha"} (default) uses the bundled, module-audited Project Alpha
#'   2026c cascade; a path, \code{data.frame}, \code{dmsa_sets} or
#'   \code{dmsa_selection} may be supplied instead, and \code{NULL} disables it
#'   (module membership then falls back to the legacy bundled module map). See
#'   \code{dmsa_sets_template()} for the schema.
#' @param outcome character vector of outcome column(s). With several, each is
#'   tested with the others as covariates (mutual adjustment). When
#'   \code{frame_role = "outcome"} these are the exposure/predictor column(s).
#' @param covariates "contract" for the Alpha student contract, or a character
#'   vector. chip is added as a fixed factor automatically when present.
#' @param random_effects column(s) defining exchangeable blocks (permutation
#'   blocks); default "cID".
#' @param chip how to handle the array/chip batch factor. \code{TRUE} (default)
#'   enters it as a fixed effect, using \code{chip_f} if present and otherwise
#'   building it from a \code{chip}/\code{chip_T*} column. \code{FALSE} leaves
#'   it out entirely - on the Alpha parent build the factor has 61 levels on
#'   about 400 arrays, so it costs 60 degrees of freedom that the control-probe
#'   surrogate variables already partly absorb, and 7 of those 61 chips hold a
#'   single array - fitted exactly by their own dummy (leverage 1, residual 0),
#'   so they are counted in \code{n} while contributing nothing to any estimate.
#'   \code{"pool"} keeps the batch adjustment but merges every level holding
#'   fewer than 3 arrays into one \code{chip_small} level, which recovers the
#'   degrees of freedom and removes the leverage-1 rows. Any other character
#'   string names the column to use. The before/after contracts forbid chip and
#'   override this argument.
#' @param chip_effect how the chip factor enters the model, once \code{chip}
#'   has named it. \code{"random"} (default, and what the Alpha covariate
#'   contract specifies) takes chip out of the fixed design and fits it as a
#'   one-way random intercept, \code{(1 | chip)}, by REML quasi-demeaning: a
#'   median of about 16 effective degrees of freedom on the parent build, and
#'   zero on the probes where the chip variance is estimated as nil.
#'   \code{"fixed"} restores the pre-1.18.0 behaviour and enters chip as a
#'   fixed factor, at one degree of freedom per chip beyond the first.
#'   \code{"none"} drops the batch adjustment altogether, whatever
#'   \code{chip} says. Validity does not rest on the variance component being
#'   right: the permutation null is computed under the same transform, so the
#'   random-intercept fit moves power, not type-I error.
#' @param outcome_levels Optional length-2 character vector naming the two
#'   levels of a two-level outcome, e.g. \code{c("T1", "T4")} or
#'   \code{c("female", "male")}. Used only in the written report, so a
#'   sentence can say which level is which instead of falling back to
#'   \code{outcome = 0} / \code{outcome = 1}. \code{NULL} (default) uses
#'   that fallback. Purely cosmetic: no level name reaches the model, the
#'   permutation or the family-wise correction.
#' @param moderation,mod,mod2 moderation switch and moderator column(s); see
#'   the layout note above.
#' @param outcome_label,predictor_labels,mod_label,mod2_label Optional display
#'   labels used in figures and in the written report instead of the raw column
#'   names. Either one per column, in order, or a named vector keyed by column
#'   name. Purely cosmetic: no label reaches the model, the permutation or the
#'   family-wise correction.
#' @param type "linear", "non-linear" (linear/quadratic/threshold arms,
#'   ACAT-combined), or "exponential" (adds an exp arm).
#' @param outcome_type "gaussian", "logistic", or "multinomial".
#' @param frame_role "predictor" (outcome ~ frame + covariates) or "outcome"
#'   (named columns predict the methylation frame).
#' @param cpg_map "confidence" or "full"; see the note above.
#' @param B,alpha,seed,correction inference settings (correction is applied
#'   within level-local families only).
#' @param progress Show a 0-100% progress bar while \code{dmsa_report()} runs.
#'   Defaults to \code{interactive()}: a bar writes carriage returns, which make
#'   a piped log unreadable, so it is off under \code{Rscript} and in tests.
#' @param beep Sound a completion signal when \code{dmsa_report()} finishes.
#'   \code{TRUE} uses \pkg{beepr} sound 8; a number picks another sound;
#'   \code{FALSE} is silent. Defaults to \code{interactive()}, so a run never
#'   makes noise during \code{R CMD check} or a testthat run. \pkg{beepr} is a
#'   Suggests: without it this is silently a no-op.
#' @param gene_models Draw real exon/intron structure under the probes in each
#'   locus panel. \code{TRUE} fetches each gene's model from Ensembl, which
#'   makes network calls, so it is off by default; you may also pass a table
#'   built once with \code{\link{dmsa_gene_model}} and reused for every gene.
#' @param palette,plots,tables,summary,plot_type,table_type,outdir output
#'   styling for \code{dmsa_report}.
#' @param autofix TRUE: the test drive corrects what is safely correctable and
#' @param weighting Character. Probe weighting engine within a unit: \code{"combined"} (default) fuses the flat and reliability statistics on one shared permutation stream, \code{"flat"} weights every usable aligned probe equally, \code{"reliability"} weights each probe by its item-rest correlation with the rest of its unit. Weights are computed from methylation alone, so the permutation null is unaffected.
#' @param w_floor Numeric. Lower bound applied to reliability weights before normalisation, so a single poorly behaved probe cannot be driven to zero influence. Ignored when \code{weighting = "flat"}.
#'   records it; FALSE: strict errors instead.
#' @return object of class \code{dmsa_frame}.
#' @examples
#' set.seed(42)
#' map <- data.frame(gene = rep(paste0("G", 1:3), each = 3), system_id = 1L,
#'                   system = "Sim system")
#' map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
#' map$column <- paste0(map$probe, "_", map$gene)
#' map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
#' map$p_plus <- ifelse(map$best_direction > 0, 0.9, 0.1)
#' d <- data.frame(out1 = rnorm(60), cov1 = rnorm(60),
#'                 cID = rep(1:30, each = 2))
#' sig <- 0.5 * outer(d$out1, map$best_direction * (map$gene == "G1"))
#' d[map$column] <- plogis(matrix(rnorm(60 * nrow(map)), 60) + sig)
#' dmsa_frame(d, map = map, outcome = "out1", covariates = "cov1",
#'            random_effects = "cID", B = 99, outdir = tempfile("dmsa_ex"))
#' @export
dmsa_frame <- function(data, methylation = NULL, map = "alpha",
                       system = TRUE, module = FALSE, gene = TRUE,
                       probe = TRUE, systems = NULL, sets = "alpha",
                       outcome,
                       covariates = "contract", random_effects = "cID",
                       chip = TRUE,
                       chip_effect = c("random", "fixed", "none"),
                       outcome_levels = NULL,
                       moderation = FALSE, mod = "", mod2 = "",
                       outcome_label = NULL, predictor_labels = NULL,
                       mod_label = NULL, mod2_label = NULL,
                       type = c("linear", "non-linear", "exponential"),
                       outcome_type = c("gaussian", "logistic", "multinomial"),
                       frame_role = c("predictor", "outcome"),
                       cpg_map = c("confidence", "full"),
                       B = 1999, alpha = 0.05, seed = 1,
                       correction = c("maxT", "minP"),
                       weighting = c("combined", "reliability", "flat"),
                       w_floor = 1.5,
                       palette = "viridis", plots = TRUE, tables = TRUE,
                       summary = TRUE, gene_models = FALSE,
                       progress = interactive(), beep = interactive(),
                       plot_type = c("png", "pdf"),
                       table_type = c("html", "docx", "rtf"),
                       outdir = "dmsa_output", autofix = TRUE) {
  chip_effect <- match.arg(chip_effect)
  type <- match.arg(type); outcome_type <- match.arg(outcome_type)
  frame_role <- match.arg(frame_role); cpg_map <- match.arg(cpg_map)
  correction <- match.arg(correction); plot_type <- match.arg(plot_type)
  weighting <- match.arg(weighting)
  table_type <- match.arg(table_type)
  data <- as.data.frame(data)
  cor <- data.frame(field = character(), issue = character(),
                    action = character(), stringsAsFactors = FALSE)

  ## ---- moderation layout validation ---------------------------------------
  if (isTRUE(moderation) && !nzchar(mod))
    stop("moderation = TRUE requires `mod` to name a column ",
         "(mod2 is optional and only together with mod)", call. = FALSE)
  if (nzchar(mod2) && !nzchar(mod))
    stop("`mod2` was set without `mod`; set `mod` first", call. = FALSE)
  if (!isTRUE(moderation) && (nzchar(mod) || nzchar(mod2))) {
    cor <- .frame_note(cor, "moderation",
                       "mod/mod2 given but moderation = FALSE",
                       "moderators ignored; set moderation = TRUE to test them")
  }
  for (mv in c(mod, mod2)) if (nzchar(mv) && !mv %in% names(data))
    stop("moderator column '", mv, "' is not in `data`", call. = FALSE)

  ## A nominal outcome cannot enter the main path. DMSA regresses methylation ON
  ## the outcome and tests ONE coefficient, so a k-level nominal variable coerced
  ## to numeric is silently analysed as an ordered score - "mode 5 is five times
  ## mode 1". That is not a defensible model and it must not run.
  if (outcome_type == "multinomial")
    stop("outcome_type = 'multinomial' is not supported by dmsa_frame(): the ",
         "main path regresses methylation ON the outcome and tests a single ",
         "1-df coefficient, so a nominal outcome would be coerced to a numeric ",
         "score and analysed as if its levels were equally spaced and ordered. ",
         "Use dmsa_outcome(), which implements a rank-based scheme for ",
         "multinomial and ordinal outcomes, or recode to a binary contrast ",
         "(e.g. one level vs the rest) and pass it as a two-level outcome.",
         call. = FALSE)

  ## ---- display labels -----------------------------------------------------
  ## Column names are the user's own and often unreadable in a figure
  ## ("Attachment_Anxiety_General_T1"). Labels are cosmetic only: nothing in the
  ## model, the permutation or the correction ever sees them, so a mislabelled
  ## run is still a correct run - but every panel and sentence can read in
  ## plain language.
  .lab_pairs <- function(cols, labs, what) {
    if (is.null(labs) || !length(labs)) return(character(0))
    labs <- as.character(labs)
    if (!is.null(names(labs)) && all(nzchar(names(labs)))) {
      unknown <- setdiff(names(labs), cols)
      if (length(unknown))
        warning("`", what, "` names not in ", what, ": ",
                paste(unknown, collapse = ", "), " - ignored", call. = FALSE)
      labs <- labs[intersect(names(labs), cols)]
      return(labs)
    }
    if (length(labs) != length(cols))
      stop("`", what, "` has ", length(labs), " label(s) for ", length(cols),
           " column(s); supply one per column or a named vector",
           call. = FALSE)
    stats::setNames(labs, cols)
  }

  ## ---- outcomes -----------------------------------------------------------
  if (missing(outcome) || !length(outcome))
    stop("`outcome` must name at least one column", call. = FALSE)
  miss <- setdiff(outcome, names(data))
  if (length(miss))
    stop("outcome column(s) not in `data`: ", paste(miss, collapse = ", "),
         call. = FALSE)

  ## ---- maps and coverage --------------------------------------------------
  maps <- .frame_maps(map)
  cs   <- .frame_sets(sets)

  ## translate the map's column names to whichever build this data frame is,
  ## BEFORE the intersection that decides whether any probe was found at all
  .avail <- if (is.matrix(methylation)) colnames(methylation)
            else if (is.character(methylation) && length(methylation) > 1)
              intersect(names(data), methylation)
            else names(data)
  .res <- .frame_resolve_columns(maps[[cpg_map]], .avail, cs,
            if (is.character(methylation) && length(methylation) > 1) NULL
            else .frame_which_build(names(data)))
  if (!is.null(.res$note)) {
    for (.k in names(maps))
      maps[[.k]]$column <- .frame_resolve_columns(maps[[.k]], .avail, cs,
        if (is.character(methylation) && length(methylation) > 1) NULL
        else .frame_which_build(names(data)))$column
    cor <- .frame_note(cor, "column names", "build-specific CpG naming",
                       .res$note)
  }

  chosen_map <- maps[[cpg_map]]
  in_data <- chosen_map$column %in% names(data)
  if (is.character(methylation) && length(methylation) > 1)
    in_data <- in_data & chosen_map$column %in% methylation
  if (is.matrix(methylation)) {
    keyed <- chosen_map$column %in% colnames(methylation) |
             chosen_map$probe  %in% colnames(methylation)
    in_data <- keyed
  }
  mp <- chosen_map[in_data, , drop = FALSE]
  if (!nrow(mp))
    stop("no mapped methylation columns found (cpg_map = '", cpg_map,
         "'); check `methylation`/`data` column names against the map",
         call. = FALSE)
  ## (untestable-gene tracking happens AFTER the systems subset below - a
  ## 1.1.0 bug counted untestable genes across ALL systems in the data)

  sys_tab <- unique(mp[, c("system_id", "system")])
  if (!is.null(systems)) {
    sys_tab <- .frame_pick_systems(systems, sys_tab, cs)
    gone <- attr(sys_tab, "absent")
    if (length(gone))
      cor <- .frame_note(cor, "systems", "selected system carries no probe here",
                         paste("no mapped probe in the data for:",
                               paste(gone, collapse = ", ")))
    mp <- mp[mp$system_id %in% sys_tab$system_id, , drop = FALSE]
  }
  mp <- mp[order(mp$system_id, mp$gene, mp$probe), , drop = FALSE]

  ## genes present in the raw map for the CHOSEN systems but with ZERO probes
  ## surviving the chosen direction bar are UNTESTABLE - tracked, reported,
  ## never silently dropped
  ## The raw map is re-read from disk, so it still carries the parent-T1
  ## column names - the build translation applied to `maps` above is NOT in it.
  ## Testing raw$column against names(data) therefore matched NOTHING in the
  ## child and before-after builds, `ut` came out empty, and every one of those
  ## reports claimed "Coverage: 0 gene(s) ... untestable" while the parent build
  ## on the SAME two systems and the SAME 18 genes correctly reported 36.
  ## Translate raw the same way before the intersection.
  raw <- .frame_read_map(map)
  raw$column <- .frame_resolve_columns(
    raw, .avail, cs,
    if (is.character(methylation) && length(methylation) > 1) NULL
    else .frame_which_build(names(data)))$column
  raw_in <- if (is.matrix(methylation))
              raw$column %in% colnames(methylation) |
              raw$probe  %in% colnames(methylation)
            else if (is.character(methylation) && length(methylation) > 1)
              raw$column %in% names(data) & raw$column %in% methylation
            else raw$column %in% names(data)
  raw <- raw[raw$system_id %in% sys_tab$system_id & raw_in, , drop = FALSE]
  ut <- setdiff(unique(raw$gene), unique(mp$gene))
  untestable <- unique(raw[raw$gene %in% ut, c("system_id", "system", "gene")])

  ## modules requested -> module map must cover the chosen systems
  mods <- NULL
  if (isTRUE(module)) {
    mm <- .frame_modules(sets)
    if (!is.null(mm) && !any(mm$gene %in% mp$gene)) mm <- NULL
    have <- if (is.null(mm)) integer(0) else
      intersect(unique(mm$system_id), sys_tab$system_id)
    lack <- setdiff(sys_tab$system_id, have)
    if (!length(have))
      stop("module = TRUE but no module map covers the chosen system(s) ",
           paste(lack, collapse = ", "),
           "; the bundled module map covers systems ",
           paste(sort(unique(mm$system_id)), collapse = ", "), call. = FALSE)
    if (length(lack)) {
      msgt <- paste("module level runs only for system(s)",
                    paste(have, collapse = ", "), "- no module map for",
                    paste(lack, collapse = ", "))
      if (!autofix) stop(msgt, call. = FALSE)
      cor <- .frame_note(cor, "module", "module map incomplete", msgt)
    }
    mods <- mm[mm$system_id %in% have, , drop = FALSE]
  }

  ## ---- covariates + autofixable derived columns ---------------------------
  if (identical(covariates, "contract")) {
    .bld <- .frame_which_build(names(data))
    covariates <- .frame_contract(names(data))
    if (!is.null(.bld)) {
      ct <- .frame_contracts[[.bld]]
      cor <- .frame_note(cor, "contract", ct$label,
                         paste0("fixed: ", paste(ct$fixed, collapse = ", "),
                                " | blocks: ", paste(ct$block, collapse = ", "),
                                if (is.null(ct$chip)) " | chip: NEVER included"
                                else paste0(" | chip: ", ct$chip)))
      ## adopt the build's blocks when the caller left the default
      if (identical(random_effects, "cID") &&
          !identical(ct$block, "cID") && all(ct$block %in% names(data))) {
        random_effects <- ct$block
        cor <- .frame_note(cor, "random_effects",
                           "repeated-measures build detected",
                           paste0("blocks set to ", paste(ct$block, collapse = ", "),
                                  " - permuting within cID alone would break the within-person pairing"))
      }
      ## the before/after builds forbid chip; honour that
      if (is.null(ct$chip)) attr(covariates, "no_chip") <- TRUE
    }
  }
  ## sex_c is DERIVED, not a column of any workbook. Deriving it only under
  ## covariates = "contract" meant a caller who typed the contract out by hand
  ## silently lost the sex adjustment - same intent, different model, different
  ## p-values. Derive it whenever it is asked for.
  if ("sex_c" %in% covariates && !"sex_c" %in% names(data) &&
      "sex" %in% names(data)) {
    sx <- suppressWarnings(as.numeric(data$sex))
    data$sex_c <- ifelse(sx == 1, -0.5, +0.5)
    cor <- .frame_note(cor, "sex_c", "absent",
                       "built from `sex` (1 -> -0.5, else +0.5)")
  }
  ## `chip` is an argument, not only an auto-detection: on the Alpha parent
  ## build chip_T1 has 61 levels on ~400 arrays (7 of them singletons), which
  ## spends 60 df on a nuisance factor that ctrlSV3/ctrlSV5 already partly
  ## absorb. Whether to pay that is a judgement the analyst has to be able to
  ## make and report, so chip = FALSE has to be sayable.
  chip_req <- chip
  chip <- NULL
  if (isTRUE(attr(covariates, "no_chip"))) {
    cor <- .frame_note(cor, "chip", "forbidden by this build's contract",
                       "not entered: the before/after builds list chip as NEVER include")
  } else if (identical(chip_req, FALSE)) {
    cor <- .frame_note(cor, "chip", "chip = FALSE",
                       paste0("not entered as a fixed effect - residual chip ",
                              "effects are left to ctrlSV3/ctrlSV5; say so in ",
                              "the write-up"))
  } else if (is.character(chip_req) && length(chip_req) == 1L && nzchar(chip_req)) {
    pool <- identical(chip_req, "pool") && !("pool" %in% names(data))
    src <- if (pool) {
      if ("chip_f" %in% names(data)) "chip_f" else
        grep("^chip(_T[0-9]+)?$", names(data), value = TRUE)[1]
    } else chip_req
    if (is.na(src) || !src %in% names(data))
      stop("chip = \"", chip_req, "\" is not a column of `data`", call. = FALSE)
    cf <- factor(data[[src]])
    if (pool) {
      ## A chip level holding a single array is fitted exactly by its own dummy:
      ## leverage 1, residual 0. That array is counted in n and contributes
      ## nothing to any estimate. On the Alpha parent build 7 of 61 chips are
      ## singletons, so 7 of 396 arrays are structurally uninformative while
      ## still being reported. Pooling the rare levels keeps the batch
      ## adjustment, recovers the degrees of freedom, and removes the
      ## leverage-1 rows.
      tb <- table(cf)
      rare <- names(tb)[tb < 3L]
      if (length(rare)) {
        lv <- as.character(cf); lv[lv %in% rare] <- "chip_small"
        cf <- factor(lv)
        cor <- .frame_note(cor, "chip_f", "chip = \"pool\"",
          sprintf(paste0("%d chip level(s) holding fewer than 3 arrays pooled ",
                         "into 'chip_small' (%d array(s)); %d level(s) remain, ",
                         "%d df instead of %d"),
                  length(rare), sum(tb[rare]), nlevels(cf), nlevels(cf) - 1L,
                  length(tb) - 1L))
      } else {
        cor <- .frame_note(cor, "chip_f", "chip = \"pool\"",
                           "no chip level held fewer than 3 arrays - nothing pooled")
      }
    } else {
      cor <- .frame_note(cor, "chip_f", paste0("chip = \"", chip_req, "\""),
                         paste0("built as factor(", src, "), fixed effect"))
    }
    data$chip_f <- cf; chip <- "chip_f"
  } else
  if ("chip_f" %in% names(data)) chip <- "chip_f"
  else {
    chip_col <- grep("^chip(_T[0-9]+)?$", names(data), value = TRUE)[1]
    if (!is.na(chip_col)) {
      data$chip_f <- factor(data[[chip_col]]); chip <- "chip_f"
      cor <- .frame_note(cor, "chip_f", "absent",
                         paste0("built as factor(", chip_col, ")"))
    }
  }

  ## The Alpha contract enters chip as a RANDOM intercept, `(1 | chip_T1)`, not
  ## as a fixed factor. Through 1.17.0 dmsa used a fixed factor: 61 levels and
  ## 60 df on the parent build, with 7 singleton chips fitted perfectly. Under
  ## `chip_effect = "random"` (the default, and what the contract specifies)
  ## chip is removed from the fixed design and handled by a quasi-demeaning GLS
  ## transform instead - a median of ~16 effective df on these data, and zero on
  ## the probes where the chip variance is estimated as nil.
  chip_random <- ""
  if (identical(chip_effect, "none")) {
    if (.nz(chip)) cor <- .frame_note(cor, "chip", "chip_effect = 'none'",
      "chip was NOT entered in any form")
    chip <- ""
  } else if (identical(chip_effect, "random") && .nz(chip)) {
    chip_random <- chip; chip <- ""
    cor <- .frame_note(cor, "chip", "entered as a RANDOM intercept",
      paste0("(1 | chip) per the Alpha covariate contract - fitted by REML ",
             "quasi-demeaning, not as a ", nlevels(factor(data$chip_f)),
             "-level fixed factor"))
  } else if (identical(chip_effect, "fixed") && .nz(chip)) {
    cor <- .frame_note(cor, "chip", "entered as a FIXED factor",
      paste0("deviation from the Alpha contract, which specifies (1 | chip); ",
             "costs ", nlevels(factor(data$chip_f)) - 1L, " df"))
  }
  miss <- setdiff(covariates, names(data))
  if (length(miss)) {
    if (!autofix)
      stop("covariate(s) not in `data`: ", paste(miss, collapse = ", "),
           call. = FALSE)
    cor <- .frame_note(cor, paste(miss, collapse = ","),
                       "covariate absent from data", "dropped from the model")
    ## a mistyped covariate name is indistinguishable from an intended one;
    ## dropping it silently produces an UNADJUSTED model that still looks
    ## adjusted in the write-up. Make it a warning, not a table footnote.
    warning("covariate(s) not found in `data` and dropped from the model: ",
            paste(miss, collapse = ", "),
            " - the reported model is NOT adjusted for them", call. = FALSE)
    covariates <- setdiff(covariates, miss)
  }

  ## A contract covariate derived from an outcome must never enter the model.
  ## Without this the frame builds sex_c from `sex`, then discovers the design
  ## is rank-deficient, then drops it - three messages that leave the user to
  ## work out for themselves that the covariate WAS the outcome.
  if (length(covariates)) {
    drop_oc <- character(0)
    for (oc0 in outcome) {
      hits <- grep(paste0("^", oc0, "(_c|_z|_f)?$"), covariates, value = TRUE)
      drop_oc <- c(drop_oc, hits)
    }
    drop_oc <- unique(setdiff(drop_oc, character(0)))
    if (length(drop_oc)) {
      covariates <- setdiff(covariates, drop_oc)
      cor <- .frame_note(cor, paste(drop_oc, collapse = ","),
                         "derived from an outcome column",
                         paste0("excluded before fitting: a covariate built ",
                                "from ", paste(outcome, collapse = "/"),
                                " cannot adjust a model that tests it"))
    }
  }

  ## ---- TEST DRIVE 1: numeric coercion with per-column loss ---------------
  numvars <- unique(c(outcome, covariates,
                      c(mod, mod2)[nzchar(c(mod, mod2))]))
  for (v in numvars) {
    old <- data[[v]]
    if (is.numeric(old)) next
    new <- suppressWarnings(as.numeric(as.character(old)))
    lost <- sum(is.na(new) & !is.na(old))
    if (!autofix && lost > 0)
      stop("column '", v, "' is non-numeric and coercion would lose ",
           lost, " values", call. = FALSE)
    data[[v]] <- new
    cor <- .frame_note(cor, v, "non-numeric",
                       sprintf("coerced to numeric%s",
                               if (lost > 0) sprintf(" (%d values -> NA)", lost)
                               else ""))
  }

  ## complete cases over model variables + blocks + methylation columns
  blockv <- intersect(random_effects, names(data))
  if (length(blockv) < length(random_effects)) {
    gone <- setdiff(random_effects, blockv)
    if (!autofix)
      stop("random_effects column(s) missing: ", paste(gone, collapse = ", "),
           call. = FALSE)
    cor <- .frame_note(cor, paste(gone, collapse = ","),
                       "random-effects column absent",
                       "dropped; permutation is unrestricted unless another block remains")
  }
  ## chip_random must be RETAINED even though it is not in the fixed design:
  ## it is the grouping vector for the random intercept. Dropping it here made
  ## ri_group NULL and the random intercept silently never applied.
  needcols <- unique(c(numvars, blockv, if (.nz(chip)) chip,
                       if (.nz(chip_random)) chip_random))
  cc <- stats::complete.cases(data[, needcols, drop = FALSE])
  loss <- vapply(needcols, function(v) sum(is.na(data[[v]])), integer(1))
  worst <- needcols[which.max(loss)]
  if (any(loss > 0))
    cor <- .frame_note(cor, "complete cases",
                       sprintf("%d of %d rows dropped", sum(!cc), length(cc)),
                       sprintf("costliest column: %s (%d missing)",
                               worst, max(loss)))
  base <- data[cc, , drop = FALSE]
  if (nrow(base) < 30)
    stop("only ", nrow(base), " complete rows remain - too few", call. = FALSE)
  if (.nz(chip)) base[[chip]] <- droplevels(factor(base[[chip]]))
  for (v in covariates) if (is.factor(base[[v]]))
    base[[v]] <- droplevels(base[[v]])

  ## methylation matrix on the complete-case rows
  if (is.matrix(methylation)) {
    key <- if (all(mp$column %in% colnames(methylation))) mp$column else mp$probe
    BE <- methylation[cc, match(key, colnames(methylation)), drop = FALSE]
  } else BE <- as.matrix(base[, mp$column, drop = FALSE])
  mode(BE) <- "numeric"
  okc <- apply(BE, 2, function(z) all(is.finite(z)))
  if (any(!okc)) {
    cor <- .frame_note(cor, "methylation",
                       sprintf("%d probe column(s) with missing values",
                               sum(!okc)), "dropped from the frame")
    mp <- mp[okc, , drop = FALSE]; BE <- BE[, okc, drop = FALSE]
  }
  ## Dropping EVERY probe column is a different situation from dropping some,
  ## and it has one usual cause: the methylation set being analysed was not
  ## measured on the same people as the covariates. Saying so beats letting the
  ## frame continue with an empty map and fail later inside the aligner.
  if (!nrow(mp))
    stop("every mapped probe column had missing values on the ", nrow(base),
         " complete-case row(s), so no probe is left to analyse. This usually ",
         "means the methylation set and the covariates were measured on ",
         "different people - e.g. analysing one family member's array with ",
         "another's cell-composition and batch covariates. Check that the ",
         "covariates belong to the same array as `methylation`.",
         call. = FALSE)
  mv <- .frame_mvalues(BE)
  if (mv$converted)
    cor <- .frame_note(cor, "methylation", "beta values detected",
                       "converted to M-values (log2 b/(1-b), eps 1e-4)")
  M <- mv$M

  ## ---- TEST DRIVE 2: design rank audit ------------------------------------
  ff <- stats::as.formula(paste("~", paste(
    c(outcome, covariates, if (.nz(chip)) chip,
      c(mod, mod2)[nzchar(c(mod, mod2))]), collapse = " + ")))
  X <- stats::model.matrix(ff, base)
  nearconst <- colnames(X)[-1][apply(X[, -1, drop = FALSE], 2,
                                     function(z) stats::sd(z) < 1e-10)]
  qr_ <- qr(X)
  if (qr_$rank < ncol(X)) {
    aliased <- colnames(X)[qr_$pivot[(qr_$rank + 1):ncol(X)]]
    keepvars <- covariates
    for (a in aliased) {
      hit <- covariates[vapply(covariates, function(v) grepl(v, a, fixed = TRUE),
                               logical(1))]
      if (length(hit)) keepvars <- setdiff(keepvars, hit)
    }
    ## compositional note: does a dropped column track the sum of the others?
    comp <- ""
    Z <- X[, -1, drop = FALSE]
    for (a in intersect(aliased, colnames(Z))) {
      others <- setdiff(colnames(Z), a)
      s <- rowSums(Z[, others, drop = FALSE])
      cc_ <- if (stats::sd(s) > 0 && stats::sd(Z[, a]) > 0)
        suppressWarnings(stats::cor(Z[, a], s)) else NA_real_
      if (is.finite(cc_) && abs(cc_) > .999)
        comp <- " (compositional group: it mirrors the sum of the others - e.g. cell fractions)"
    }
    if (!autofix)
      stop("design is rank-deficient; aliased: ",
           paste(aliased, collapse = ", "), call. = FALSE)
    cor <- .frame_note(cor, paste(aliased, collapse = ","),
                       "aliased (rank-deficient design)",
                       paste0("dropped from the model", comp))
    covariates <- keepvars
  }
  if (length(nearconst)) {
    if (!autofix) stop("near-constant column(s): ",
                       paste(nearconst, collapse = ", "), call. = FALSE)
    hit <- covariates[vapply(covariates, function(v)
      any(grepl(v, nearconst, fixed = TRUE)), logical(1))]
    if (length(hit)) {
      covariates <- setdiff(covariates, hit)
      cor <- .frame_note(cor, paste(hit, collapse = ","), "near-constant",
                         "dropped from the model")
    }
  }
  Xs <- scale(X[, -1, drop = FALSE])
  Xs <- Xs[, apply(Xs, 2, function(z) all(is.finite(z))), drop = FALSE]
  kappa_ <- if (ncol(Xs)) kappa(crossprod(Xs) / (nrow(Xs) - 1)) else 1

  ## ---- TEST DRIVE 3: random-effects / block audit -------------------------
  blk <- if (length(blockv)) interaction(base[, blockv, drop = FALSE],
                                         drop = TRUE) else
    factor(seq_len(nrow(base)))
  bt <- table(table(blk))
  sizes <- as.integer(names(bt))
  ## permutation count within equal-size strata (couples + singletons lesson)
  logperm <- 0
  for (i in seq_along(bt)) logperm <- logperm + lfactorial(bt[i]) / log(2)
  singles <- if ("1" %in% names(bt)) bt[["1"]] else 0
  if (all(sizes == 1))
    cor <- .frame_note(cor, paste(blockv, collapse = ","),
                       "all blocks are singletons",
                       "reduces to unrestricted permutation")
  lost_strata <- sum(bt == 1 & sizes > 1)
  if (lost_strata > 0)
    cor <- .frame_note(cor, paste(blockv, collapse = ","),
                       sprintf("%d stratum/strata hold a single block", lost_strata),
                       "those blocks cannot move; counted as lost exchangeable units")
  ## block confounded with a fixed covariate?
  for (v in intersect(covariates, names(base))) {
    if (!is.numeric(base[[v]])) next
    wv <- tapply(base[[v]], blk, function(z) stats::var(z))
    if (all(!is.finite(wv) | wv < 1e-12) && stats::sd(base[[v]]) > 0 &&
        nlevels(blk) < nrow(base))
      cor <- .frame_note(cor, v, "constant within every block",
                         "block and covariate are confounded; permutation respects both")
  }
  if (logperm < 10)
    stop(sprintf(paste0("only ~2^%.1f distinct block permutations are ",
                        "available - too few for stable p-values. Merge ",
                        "blocks or drop the random term."), logperm),
         call. = FALSE)

  ## ---- TEST DRIVE 4: outcome audit ----------------------------------------
  for (oc in outcome) {
    u <- unique(base[[oc]][is.finite(base[[oc]])])
    if (length(u) <= 1) stop("outcome '", oc, "' is constant", call. = FALSE)
    if (length(u) == 2 && outcome_type == "gaussian") {
      if (autofix) {
        outcome_type <- "logistic"
        ## "switched to logistic" reads as though a logistic regression is
        ## fitted. None is. Methylation is the RESPONSE and the outcome is a
        ## predictor of it, so a two-level outcome is a group contrast on a
        ## continuous response - and, in the before/after builds, a
        ## within-person one.
        .k <- .frame_outcome_kind(base, oc, blockv,
                                  .frame_which_build(names(data)))
        cor <- .frame_note(cor, oc, "two-level outcome",
          if (identical(.k$kind, "wave"))
            paste0("treated as a WITHIN-PERSON contrast between the two waves ",
                   "(", oc, " = ", .k$lo, " vs ", .k$hi, "); methylation is the ",
                   "response, so no logistic model is fitted")
          else
            paste0("treated as a group contrast (", oc, " = ", .k$lo, " vs ",
                   .k$hi, ") on a continuous response; no logistic model is ",
                   "fitted"))
      } else stop("outcome '", oc, "' is binary; set outcome_type",
                  call. = FALSE)
    }
    if (length(u) > 2 && length(u) <= 6 && outcome_type == "gaussian" &&
        all(u == round(u)))
      cor <- .frame_note(cor, oc,
                         sprintf("%d distinct integer levels", length(u)),
                         "kept gaussian and analysed as an ORDERED score - if the levels are unordered this is wrong; recode to a binary contrast or use dmsa_outcome()")
  }

  ## ---- alignment per system under the chosen map --------------------------
  sets <- list()
  for (sid in sys_tab$system_id) {
    s <- mp[mp$system_id == sid, , drop = FALSE]
    ## A selected system with no surviving row used to reach dmsa_align() with a
    ## zero-row table and die inside it on "replacement has 1 row, data has 0",
    ## which says nothing about what the user did wrong.
    if (!nrow(s))
      stop("system '", sys_tab$system[match(sid, sys_tab$system_id)][1],
           "' (id ", sid, ") has no direction-called probe among the ",
           "methylation columns being analysed. ",
           nrow(mp), " mapped probe(s) remain overall, covering system id(s) ",
           paste(unique(mp$system_id), collapse = ", "),
           ". If `methylation` was given explicitly, it may not contain this ",
           "system's columns.", call. = FALSE)
    al <- dmsa_align(data.frame(cpg = s$probe, d = s$best_direction,
                                p_plus = s$p_plus),
                     genes = s$gene, level = "gene")
    sets[[as.character(sid)]] <-
      list(system_id = sid, system = s$system[1], map = s, alignment = al,
           columns = s$column, genes = unique(s$gene))
  }

  ## ---- map disagreement check (both maps, always) -------------------------
  conflicts <- .frame_conflicts(maps, sys_tab$system_id)
  if (nrow(conflicts))
    warning("direction map disagreement (confidence vs full) for: ",
            paste(sprintf("%s '%s'", conflicts$level, conflicts$unit),
                  collapse = "; "),
            " - analysed under cpg_map = '", cpg_map,
            "'; see frame$map_conflicts", call. = FALSE)

  shared_pairs <- tryCatch(.frame_shared_probes(mp, cs),
                           error = function(e) NULL)
  if (!is.null(shared_pairs) && nrow(shared_pairs))
    cor <- .frame_note(cor, "shared probes",
      sprintf("%d gene pair(s) share probe columns inside a system",
              nrow(shared_pairs)),
      paste0("the same measurement appears under both names, so they are NOT ",
             "independent units in their family: ",
             paste(sprintf("%s (%d)", shared_pairs$pair, shared_pairs$n_shared),
                   collapse = "; ")))

  frame <- structure(list(
    data = base, M = M, map = mp, sets = sets, systems = sys_tab,
    modules = mods, selection = cs,
    polarity_table = if (is.null(cs)) NULL else cs$polarity,
    module_evidence = if (is.null(cs)) NULL else
      cs$modules[as.character(cs$modules$system_id) %in%
                   as.character(sys_tab$system_id), , drop = FALSE],
    outcome = outcome, covariates = covariates, chip = chip,
    chip_random = chip_random, chip_effect = chip_effect,
    ## names for the two levels of a two-level outcome, e.g. c("T1", "T4") or
    ## c("female", "male"). Without them the sentence has to fall back to
    ## "time = 0" / "time = 1", which is unambiguous but ugly.
    outcome_levels = if (is.null(outcome_levels)) NULL else {
      .ol <- as.character(outcome_levels)
      if (length(.ol) != 2L)
        stop("`outcome_levels` must be exactly two names for the two levels ",
             "of a two-level outcome", call. = FALSE)
      .ol
    },
    outcome_kind = stats::setNames(lapply(outcome, function(.o)
      .frame_outcome_kind(base, .o, blockv, .frame_which_build(names(data)))),
      outcome),
    shared_probes = shared_pairs,
    labels = c(.lab_pairs(outcome, outcome_label, "outcome_label"),
               .lab_pairs(covariates, predictor_labels, "predictor_labels"),
               if (nzchar(mod) && !is.null(mod_label))
                 stats::setNames(as.character(mod_label)[1], mod) else NULL,
               if (nzchar(mod2) && !is.null(mod2_label))
                 stats::setNames(as.character(mod2_label)[1], mod2) else NULL),
    block = blk, block_cols = blockv,
    levels = c(system = isTRUE(system), module = isTRUE(module),
               gene = isTRUE(gene), probe = isTRUE(probe)),
    moderation = isTRUE(moderation), mod = mod, mod2 = mod2,
    type = type, outcome_type = outcome_type, frame_role = frame_role,
    cpg_map = cpg_map, B = as.integer(B), alpha = alpha, seed = seed,
    correction = correction, weighting = weighting, w_floor = w_floor,
    palette = palette,
    plots = plots, tables = tables, summary = summary,
    plot_type = plot_type, table_type = table_type, outdir = outdir,
    gene_models = gene_models, progress = progress, beep = beep,
    corrections = cor, map_conflicts = conflicts, untestable = untestable,
    kappa = kappa_, log2_permutations = as.numeric(logperm),
    block_sizes = bt, pilot = NULL), class = "dmsa_frame")

  ## ---- TEST DRIVE 5: B = 49 pilot on the smallest system ------------------
  smallest <- names(sets)[which.min(vapply(sets, function(s)
    length(s$columns), integer(1)))]
  t0 <- Sys.time()
  pr <- tryCatch(
    .report_gene_level(frame, outcome[1], smallest, B = 49L, quiet = TRUE),
    error = function(e) e)
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(pr, "error"))
    stop("pilot run (B = 49, system ", smallest, ") failed: ",
         conditionMessage(pr), call. = FALSE)
  if (!all(is.finite(pr$p_omnibus[pr$n_probes > 0])))
    stop("pilot run produced non-finite p-values; inspect the frame",
         call. = FALSE)
  units_total <- sum(vapply(sets, function(s) length(unique(s$map$gene)),
                            integer(1)))
  eta_min <- dt / 49 * B * length(outcome) *
    max(1, length(sets)) / max(1, 1) / 60
  frame$pilot <- list(seconds_b49 = dt, eta_minutes = round(eta_min, 1),
                      system = smallest, ok = TRUE)
  frame
}

#' Re-run the frame's test drive
#'
#' @param frame a \code{dmsa_frame}.
#' @return the frame's corrections table, invisibly (printed first).
#' @examples
#' set.seed(1)
#' map <- data.frame(gene = "NR3C1", system_id = 1L, system = "HPA axis",
#'                   probe = c("cg01", "cg02"), column = c("cg01", "cg02"),
#'                   best_direction = c(-1, 1), p_plus = c(.1, .9))
#' dat <- data.frame(anx = rnorm(40), cov1 = rnorm(40), cID = rep(1:20, each = 2),
#'                   cg01 = plogis(rnorm(40)), cg02 = plogis(rnorm(40)))
#' fr <- dmsa_frame(dat, map = map, outcome = "anx", covariates = "cov1",
#'                  B = 49, seed = 1, outdir = tempfile())
#'
#' ## beta values were detected and converted; every such fix is on the record
#' dmsa_test_drive(fr)
#' @export
dmsa_test_drive <- function(frame) {
  stopifnot(inherits(frame, "dmsa_frame"))
  cat("Test drive -", nrow(frame$corrections), "recorded action(s)\n")
  if (nrow(frame$corrections)) print(frame$corrections, row.names = FALSE)
  cat(sprintf("design condition number: %.1f\n", frame$kappa))
  cat(sprintf("block permutations: ~2^%.0f\n", frame$log2_permutations))
  if (!is.null(frame$pilot))
    cat(sprintf("pilot (B = 49, system %s): %.1fs -> full run ETA ~%.0f min\n",
                frame$pilot$system, frame$pilot$seconds_b49,
                frame$pilot$eta_minutes))
  invisible(frame$corrections)
}

#' @export
print.dmsa_frame <- function(x, ...) {
  cat("DMSA frame -", nrow(x$data), "rows,", nrow(x$map),
      "mapped probes,", nrow(x$systems), "system(s)\n")
  lv <- names(x$levels)[x$levels]
  cat("  levels:   ", paste(lv, collapse = " > "), "\n")
  fam <- vapply(x$sets, function(s) length(unique(s$map$gene)), integer(1))
  cat("  families: ", paste(sprintf("%s (%d genes, %d probes)",
      vapply(x$sets, `[[`, "", "system"), fam,
      vapply(x$sets, function(s) length(s$columns), integer(1))),
      collapse = "; "), "\n")
  cat("  outcome(s):", paste(x$outcome, collapse = ", "),
      if (x$frame_role == "outcome") " [frame_role = outcome: these PREDICT the frame]" else "",
      "\n")
  cat("  model:    ", x$outcome_type, x$type, "| outcome enters as a PREDICTOR of methylation",
      if (x$outcome_type != "gaussian")
        "\n             (so a binary/categorical outcome is a group contrast on a continuous response - no link function is needed or used)" else "",
      "\n            ",
      if (x$moderation) sprintf("| moderation: frame x %s%s", x$mod,
        if (nzchar(x$mod2)) paste0(" x ", x$mod2) else "") else "",
      "\n")
  if (!is.null(x$module_evidence) && nrow(x$module_evidence)) {
    cat("  sets:     ", x$selection$name, "\n", sep = "")
    .cas_evidence_banner(x$module_evidence, indent = "  ")
    .cas_polarity_banner(list(polarity = x$polarity_table), indent = "  ")
  }
  cat("  map:      cpg_map =", x$cpg_map,
      if (nrow(x$map_conflicts)) sprintf("(%d unit(s) DISAGREE between maps - see $map_conflicts)",
                                         nrow(x$map_conflicts)) else "(maps agree)", "\n")
  cat("  engine:   ", (x$weighting %||% "combined"),
      switch(x$weighting %||% "combined",
        combined = "(flat + reliability fused per lens; joint-null maxT)",
        reliability = "(centrality-weighted)", flat = "(equal weights)"), "\n")
  cat("  inference: B =", x$B, "|", x$correction,
      "within level-local families | blocks:",
      if (length(x$block_cols)) paste(x$block_cols, collapse = ", ") else "none",
      sprintf("(~2^%.0f permutations)\n", x$log2_permutations))
  if (!is.null(x$pilot))
    cat(sprintf("  pilot:     B = 49 in %.1fs -> ETA ~%.0f min for the full report\n",
                x$pilot$seconds_b49, x$pilot$eta_minutes))
  if (nrow(x$corrections)) {
    cat("  corrections (test drive):\n")
    for (i in seq_len(nrow(x$corrections)))
      cat(sprintf("   - %s: %s -> %s\n", x$corrections$field[i],
                  x$corrections$issue[i], x$corrections$action[i]))
  }
  invisible(x)
}
