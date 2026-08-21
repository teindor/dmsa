# dmsa design layer
# --------------------------------------------------------------------------
# The covariate contract and the dependence structure are declared ONCE, as
# an object, and validated against the data before anything is fitted. This
# exists because a hand-written covariate list silently omitted a required
# batch random effect across an entire analysis family; the contract must be
# machine-checked, not remembered.

#' Declare a DMSA design: covariates, dependence, and permutation blocks
#'
#' @param focal Character. The focal term(s) whose pooled aligned effect is
#'   tested. May be a main effect (\code{"IC_c"}), an interaction written with
#'   a colon (\code{"time:BSI_c"}), or several terms, in which case the first
#'   is the tested one and the rest are mutually-adjusted co-focal terms.
#' @param fixed Character vector of fixed-effect covariate columns. These are
#'   the build's contract, in full - no defaults are supplied, because a
#'   default is exactly what goes stale.
#' @param random Grouping factors for random intercepts: a character vector
#'   (\code{c("chip","cID")}, expanded to \code{(1|chip) + (1|cID)}) or a
#'   one-sided formula (\code{~ (1|chip) + (1|cID)}). \code{NULL} declares
#'   independent observations - which \code{dmsa_check_design()} will
#'   challenge if it finds repeated block IDs.
#' @param exchangeable Character. The column defining the outermost
#'   independent unit for permutation: couples share a family, repeated
#'   measures share a person, and a person sits inside a family - so the
#'   exchangeable unit is the family in all three cases. Rows inside a block
#'   move together and keep their order, which is what makes the same
#'   mechanism valid cross-sectionally and longitudinally.
#' @param forbid Character vector of columns that must NOT enter the model
#'   (the build's never-list). Supplying these turns a documentation rule
#'   into an error.
#' @param label Optional name for printing.
#' @return An object of class \code{dmsa_design}.
#' @examples
#' dmsa_design(focal = "IC_c",
#'             fixed = c("sex_c", "age", "Fib", "ctrlSV3", "ctrlSV5"),
#'             random = c("chip", "cID"), exchangeable = "cID",
#'             forbid = c("ctrlSV4", "plate"))
#' @export
dmsa_design <- function(focal, fixed, random = NULL, exchangeable = NULL,
                        forbid = character(), label = NULL) {
  if (!length(focal)) stop("focal must name at least one term", call. = FALSE)
  focal <- as.character(focal)
  fixed <- if (is.null(fixed)) character() else as.character(fixed)

  rg <- character()
  if (!is.null(random)) {
    if (inherits(random, "formula")) {
      tl <- attr(stats::terms(random), "term.labels")
      bars <- grep("\\|", tl, value = TRUE)
      if (!length(bars))
        stop("random formula must use (1|group) form, e.g. ~ (1|chip) + (1|cID)",
             call. = FALSE)
      rg <- trimws(sub(".*\\|", "", gsub("[()]", "", bars)))
    } else {
      rg <- as.character(random)
    }
  }

  ## every variable the design touches, with interaction terms split out
  focal_vars <- unique(unlist(strsplit(focal, ":", fixed = TRUE)))
  vars <- unique(c(focal_vars, fixed, rg, exchangeable))

  bad <- intersect(forbid, c(focal_vars, fixed))
  if (length(bad))
    stop("these terms are on the never-list for this design but appear as ",
         "focal or fixed: ", paste(bad, collapse = ", "), call. = FALSE)

  structure(list(
    focal = focal, focal_test = focal[1], focal_vars = focal_vars,
    fixed = fixed, random_groups = rg, exchangeable = exchangeable,
    forbid = as.character(forbid), label = label,
    vars = vars
  ), class = "dmsa_design")
}

#' @export
print.dmsa_design <- function(x, ...) {
  cat("dmsa design", if (!is.null(x$label)) paste0(": ", x$label), "\n", sep = "")
  cat("  focal      ", paste(x$focal, collapse = " + "),
      if (length(x$focal) > 1) paste0("   (tested: ", x$focal_test, ")"), "\n", sep = "")
  cat("  fixed      ", if (length(x$fixed)) paste(x$fixed, collapse = ", ") else "(none)", "\n")
  cat("  random     ", if (length(x$random_groups))
        paste0("(1|", x$random_groups, ")", collapse = " + ") else "(none declared)", "\n")
  cat("  exchange   ", if (is.null(x$exchangeable)) "(none - permutation will assume independence)"
                       else x$exchangeable, "\n")
  if (length(x$forbid))
    cat("  never      ", paste(x$forbid, collapse = ", "), "\n")
  invisible(x)
}

#' Project Alpha build contracts
#'
#' The four 2026 Project Alpha builds, transcribed from the student pack's
#' \code{covariate_sets.csv} including its never-list. Supply \code{focal} for
#' the analysis at hand; everything else is fixed by the build.
#'
#' Build 2 and 3 are longitudinal: their contract asks for the
#' \code{time x predictor} interaction, so pass e.g.
#' \code{focal = "time:BSI_Total_c"}.
#'
#' @param build 1 (parents T1), 2 (T1 to T4), 3 (October 7), or 4 (children).
#' @param focal Character focal term(s), passed to \code{dmsa_design()}.
#' @param drop Optional character vector of contract covariates to drop, for
#'   a declared deviation (e.g. dropping \code{Epi_T1} when the exposure is
#'   the immune-cell fraction, its compositional complement). Recorded in the
#'   returned object so it shows up in print and cannot be silent.
#' @return A \code{dmsa_design}.
#' @examples
#' d <- alpha_design(1, focal = "IC_T1_c")
#' d
#' # builds 2 and 3 are longitudinal: the contract asks for the time interaction
#' alpha_design(2, focal = "time:BSI_Total_c")
#' # a deviation from the contract is recorded on the object, never silent
#' alpha_design(1, focal = "IC_T1_c", drop = "Epi_T1")$dropped
#' @export
alpha_design <- function(build, focal, drop = character()) {
  build <- as.integer(build)
  spec <- switch(as.character(build),
    "1" = list(
      fixed  = c("sex_c", "age_at_array_T1", "Epi_T1", "Fib_T1", "ctrlSV3_T1", "ctrlSV5_T1"),
      random = c("chip_T1", "cID"), exch = "cID",
      forbid = c("ctrlSV1_T1", "ctrlSV2_T1", "ctrlSV4_T1", "plate_T1",
                 "submission_T1", "submission_kit_T1", "chip_row_T1"),
      label  = "Alpha build 1 (parents T1)"),
    "2" = list(
      fixed  = c("sex_c", "Epi_array", "Fib_array", "ctrlSV5_array"),
      random = c("ID", "cID"), exch = "cID",
      forbid = c("submission_array", "plate_array", "chip_array", "ctrlSV1_array",
                 "ctrlSV2_array", "ctrlSV3_array", "ctrlSV4_array",
                 "age_at_array_array", "chip_row_array"),
      label  = "Alpha build 2 (T1 to T4)"),
    "3" = list(
      fixed  = c("sex_c", "Epi_array", "Fib_array", "ctrlSV5_array"),
      random = c("ID", "cID"), exch = "cID",
      forbid = c("submission_array", "plate_array", "chip_array", "ctrlSV1_array",
                 "ctrlSV2_array", "ctrlSV3_array", "ctrlSV4_array",
                 "age_at_array_array", "chip_row_array", "wave"),
      label  = "Alpha build 3 (October 7)"),
    "4" = list(
      fixed  = c("sex_c", "birth_week", "birth_weight", "Epi_T4"),
      random = "chip_T4", exch = "cID",
      forbid = c("Fib_T4", "age_at_array_T4", "ctrlSV1_T4", "ctrlSV2_T4",
                 "ctrlSV3_T4", "ctrlSV4_T4", "ctrlSV5_T4", "plate_T4",
                 "submission_T4", "submission_kit_T4", "chip_row_T4", "preterm"),
      label  = "Alpha build 4 (children)"),
    stop("build must be 1, 2, 3 or 4", call. = FALSE))

  keep <- setdiff(spec$fixed, drop)
  dropped <- intersect(spec$fixed, drop)
  out <- dmsa_design(focal = focal, fixed = keep, random = spec$random,
                     exchangeable = spec$exch, forbid = spec$forbid,
                     label = paste0(spec$label,
                       if (length(dropped)) paste0(" [declared deviation: -",
                         paste(dropped, collapse = ", "), "]") else ""))
  out$dropped <- dropped
  out
}

#' Validate a design against the data before fitting
#'
#' Checks the things that silently ruin a methylation set analysis: a
#' covariate that is constant in this subsample, a covariate that proxies the
#' focal term, a grouping factor with too few levels to estimate, a grouping
#' factor aliased with the focal term, and declared independence in data that
#' plainly repeats.
#'
#' @param design A \code{dmsa_design}.
#' @param data A data.frame.
#' @param strict If TRUE (default) problems are errors; if FALSE they are
#'   warnings and the checks are still returned.
#' @return Invisibly, a list with \code{n}, \code{problems}, \code{notes} and
#'   the per-term diagnostics.
#' @examples
#' set.seed(1)
#' d <- data.frame(exposure = rnorm(80), sex_c = rep(c(-0.5, 0.5), 40),
#'                 age = rnorm(80, 40, 5), Fib = 0,
#'                 chip = rep(1:10, each = 8), cID = rep(1:40, each = 2))
#' des <- dmsa_design("exposure", c("sex_c", "age"),
#'                    random = c("chip", "cID"), exchangeable = "cID")
#' dmsa_check_design(des, d)$n
#' # Fib is constant in this subsample: absorbed silently by lm, caught here
#' bad <- dmsa_design("exposure", c("sex_c", "Fib"), random = "chip",
#'                    exchangeable = "cID")
#' dmsa_check_design(bad, d, strict = FALSE)$problems
#' @export
dmsa_check_design <- function(design, data, strict = TRUE) {
  stopifnot(inherits(design, "dmsa_design"))
  data <- as.data.frame(data)
  problems <- character(); notes <- character()

  miss <- setdiff(design$vars, names(data))
  if (length(miss))
    stop("design refers to columns absent from data: ", paste(miss, collapse = ", "),
         call. = FALSE)

  hit <- intersect(design$forbid, names(data))
  used <- intersect(hit, c(design$focal_vars, design$fixed))
  if (length(used))
    problems <- c(problems, paste0("never-list term in the model: ",
                                   paste(used, collapse = ", ")))

  cc  <- stats::complete.cases(data[, design$vars, drop = FALSE])
  dat <- data[cc, , drop = FALSE]
  n <- nrow(dat)
  if (n < 10) problems <- c(problems, sprintf("only %d complete cases", n))

  ## constant covariates -- a constant cannot be a covariate
  const <- character()
  for (v in c(design$fixed, design$focal_vars)) {
    x <- dat[[v]]
    if (is.numeric(x)) { if (isTRUE(all.equal(stats::sd(x), 0))) const <- c(const, v) }
    else if (length(unique(x)) < 2) const <- c(const, v)
  }
  if (length(const))
    problems <- c(problems, paste0("constant in this subsample (drop it): ",
                                   paste(const, collapse = ", ")))

  ## collinearity among the fitted fixed part
  vif <- NULL
  num <- setdiff(c(design$focal_vars, design$fixed), const)
  num <- num[vapply(num, function(v) is.numeric(dat[[v]]), logical(1))]
  if (length(num) >= 2) {
    vif <- vapply(num, function(v) {
      rhs <- setdiff(num, v)
      r2 <- summary(stats::lm(stats::reformulate(rhs, response = v), data = dat))$r.squared
      1 / max(1 - r2, 1e-12)
    }, numeric(1))
    bad <- names(vif)[vif > 5]
    if (length(bad))
      problems <- c(problems, sprintf("VIF > 5: %s",
        paste(sprintf("%s=%.1f", bad, vif[bad]), collapse = ", ")))
    warn <- names(vif)[vif > 2.5 & vif <= 5]
    if (length(warn))
      notes <- c(notes, sprintf("VIF 2.5-5 (watch): %s",
        paste(sprintf("%s=%.1f", warn, vif[warn]), collapse = ", ")))
  }

  ## random-effect grouping factors
  grp <- list()
  for (g in design$random_groups) {
    lv <- length(unique(dat[[g]]))
    grp[[g]] <- c(levels = lv, mean_size = n / lv)
    if (lv < 3)
      problems <- c(problems, sprintf("(1|%s) has %d level(s) - not estimable", g, lv))
    else if (lv >= n)
      problems <- c(problems, sprintf("(1|%s) has one level per row - not estimable", g))
    ## aliasing: is the focal term constant within every level?
    ft <- design$focal_vars[1]
    if (is.numeric(dat[[ft]])) {
      wsd <- stats::aggregate(dat[[ft]], list(dat[[g]]), function(z) stats::sd(z))
      if (all(is.na(wsd$x) | wsd$x < 1e-9))
        problems <- c(problems, sprintf(
          "focal '%s' is constant within every level of %s - the two are aliased", ft, g))
    }
  }

  ## declared independence in obviously repeated data
  if (!length(design$random_groups) && !is.null(design$exchangeable)) {
    r <- max(table(dat[[design$exchangeable]]))
    if (r > 1)
      problems <- c(problems, sprintf(
        "no random effects declared, but %s has up to %d rows per level",
        design$exchangeable, r))
  }

  ## permutation blocks
  blocks <- NULL
  if (!is.null(design$exchangeable)) {
    tb <- table(dat[[design$exchangeable]])
    blocks <- c(n_blocks = length(tb), max_size = max(tb))
    sizes <- table(tb)
    lone <- sum(tb == 1)
    if (length(tb) < 10)
      problems <- c(problems, sprintf("only %d exchangeable blocks - permutation has no resolution",
                                      length(tb)))
    if (any(as.integer(names(sizes)) > 1 & sizes == 1))
      notes <- c(notes, "some block sizes occur once; those blocks cannot be swapped")
    if (lone) notes <- c(notes, sprintf("%d singleton block(s)", lone))
  } else {
    notes <- c(notes, "no exchangeable unit declared: permutation will treat rows as independent")
  }

  res <- list(n = n, problems = problems, notes = notes,
              vif = vif, groups = grp, blocks = blocks, constant = const)
  if (length(problems)) {
    msg <- paste0("design check failed:\n  - ", paste(problems, collapse = "\n  - "))
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  if (length(notes)) message("dmsa design notes:\n  - ", paste(notes, collapse = "\n  - "))
  invisible(res)
}
