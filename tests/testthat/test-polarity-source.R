## Spec 39/40/41: polarity comes from the ACTIVE reference; it is never
## silently replaced by all +1, and a malformed user polarity is never
## silently discarded.

.pol_cascade <- function(w_g = NULL) {
  cas <- data.frame(system_id = 1, system = "S", module_id = "1.1",
                    module = "M", gene = c("A", "B"), cpg = c("cg1", "cg2"),
                    stringsAsFactors = FALSE)
  if (!is.null(w_g)) cas$w_g <- w_g
  cas
}

test_that("spec 41: a malformed `polarity` argument hard-errors", {
  e <- expect_error(dmsa_sets(.pol_cascade(), polarity = data.frame(x = 1:2)))
  expect_match(conditionMessage(e), "could not be read")
  expect_match(conditionMessage(e), "weights every gene \\+1")
})

test_that("spec 41: a malformed w_g column in the cascade hard-errors", {
  e <- expect_error(dmsa_sets(.pol_cascade(w_g = c("up", "down"))))
  expect_match(conditionMessage(e), "could not be read as polarity")
})

test_that("a valid w_g column still becomes the cascade's polarity", {
  s <- dmsa_sets(.pol_cascade(w_g = c(1, -1)))
  expect_equal(nrow(s$polarity$polarity), 2L)
})

## self-contained sim: a user-supplied map with NO polarity anywhere
.pol_sim <- function(n = 120, seed = 42) {
  set.seed(seed)
  map <- expand.grid(gene = paste0("G", 1:6), k = 1:4,
                     stringsAsFactors = FALSE)
  map$system_id <- ifelse(map$gene %in% paste0("G", 1:3), 1L, 2L)
  map$system <- ifelse(map$system_id == 1, "Sim system one", "Sim system two")
  map$probe <- sprintf("cg%07d", seq_len(nrow(map)))
  map$column <- paste0(map$probe, "_", map$gene)
  map$best_direction <- rep(c(-1, 1), length.out = nrow(map))
  map$p_plus <- ifelse(map$best_direction > 0, .9, .1)
  map$best_tier <- "A"; map$smr_tier <- ""
  d <- data.frame(out1 = stats::rnorm(n), cov1 = stats::rnorm(n),
                  cov2 = stats::rnorm(n), cID = rep(seq_len(n / 2), each = 2))
  for (i in seq_len(nrow(map))) {
    sig <- if (map$gene[i] == "G1")
      0.35 * d$out1 * sign(map$best_direction[i]) else 0
    d[[map$column[i]]] <- stats::plogis(stats::rnorm(n, sd = 1) + sig)
  }
  list(data = d, map = map)
}

test_that("spec 39: the system level is SKIPPED, not weighted +1, with no polarity", {
  s <- .pol_sim()
  od <- file.path(tempdir(), paste0("dmsa_pol_", as.integer(runif(1, 1, 1e6))))
  fr <- dmsa_frame(s$data, map = s$map, outcome = "out1",
                   covariates = c("cov1", "cov2"), random_effects = "cID",
                   B = 99, seed = 1, outdir = od, plot_type = "png")
  expect_message(r <- dmsa_report(fr), "system level SKIPPED")
  ## the gene level is unaffected - it aligns on d alone and needs no w_g
  expect_s3_class(r, "dmsa_report")
  expect_true(any(r$results$level == "gene"))
  expect_false(any(r$results$level == "system"))
})
