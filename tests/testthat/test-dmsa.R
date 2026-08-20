test_that("NR3C1 chain: methylation up -> GR down -> HPA tone UP", {
  # Exposure increases NR3C1 promoter methylation: b = +0.3
  # cpgdirection: promoter CpG silences GR -> d = -1 (conf .8)
  # polarity: GR is the HPA brake -> w_g = -1
  # so s = d * w_g = +1: the probe VOTES FOR higher axis tone
  dir <- data.frame(cpg = "cg_nr3c1", d = -1, p_plus = 0.1)  # P(d=+1)=.1
  al  <- dmsa_align(dir, genes = "NR3C1", level = "system",
                    polarity = data.frame(gene = "NR3C1", w_g = -1))
  expect_equal(al$s, +1)
  expect_equal(al$p_s_plus, 0.9)          # 0.1*0 + 0.9*1
  tt <- dmsa_test(b = 0.3, se = 0.1, alignment = al, method = "fixed")
  expect_gt(tt$estimate, 0)               # contributes POSITIVE axis-tone evidence
})

test_that("gene level ignores polarity; system level applies it", {
  dir <- data.frame(cpg = c("p1", "p2"), d = c(-1, -1), p_plus = c(0.2, 0.2))
  g   <- c("CRH", "NR3C1")
  pol <- data.frame(gene = c("CRH", "NR3C1"), w_g = c(1, -1))
  gene_al <- dmsa_align(dir, g, level = "gene")
  sys_al  <- dmsa_align(dir, g, level = "system", polarity = pol)
  expect_equal(gene_al$s, c(-1, -1))      # both: expression down
  expect_equal(sys_al$s,  c(-1, +1))      # CRH down = tone down; GR down = tone UP
})

test_that("missing polarity errors and names the gene (ask-the-user contract)", {
  dir <- data.frame(cpg = "p1", d = 1, p_plus = 0.9)
  expect_error(
    dmsa_align(dir, genes = "MYSTERY1", level = "system",
               polarity = data.frame(gene = "OTHER", w_g = 1)),
    "MYSTERY1")
})

test_that("missing polarity 'zero' keeps probe but marks off-axis", {
  dir <- data.frame(cpg = c("p1","p2"), d = c(1,-1), p_plus = c(.9,.2))
  al <- suppressWarnings(
    dmsa_align(dir, genes = c("KNOWN","UNKNOWN"), level = "system",
               polarity = data.frame(gene = "KNOWN", w_g = 1),
               missing_polarity = "zero"))
  expect_equal(al$usable, c(TRUE, FALSE))
  expect_equal(al$reason[2], "off_axis_gene")
})

test_that("expected-sign shrinks uncertain probes toward zero contribution", {
  dir <- data.frame(cpg = c("sure","unsure"), d = c(1, 1), p_plus = c(0.99, 0.55))
  al  <- dmsa_align(dir, genes = c("A","B"), level = "system",
                    polarity = data.frame(gene = c("A","B"), w_g = c(1,1)))
  tt  <- dmsa_test(b = c(0.5, 0.5), se = c(0.1, 0.1), al, method = "expected")
  m   <- tt$table$multiplier
  expect_gt(m[1], 0.9); expect_lt(m[2], 0.2)   # 2p-1: .98 vs .10
})

test_that("abstained probe (d = NA) is excluded and reported", {
  dir <- data.frame(cpg = c("p1","p2"), d = c(1, NA), p_plus = c(.9, NA))
  al <- dmsa_align(dir, genes = c("A","A"), level = "system",
                   polarity = data.frame(gene = "A", w_g = 1))
  expect_equal(al$reason[2], "no_direction_call")
  tt <- dmsa_test(b = c(.4, .4), se = c(.1, .1), al, method = "fixed")
  expect_equal(tt$n_used, 1)
})

test_that("balance table catches one-sided coverage", {
  dir <- data.frame(cpg = paste0("p", 1:4), d = c(1,-1,1,NA), p_plus = c(.9,.1,.8,NA))
  g   <- c("ACT1","ACT2","OFF","BRAKE")
  pol <- data.frame(gene = g, w_g = c(1, 1, 0, -1))
  al  <- dmsa_align(dir, g, level = "system", polarity = pol)
  bal <- dmsa_balance(al)
  expect_equal(bal$activation_probes, 2)
  expect_equal(bal$brake_probes, 0)       # the brake gene's probe abstained
  expect_equal(bal$off_axis, 1)
})

test_that("permutation p has the +1 floor", {
  expect_equal(dmsa_perm_pvalue(10, rep(0, 499)), 1/500)
  expect_gt(dmsa_perm_pvalue(0.1, rnorm(499)), 0.5)
})

test_that("bundled Alpha tables load", {
  gs <- alpha_gene_systems()
  expect_true(all(c("system_id","system","gene") %in% names(gs)))
  expect_equal(length(unique(gs$system_id)), 30)
  pol <- alpha_polarity()
  expect_true(all(c("gene","w_g") %in% names(pol)))
  expect_equal(pol$w_g[pol$gene == "NR3C1" & pol$system_id == 2], -1)
  expect_equal(pol$w_g[pol$gene == "FKBP5" & pol$system_id == 2], +1)  # brake-of-brake
  expect_equal(pol$w_g[pol$gene == "AVPR2" & pol$system_id == 1], 0)   # specificity control
})

test_that("abstained probe cannot contribute via p_plus leakage (expected mode)", {
  dir <- data.frame(cpg = c("p1","p2"), d = c(1, NA), p_plus = c(.9, .9))
  al <- dmsa_align(dir, genes = c("A","A"), level = "system",
                   polarity = data.frame(gene = "A", w_g = 1))
  tt <- dmsa_test(b = c(.4, .4), se = c(.1, .1), al, method = "expected")
  expect_equal(tt$n_used, 1)
})
