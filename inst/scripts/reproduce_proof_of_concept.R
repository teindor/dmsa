## Proof-of-concept analyses reported in the manuscript.
##
## Participant data are held under IRB approval (Reichman University 5-2020;
## Israel Ministry of Health Helsinki Committee) and are NOT distributed with
## this package. `alpha` below is the merged Project Alpha parent build,
## available from the corresponding author under a data-use agreement.
##
## chip_effect = "fixed" is pinned deliberately. The study covariate contract
## specifies a random intercept, `(1 | chip_T1)`, and from 1.18.0 that is the
## package default; the analyses reported in the manuscript were run with chip
## as a FIXED factor, as the Methods state. Pinning it here reproduces the
## published table exactly. Both findings survive under either specification.

library(dmsa)

fr <- dmsa_frame(
  data        = alpha,
  systems     = c("oxytocin", "imprinted"),
  outcome     = c("Attachment_Anxiety_General_T1",
                  "Attachment_Avoidance_General_T1"),
  covariates  = "contract",
  module      = TRUE,
  cpg_map     = "full",
  weighting   = "flat",
  chip_effect = "fixed",
  B           = 1999,
  seed        = 1,
  outcome_label = c("Attachment anxiety (T1)", "Attachment avoidance (T1)"),
  outdir      = "proof_of_concept")

rp <- dmsa_report(fr)

## Expected (manuscript Table 1):
##   AVP    coherence family-adjusted p = .0120   survives
##   PLAGL1 composite family-adjusted p = .0080   survives
##   PLAGL1 coherence family-adjusted p = .0425
##   system diffuse   family-adjusted p = .0095
##   OXTR   does not survive its family
stopifnot(nrow(fr$data) == 396L)
