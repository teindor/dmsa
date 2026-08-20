# ============================================================================
# ARTIFACT-PROBE MASK for EPIC probes carrying a reported finding.
#
# Any single-probe or two-probe finding on a methylation array has two standard
# artifact explanations that must be excluded before it is read biologically:
#
#   (1) a polymorphism under the probe, so the array reads genotype as
#       methylation - the decisive alternative whenever a finding could be
#       inherited rather than acquired;
#   (2) a probe that cross-hybridises elsewhere in the genome, so the signal
#       does not come from the annotated locus at all.
#
# This script checks a set of probes against five published lists covering both.
#
# SOURCES. Neither publisher exposes the supplements at a stable direct URL,
# but both are redistributed inside public repositories:
#
#   Pidsley et al. (2016) Genome Biology 17:208
#     MOESM1  cross-reactive EPIC probes            (43,254)
#     MOESM4  SNP at the target CpG                 (12,378 probes)
#     MOESM5  INDEL at the target CpG                  (413 probes)
#     MOESM6  1000G variant anywhere in probe body  (97,345 probes)
#     -> github.com/sirselim/illumina450k_filtering  (EPIC/)
#
#   McCartney et al. (2016) Genomics Data 9:22-24
#     mmc2 CpG-targeting cross-hybridising probes
#     mmc3 non-CpG-targeting cross-hybridising probes   (44,210 combined)
#     -> github.com/markgene/maxprobes                (inst/extdata/)
#
# USAGE
#   git clone --depth 1 https://github.com/sirselim/illumina450k_filtering.git
#   git clone --depth 1 https://github.com/markgene/maxprobes.git
#   Rscript artifact_probe_mask.R  <dir_filtering> <dir_maxprobes>  cg1 cg2 ...
#
# Reports one line per probe per criterion. "clean" means absent from that
# list; it does NOT mean the probe is above suspicion, only that the documented
# artifact classes do not apply. A variant absent from the 1000 Genomes
# reference cannot be excluded this way.
# ============================================================================
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: <illumina450k_filtering dir> <maxprobes dir> <probe> [probe ...]")
FDIR   <- file.path(args[1], "EPIC")
MDIR   <- file.path(args[2], "inst", "extdata")
probes <- args[-(1:2)]

need <- c(file.path(FDIR, c("13059_2016_1066_MOESM1_ESM.csv",
                            "13059_2016_1066_MOESM4_ESM.csv",
                            "13059_2016_1066_MOESM5_ESM.csv",
                            "13059_2016_1066_MOESM6_ESM.csv")),
          file.path(MDIR, c("1-s2.0-S221359601630071X-mmc2.txt",
                            "1-s2.0-S221359601630071X-mmc3.txt")))
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing list files:\n  ", paste(miss, collapse = "\n  "))

pid_cr <- rownames(read.csv(need[1], row.names = 1))
snp_at <- fread(need[2]); ind_at <- fread(need[3]); var_body <- fread(need[4])
mcc    <- unique(c(readLines(need[5]), readLines(need[6])))

cat(sprintf("lists: cross-reactive %d | cross-hybridising %d | SNP-at-target %d | INDEL-at-target %d | variant-in-body %d\n\n",
            length(pid_cr), length(mcc), uniqueN(snp_at$PROBE),
            uniqueN(ind_at$PROBE), uniqueN(var_body$PROBE)))

flag <- function(x, set) if (x %in% set) "FLAGGED" else "clean"
for (p in probes) {
  b <- var_body[PROBE == p]
  cat(sprintf("%s\n", p))
  cat(sprintf("   cross-reactive (Pidsley MOESM1)     : %s\n", flag(p, pid_cr)))
  cat(sprintf("   cross-hybridising (McCartney)       : %s\n", flag(p, mcc)))
  cat(sprintf("   SNP at target CpG (MOESM4)          : %s\n", flag(p, snp_at$PROBE)))
  cat(sprintf("   INDEL at target CpG (MOESM5)        : %s\n", flag(p, ind_at$PROBE)))
  if (nrow(b)) {
    cat(sprintf("   variant in probe body (MOESM6)      : FLAGGED - %d variant(s), max AF %.3f\n",
                nrow(b), max(b$AF, na.rm = TRUE)))
    print(b[order(-AF)][seq_len(min(5, nrow(b))),
                        .(VARIANT_ID, VARIANT_TYPE, REF, ALT, AF, EUR_AF)])
  } else cat("   variant in probe body (MOESM6)      : clean\n")
  cat("\n")
}
