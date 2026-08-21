# Load a gene-to-system polarity table

Load a gene-to-system polarity table

## Usage

``` r
dmsa_polarity(x = "alpha", sets = NULL)
```

## Arguments

- x:

  `"alpha"` for the bundled Project Alpha 2026c polarity table, a path
  to a CSV in the same shape, a `data.frame`, or a `dmsa_sets` /
  `dmsa_selection` / `dmsa_frame` object whose attached polarity should
  be used - so a narrowed selection reviews only its own genes. `NULL`
  returns `NULL`, which is how a caller opts out of polarity entirely
  (system-level scores then weight every gene +1 and say so).

- sets:

  Optional `dmsa_sets` or `dmsa_selection`. When given, coverage is
  reported against it: which of its genes have a sign, which do not, and
  at what evidence grade.

## Value

An object of class `dmsa_polarity`, or `NULL`.

## Details

Signs are read against the system's declared activation tone, not
against a gene's immediate partner. A gene that inhibits a brake
therefore carries `w_g = +1`: FKBP5 restrains GR, GR mediates HPA
negative feedback, so more FKBP5 means more axis drive. The `role`
vocabulary names this explicitly (`brake-of-brake`, `feedback-enabler`)
because it is the error that inverts a module.

## See also

[`dmsa_polarity_review()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_review.md)
for the rows needing adjudication,
[`dmsa_polarity_fetch()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_fetch.md)
to draft signs for a panel of your own,
[`dmsa_polarity_check()`](https://teindor.github.io/dmsa/reference/dmsa_polarity_check.md)
to validate a candidate table.

## Examples

``` r
# \donttest{
pol <- dmsa_polarity()
pol
#> dmsa gene-to-system polarity: Project Alpha 2026c (audited) 
#>   1234 gene-system rows | 30 systems | +1: 583, -1: 291, 0: 360 | 98 anchor(s)
#>   evidence grade: curated 112 (9%) | database 191 (15%) | literature 228 (18%) | heuristic 676 (55%) | none 27 (2%) 
#>   190 row(s) flagged for review - dmsa_polarity_review()
#>   NOTE 55% of signs rest on a functional-class heuristic rather than a database or paper.
#>        Treat system-level direction claims in those systems as provisional.
dmsa_polarity_review(pol)
#> polarity rows needing a decision: 190
#>   low_confidence_signed                  153
#>   pi_adjudicate_disagreement             1
#>   pi_framing_choice                      8
#>   unresolved_needs_evidence              28
#> 
#> the ones that change a sign if you rule the other way:
#>  FOXL2     sys 5   w_g=+0  literature
#>      RESOLVED SET_ZERO: No single sign vs GnRH/gonadotropin tone: pituitary FOXL2 is required for FSHB and GNRHR (+1) while granulosa FOXL2 represses STAR/
#>  ADORA2A   sys 8   w_g=-1  curated
#>      DISAGREE: PI sign kept. A2A is Gs and raises cAMP, so a strict transduction reading of a D1-like/cAMP tone gives +1
#>  ADORA2B   sys 8   w_g=+0  literature
#>      RESOLVED SET_ZERO: Gs coupling is real but sits in low-affinity astroglial A2B, not the striatal cAMP-DARPP-32 output the tone names; the A1/A2A class
#>  REST      sys 9   w_g=+0  literature
#>      RESOLVED KEEP_CURRENT: Signed vs serotonergic tone: TRRUST's direction is wrong - REST represses TPH2 (DN-REST raises it), giving -1, which collides w
#>  PHOX2A    sys 10  w_g=+1  literature
#>      RESOLVED ADOPT_DATABASE: Signed vs noradrenergic output: in adult LC neurons that already exist, raising Phox2a raises DBH/NET, and endogenous Phox2a 
#>  PHOX2B    sys 10  w_g=+1  literature
#>      RESOLVED ADOPT_DATABASE: Same adult LC test as PHOX2A: level moves DBH/NET per cell, so it fails CONVENTION's 'effect via fate not signalling level' t
#>  EP300     sys 20  w_g=+0  literature
#>      RESOLVED SET_ZERO: Signed vs 5mC-writing tone. The EP300->DNMT1 edge is real (19275888) but p300 writes acetyl not 5mC; an edge onto an anchor is not 
#>  FOXP2     sys 25  w_g=+0  literature
#>      RESOLVED SET_ZERO: Signed vs synapse-organising tone: direction is region-dependent (striatal Foxp2-|Mef2c gives +1) and the CNTNAP2 edge is Dual and 
#>  PAX8      sys 26  w_g=+0  rule
#>      RESOLVED KEEP_CURRENT: Signed vs T3-signalling tone. PAX8 does transactivate TG/NIS/TSHR (17614769;14630715) but is the convention's named lineage-det
#> 
#> first 20 of 181 lower-priority rows:
#>  ORAI3     sys 1   w_g=+0 low      unresolved no isoform-specific evidence retrieved on whether higher ORA
#>  PPID      sys 2   w_g=+0 low      unresolved TPR immunophilin competing with FKBP4/FKBP5 on Hsp90; no evi
#>  HSD17B11  sys 4   w_g=+0 low      unresolved reported 17-beta-HSD activity on 3-alpha-androstanediol but 
#>  HSD17B6   sys 4   w_g=+0 low      unresolved carries both retinol-dehydrogenase and oxidative 3-alpha-HSD
#>  ACVR1     sys 5   w_g=+0 low      unresolved ALK2 serves both BMP and AMH signalling with opposite conseq
#>  DHRS9     sys 5   w_g=+0 low      unresolved 3-alpha-HSD/retinol dehydrogenase with no established direct
#>  LGR4      sys 5   w_g=+0 low      unresolved R-spondin receptor with no established directional relation 
#>  VN1R1     sys 6   w_g=+0 low      unresolved no established human ligand or physiological output; deorpha
#>  YY1       sys 9   w_g=+0 low      unresolved no signed HTR1A or SLC6A4 edge in OmniPath and no direction-
#>  TAAR6     sys 11  w_g=+0 low      unresolved no established ligand or transduction direction for human TA
#>  ECE1      sys 12  w_g=+0 low      unresolved M13 peptidase whose canonical substrate is big-endothelin; n
#>  RNPEP     sys 12  w_g=+0 low      unresolved arginyl aminopeptidase; the secretory-vesicle basic-residue 
#>  GPR18     sys 14  w_g=+0 medium   gtopdb     deorphanisation contested (NAGly and THC agonism disputed); 
#>  GRINA     sys 16  w_g=+0 low      unresolved GRINA is TMBIM3/Lifeguard (ER Ca and anti-apoptotic) not an 
#>  BRSK2     sys 18  w_g=+0 low      unresolved SAD-B is redundant with SAD-A for polarity/vesicle clusterin
#>  CEBPA     sys 18  w_g=+0 low      unresolved no signed edge to BDNF/CREB1; OmniPath targets are myeloid l
#>  PHF8      sys 20  w_g=+0 low      unresolved H3K9me1/2 and H4K20me1 demethylase; the principal in vivo su
#>  UHRF2     sys 20  w_g=+0 low      unresolved UHRF1 paralogue binding 5hmC/hemimethylated DNA; whether it 
#>  AHCYL1    sys 22  w_g=+0 low      unresolved AHCY-like paralogue (IRBIT) acting mainly as an IP3-receptor
#>  HLA-F-AS1 sys 24  w_g=+0 low      unresolved antisense lncRNA at HLA-F; no directional evidence vs inflam
# }
```
