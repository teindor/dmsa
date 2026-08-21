# Selection sets: system \> module \> gene \> probe

DMSA corrects multiplicity inside **declared** families — the genes of a
named system, the probes of a named gene. The selection set is where
those families are declared, and it is the reason a DMSA p-value does
not move when you load a different set of probes: the family comes from
this file, not from your session.

The package ships one, so most users never build their own.

------------------------------------------------------------------------

## Using the bundled set

``` r

library(dmsa)

dmsa_systems()          # the 30 systems and the short name for each
```

      1  oxytocin         Oxytocin, vasopressin & CICR                   9 mod    90 gene   1690 cpg
      2  hpa              HPA axis & glucocorticoid signalling           8 mod    49 gene    538 cpg
      3  steroidogenesis  Steroidogenesis (adrenal & gonadal)            6 mod    35 gene    375 cpg
     ...
     24  immune           Immune, inflammation & HLA                    13 mod   127 gene   1084 cpg
     30  pharmacogenes    Pharmacogenes                                  3 mod    12 gene     74 cpg

Then name the system or systems you are asking about:

``` r

frame <- dmsa_frame(
  data     = dat,
  systems  = c("hpa", "oxytocin"),     # <- short names
  outcome  = "pcl5",
  module   = TRUE
)
```

**Everything below the system defaults to full.** Choosing `"hpa"` takes
all 8 modules, all 49 genes and all 538 CpGs the set assigns to the HPA
axis. You do not list them, and you cannot forget one.

`systems =` accepts, in this order: the short name (`"hpa"`), the system
id (`2` or `"2"`), the full name
(`"HPA axis & glucocorticoid signalling"`), or a unique prefix of a
short name (`"oxy"`). Matching ignores case. Two things always error
rather than guess:

``` r

dmsa_select(systems = "ox")
#> Error: `systems = "ox"` is ambiguous - it matches 'oxytocin', 'oxidative';
#>   use one of those short names

dmsa_select(systems = "cortisol")
#> Error: `systems` did not match: 'cortisol'
#>   available short names: oxytocin, hpa, steroidogenesis, ...
```

### Looking before you run

``` r

sel <- dmsa_select(systems = c("hpa", "oxytocin"))
sel
```

    dmsa selection from: Project Alpha 2026c (module-audited)
      2 system(s) selected | modules full | genes full | probes full
       hpa              HPA axis & glucocorticoid signalling      8 mod    49 gene    538 cpg
       oxytocin         Oxytocin, vasopressin & CICR              9 mod    90 gene   1690 cpg
      total: 17 modules, 139 genes, 2228 CpGs
      module evidence: 17 High, 0 Moderate

### Narrowing a level

Allowed, but understand what it does: narrowing changes the declared
family and therefore the multiplicity toll. Declare a narrower question
in advance if you intend to ask one.

``` r

dmsa_select(systems = "hpa", modules = c("2.6", "2.8"))
dmsa_select(systems = "immune", genes = c("IL6", "TNF", "CRP"))

# pass a narrowed selection to the frame
frame <- dmsa_frame(dat, sets = dmsa_select(systems = "hpa", modules = "2.6"),
                    outcome = "pcl5")
```

### Keying the set to your data columns

The bundled set carries the data-matrix column name for each probe under
five analyses (`col_parent_T1`, `col_long`, `col_child_T4`,
`col_child_maternal_T1`, `col_child_paternal_T1`). Let DMSA work out
which one your matrix uses:

``` r

dmsa_select(systems = "hpa", columns = "auto", data = dat)$column_key
#> [1] "col_parent_T1"
```

------------------------------------------------------------------------

## Module evidence is part of the result

Every module label in the bundled set was checked against the
literature. The check recorded whether the label survived, whether the
locked gene membership is homogeneous, how strong the evidence for the
*definition* is, and what it rests on. DMSA reports that next to your
result, because a module-level finding is only as good as the module’s
definition.

``` r

dmsa_evidence(sel, which = "flagged")
```

    module evidence (3 modules)
     24.11   Moderate T-cell receptor loci and CCL2 chemotaxis    [Loz07]; [New12]
              status: renamed_for_mechanistic_precision
     24.12   Moderate Leukocyte activation and inflammatory sign  [Loz07]; [Pla18]
              status: renamed_from_association_label_to_process_label
     24.13   Moderate CCL3 chemokine signaling                    [Loz07]; [New12]
              status: renamed_for_mechanistic_precision
     evidence searches:
       https://app.undermind.ai/projects/...

- `which = "all"` — every module in the selection
- `which = "moderate"` — Moderate evidence only
- `which = "flagged"` — Moderate, plus modules whose locked membership
  is heterogeneous or whose grouping is measurement-defined rather than
  mechanistic

`evidence_strength` is **High** when the locked genes form one coherent
subprocess with direct literature support, **Moderate** when the
grouping is defensible but the members are mechanistically mixed or the
supporting evidence is thinner. `audit_status` records what happened to
the label: `supported_as_named`, one of the `renamed_*` values, or a
membership flag such as `heterogeneous_locked_membership` or
`measurement_defined_module`.

The counts appear in [`print()`](https://rdrr.io/r/base/print.html) for
the set, the selection and the frame, and `summary.md` from
[`dmsa_report()`](https://teindor.github.io/dmsa/reference/dmsa_report.md)
carries a **Module evidence** table for the modules the run actually
reported on.

------------------------------------------------------------------------

## Bringing your own set

``` r

dmsa_sets_template("my_sets.csv")   # writes the CSV skeleton + these instructions
dmsa_sets_check("my_sets.csv")      # validate before you trust it
cas <- dmsa_sets("my_sets.csv")     # load; prints the short names it found
dmsa_systems(cas)

frame <- dmsa_frame(dat, sets = cas, systems = c("my_system"), outcome = "y")
```

One row per CpG. The file answers four questions at once: which system,
which module inside it, which gene, which CpG.

### Required columns

| column | what it holds |
|----|----|
| `system_id` | stable id, e.g. `2`. Held as **text** — `24.10` must not become `24.1` |
| `system` | the printed system name |
| `module_id` | stable id nested under the system, e.g. `2.6` |
| `module` | the printed module name |
| `gene` | gene symbol |
| `cpg` | CpG id, e.g. `cg04712664` |

### Optional, and what each one buys you

| column | what it buys |
|----|----|
| `system_short` | the handle you type in `systems = c("hpa")`. Derived from `system` if absent, but the derivation of a long name is clumsy — set it yourself |
| `probe_id` | when one CpG is measured by more than one probe. Defaults to `cpg` |
| `col_*` | one column per analysis holding the **data-matrix** column name for that probe. Lets `columns = "auto"` key the set to your matrix instead of you renaming columns |
| `evidence_strength` | `High` or `Moderate` — how well supported the module’s *definition* is |
| `audit_status` | `supported_as_named`, `renamed_*`, `heterogeneous_locked_membership`, `measurement_defined_module` |
| `citation_keys` | the references behind the label, e.g. `[Fri17]; [Fad23]` |
| `evidence_note` | one sentence a reader can check |
| `deep_search_url` | where the literature check lives |

The evidence columns may instead live in a separate one-row-per-module
table:

``` r

dmsa_sets("my_sets.csv", audit = "my_module_evidence.csv")
```

Either way DMSA surfaces them wherever a module-level result is printed.
A set with no evidence columns works fine — the banner then says the
modules are not annotated, which is honest rather than silent.

### What the validator enforces

1.  Every `module_id` belongs to exactly one `system_id`.
2.  Every `gene` belongs to exactly one `module_id`. A gene in two
    modules makes its family ambiguous, and there is no principled way
    to charge the toll.
3.  `system_short` is unique across systems.
4.  Module metadata is constant within a `module_id`.
5.  No duplicate `module_id` + `gene` + `cpg` rows.

A CpG **may** appear under two genes — overlapping annotation is real,
and the bundled set has 228 such CpGs. The validator reports them and
lets you decide.

### This is not the direction map

Three files, three questions, three kinds of evidence:

| file | question | argument |
|----|----|----|
| selection set | what belongs together? | `sets =` |
| direction map | which way does each probe point? | `map =` |
| polarity table | how does each gene relate to its system’s tone? | in the reference bundle |

Keep them separate. They are curated from different evidence and revised
on different schedules, and conflating them is how a set-level p-value
starts depending on things it should not.
