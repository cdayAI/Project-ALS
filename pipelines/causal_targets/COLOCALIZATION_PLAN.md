# Colocalization plan: ALS GWAS loci × Project MinE rare-variant burden

Status of this document: a **plan**, not a result -- per the task board,
this is "Stream A scaffolding," and Project MinE access is still an
application-pending controlled-access resource (`research/03_data_resources.md`
§1.1, EGA Data Access Committee). No MinE genotype or summary-statistics
data has been requested or obtained as of this writing. This document exists
so that the moment access is granted, there is a concrete, pre-registered
analysis plan to run -- not a design exercise done after the data arrives
(which risks post-hoc rationalization of whatever comes out).

## Why this is a separate axis from `fetch_gwas_targets.py`, not a duplicate

`fetch_gwas_targets.py` surfaces genes with **common-variant** GWAS
association signal (SOD1, C9orf72, UNC13A, TBK1, ...). Running it live and
checking its own positive control (see `outputs/positive_control.json`)
confirmed TARDBP and FUS -- both core ALS genes per
`research/01_biology_and_therapeutics.md` -- have **zero** GWAS Catalog
association signal for ALS, because they were discovered via family-linkage
sequencing of rare, highly-penetrant variants (TARDBP ~3-5% of familial ALS,
FUS ~0.3-0.9% of all ALS), not population case-control GWAS, which needs
common allele frequencies for statistical power. This is expected, not a
bug -- but it means a GWAS-only pipeline is structurally blind to exactly
the kind of gene MinE's whole-genome rare-variant burden testing exists to
find. The two resources are complementary, not redundant:

| | Common-variant GWAS (`fetch_gwas_targets.py`) | Project MinE rare-variant burden |
|---|---|---|
| Variant frequency | Common (population-scale MAF) | Rare, often private/family-specific |
| Discovers | SOD1, C9orf72, UNC13A, TBK1, MOBP, ... | TARDBP, FUS, and novel rare-variant genes |
| Access | Open, no auth (GWAS Catalog REST API) | Controlled (EGA DAC) |
| This repo's status | **Done** -- `outputs/gwas_druggable_targets.csv` | **Pending** -- this plan |

## Analysis plan (to run once MinE summary stats / genotypes are available)

1. **Input**: Project MinE WGS-derived rare-variant burden test results
   (gene-level, e.g. SKAT-O or burden p-values) or raw genotype VCFs if
   summary stats aren't released directly -- confirm which access tier MinE
   grants before assuming either.
2. **Rare-variant gene list**: genes with burden-test significance below a
   pre-registered threshold (e.g. genome-wide-corrected p < 2.5e-6 for
   ~20,000 genes, Bonferroni; document the actual threshold used at
   analysis time, per AGENTS.md rule 6).
3. **Cross-reference against `outputs/gwas_druggable_targets.csv`**:
   - Genes appearing in BOTH lists (common + rare variant signal at the
     same locus) get the strongest causal-evidence tier -- this is the
     GWAS+rare-variant convergence pattern that validated SOD1 historically.
   - Genes appearing ONLY in the MinE rare-variant list (expected: TARDBP,
     FUS, and any novel genes) are the primary yield of this step --
     exactly what the GWAS-only pipeline structurally cannot find.
4. **Tractability overlay**: run the same Open Targets tractability lookup
   (`resolve_ensembl_id` / `tractability_summary` in `fetch_gwas_targets.py`
   -- reuse directly, don't reimplement) on the MinE-only gene list, so the
   rare-variant output is druggability-annotated the same way as the GWAS
   output.
5. **Positive control for this step specifically**: TARDBP and FUS
   themselves. If a MinE rare-variant burden run does NOT recover TARDBP
   and FUS as significant, the burden test's parameters (variant
   inclusion criteria, functional annotation filter, cohort ancestry
   composition) need review before trusting any novel gene it nominates --
   per AGENTS.md rule 5, this is non-negotiable, not optional.
6. **Adversarial review requirements** (AGENTS.md rule 4) specific to this
   analysis: multiple-testing correction choice, population-stratification
   /ancestry confounding (MinE is European-ancestry-dominant per
   `research/03_data_resources.md` §1.7's diversity caveat), and whether
   any novel hit replicates in NYGC ALS Consortium WGS
   (`research/03_data_resources.md` §1.5, BioProject PRJNA573105 -- already
   openly accessible, no DAC needed, and a natural independent replication
   cohort for whatever MinE nominates).

## What can happen before MinE access arrives

- Apply to the Project MinE Data Access Committee now (`research/03
  §1.1`'s own guidance: "DAC approval typically takes weeks-months; budget
  that time up front"). This repo has not yet recorded that an application
  was submitted -- doing so is the actual highest-value next action for
  this plan, not further pipeline engineering.
- The NYGC ALS Consortium WGS data (open, BioProject PRJNA573105) could
  serve as a stand-in rare-variant cohort to prototype and validate steps
  2-5 of this plan end-to-end before MinE access lands, so the analysis
  code is tested and the positive control (step 5) is already known to
  pass on a real dataset before it matters on the real target cohort.
