# DRAFT — PRO-ACT Data Access Application (H-008 / exp002)
**STATUS: DRAFT — not yet submitted. Owner: biomarker-factory agent (pipeline/biomarker-enrichment).**
Last edited: see git log of this branch. Human co-signatory required before submission.

---

## 1. Requested dataset

- **PRO-ACT Database** (Pooled Resource Open-Access ALS Clinical Trials).
  Aggregated placebo-arm + treatment-arm longitudinal clinical data from
  17(+) completed ALS phase II/III trials; >10,000 subjects, >8,600 with
  sufficient ALSFRS-R follow-up.
- Requested tables: `ALSFRS-R` (all visits), `Demographics`, `Lab` (baseline),
  `Vital Signs`, `Treatment` (arm codes only), `Death/censoring dates`,
  optional `RNA/SOD1` flags if released in current version.
- Format: current PRO-ACT CSV distribution (via the PRO-ACT data server /
  ALS TDI). Local storage will be `data/proact/` (gitignored; never committed,
  per repo rule 1).

## 2. Applicant information

> [TO FILL BEFORE SUBMISSION] Principal investigator name, affiliation,
> ORCID, contact e-mail, research-supervisor attestation if required by the
> current PRO-ACT terms.

## 3. Intended use statement (for the application form)

We request PRO-ACT access for a non-commercial, publish-openly research
project testing whether **trial enrollment enriched on NfL-trajectory features
(baseline NfL + early slope), modeled with machine learning, reduces required
sample size >=25% at equal power relative to conventional ALSFRS-R-slope
enrichment and to unenriched designs** (pre-registered hypothesis H-008,
repository `Project-ALS`, hypotheses/H-008.md).

Specifically we will:
1. Fit prognostic models of individual progression rate on PRO-ACT
   clinical covariates (ENCALS-style feature set, below) as comparator arms.
2. Where PRO-ACT contains plasma/CSF biomarker or validated surrogate fields,
   substitute them for our simulation's NfL blocks; where it does not, keep
   progression-rate surrogates explicitly calibrated to published NfL
   literature (PMIDs 30014505, 31432691, 34690913, 38674431, 35585374) —
   this substitution boundary is declared in every output file.
3. Estimate sample size for 80% power via bootstrap trial resampling under
   treatment-effect scenarios 0.15–0.50 ALSFRS-R points/month.
4. Publish all code and aggregate results openly; no individual-level data
   will be redistributed.

No commercial use. No re-identification attempts. No sharing of raw records
beyond the terms of service.

## 4. ENCALS model parameters needed (comparator arm)

The ENCALS survival model is our pre-registered clinical-only comparator
(PMID 29598923). We require (or will reconstruct from PRO-ACT fields) these
inputs, with their published coefficient structure:

| # | Feature | Source field(s) | Notes |
|---|---------|-----------------|-------|
| 1 | Age at onset | Demographics | years |
| 2 | Site of onset (bulbar/spinal/other) | Demographics | categorical |
| 3 | Time since onset to enrolment | onset date vs first visit | months |
| 4 | Baseline ALSFRS-R total | ALSFRS-R | 0-48 |
| 5 | ALSFRS-R decline rate pre-enrolment | derived | ΔFS/month |
| 6 | Baseline FVC (% predicted) | Vital Signs/Lab | sitting |
| 7 | BMI / weight at enrolment | Vital Signs | |
| 8 | C9orf72 repeat-expansion status | genetics table if present | hexanucleotide carrier flag |
| 9 | Survival time / censoring | follow-up tables | endpoint anchor |

If the current PRO-ACT release lacks any ENCALS inputs (notably C9orf72
status and FVC coverage), the comparator arm will be run on the available
subset and that restriction will be reported in outputs.

## 5. NfL-specific request

PRO-ACT predates routine NfL sampling; if the current release contains no
NfL fields, we will (a) proceed with progression-rate surrogates calibrated
to published NfL kinetics (declared simulated-vs-real boundary per H-008
falsification criterion 1), and (b) separately seek an external NfL-
longitudinal cohort (candidates: ATLAS-type SOD1-carrier datasets PMID
35585374; ENALS/NfL biobank cohorts surfaced in PMIDs 38609644, 38859579)
before any headline claim is finalized.

## 6. Pre-submission checklist

- [ ] PI identity block filled (section 2)
- [ ] Terms-of-service read by a human; this draft matches current terms
- [ ] Adversarial reviewer assigned (`reviews/`) before any verdict is written
- [ ] On approval: download to `data/proact/raw/` (gitignored), record SHA-256
      checksums in `data/proact/inventory.md`, then wire real data behind the
      SAME arm definitions as `pipelines/biomarker_enrichment/run.py`
