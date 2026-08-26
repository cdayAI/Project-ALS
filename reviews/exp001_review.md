# Adversarial review: exp001_repurposing

Reviewer: review-exp001-1 (independent adversarial reviewer)
Branch: `reviews/exp001` | Date: 2026-08-25
Subject: `experiments/exp001_repurposing` (GSE124439 bulk post-mortem CNS signature -> LINCS L1000 Level5 reversal screen, 19,811 small-molecule perturbagens)

## VERDICT: **KILLED** — confidence: HIGH

The current candidate ranking is invalid. The positive-control gate fails once the
ranking's mechanical bias is removed, the primary metric is dominated by a
signature-count artifact, top candidates do not replicate in an independent
ALS dataset, and the annotation step silently cached failed API calls. Per
AGENTS.md rule 5, a pipeline that cannot rediscover known truth produces no
verdict at all. No candidate from this run advances to H-001/H-002 refinement.

This kills the *run*, not necessarily the approach. A pre-registered revision
(next section) is cheap and worth one retry.

## Findings

### 1. Multiple testing / ranking-metric stability — FAIL

- **No p-values or FDR exist for the perturbagen comparison.** BH correction is
  applied only inside the disease DE step. The drug ranking is a raw point
  estimate (-cosine) compared across 19,811 drugs with no null model.
- The primary metric, `best_score = min(score over all signatures of a drug)`,
  is an order statistic: drugs profiled in more signatures are mechanically
  pushed to more extreme minima.
  **Spearman(n_sigs, best_score) = -0.63 (p<1e-300)** over the 19,811 drugs.
  Mean (median - best) gap grows monotonically from 0.00 at n_sigs=1 to 0.26 at
  n_sigs>100 — exactly the pattern of pure exposure bias.
- The `drug_ranking_nsig3_sensitivity.csv` file does NOT control this: its top-10
  is identical to the primary ranking (all members already have n_sigs>=4).
  Filtering n_sigs>=3 removes single-signature drugs but not the min-of-N bias
  among well-profiled drugs.
- Re-ranking by `median_score` (unbiased w.r.t. n_sigs,
  Spearman(n_sigs, median)= -0.01) produces a completely different list, itself
  dominated by 2-4-signature BRD probe codes — i.e., no stable, defensible
  primary list exists in the outputs.
- **Exposure-matched positive controls expose the artifact decisively**
  (percentile among drugs with n_sigs within 2x):

  | compound | raw pctile (best-agg) | exposure-matched pctile |
  |---|---|---|
  | riluzole | 13.9 | **54.2** (n=1,435 matched) |
  | edaravone | 32.9 | 51.0 |
  | dexamethasone | 6.9 | 58.3 |
  | tacrolimus | 4.0 | 26.3 |

  Riluzole's reported "13.9th percentile" is therefore entirely explained by its
  signature count; conditional on exposure it sits at chance.

### 2. Batch effects / confounding — SERIOUS, UNMODELED

- No covariate adjustment for RIN, PMI, age, sex, or cortex subregion (the
  subregion field is parsed into `sample_sheet.csv` but never used in the design
  matrix, which is just `[1, ALS]`). 145-vs-17 imbalance with no variance-partition
  check.
- Cell-composition markers in the shipped disease signature argue the contrast
  is substantially composition/tissue-quality, not neuronal degeneration:
  oligodendrocyte panel mean t = +1.15 (7/7 up), endothelial = -2.27,
  interferon-stimulated genes = -1.76, microglia slightly DOWN (-0.37), neurons ~0.
  Gliosis-driven immune upregulation — the expected ALS post-mortem signal — is absent.
- The top-300 UP enrichment (Myc Targets V1, mTORC1, protein secretion, ER/UPR)
  is equally consistent with a proliferative/biosynthetic and agonal-state signal.
  This matters downstream: anti-proliferative compounds (see finding 5) reverse
  such signatures generically.
- Required check not performed: cell-type deconvolution of the query against a
  spinal-cord/cortex atlas marker matrix (already pre-registered as
  hypotheses/H-001). Until that runs, "reverses the ALS signature" has no
  cell-type interpretation.

### 3. Positive-control adjudication — GATE FAILED

- Reported: riluzole 2758/19,811 ("neutral"), edaravone neutral, dexamethasone
  "reverses", minocycline "mimics".
- Under median aggregation every positive control falls to the 33rd-83rd
  percentile; under exposure matching riluzole/dexamethasone sit at ~54/58th
  percentile. **No compound with a plausible ALS mechanism shows reversal above
  chance in this pipeline.** The apparent top-11% cluster is composed almost
  exclusively of immunosuppressants/anti-proliferatives (tacrolimus,
  dexamethasone, methylprednisolone, ibuprofen, sirolimus).
- Could "riluzole need not reverse dead tissue" rescue this? Partially — but the
  same insensitivity applies to the whole panel, and dexamethasone/minocycline
  were included precisely because their transcriptional footprints should be
  recoverable. When the two "hits" (dex, tacrolimus/sirolimus) are exactly the
  compounds whose L1000 footprints suppress Myc/mTORC1/protein-synthesis programs
  — the dominant enriched axis of our query — the parsimonious explanation is
  assay-class artifact, not neuroprotection detection.
- Against the standard set by `research/02_ai_methods.md` sec 4-5: clemastine and
  ropinirole succeeded when signature logic was fused to a functional assay;
  bexarotene shows how fast signature-plausibility stories collapse without one.
  This run has neither functional support nor a positive control that survives
  fair scoring. It cannot generate triage-worthy candidates as-is.

### 4. Replication attempt vs GSE255602 (iPSC motor neurons, 27 ALS / 44 ctrl)

I recomputed the full pipeline independently on `data/gse255602/` (salmon gene
counts -> CPM filter -> log2 -> moderated t -> map to 9,545 L1000 genes ->
rescored all 473,647 Level5 signatures -> per-drug best/median). Results:

- **Global concordance is moderate**: Spearman rho = 0.70 (best-score) / 0.64
  (median-score) across all 19,811 overlapping drugs. The geometry is stable.
- **Top hits do not replicate**: overlap of the top-100 lists is 3 drugs
  (~1 expected by chance); top-500 overlap 54 (~13 expected). The replication
  screen's own top ranks (tretinoin, apicidin, batimastat, teniposide...) share
  essentially nothing with exp001's (BRD probes, tw-37, trifluoperazine...).
- Riluzole again lands at chance-ish neutral (15.6th percentile raw best-agg) in
  the replication — consistent with finding 3: both signatures fail the same gate.
- Immunosuppressants (sirolimus 0.2%, tacrolimus 1.0%, dexamethasone 0.5%) top
  the iPSC-MN ranking too — confirming these are generic footprint-reversers,
  not ALS-specific signals.

Honest summary: rankings correlate genome-wide because both queries share weak
global structure and the same n_sigs artifact; the actual candidates are
dataset-specific noise.

### 5. Top-candidate sanity — FAIL

Of exp001's top-50: ~24 are unannotated BRD/SA/MW probe codes; several others
are classic assay-artifact classes (tw-37 Bcl-2-family cytotoxic;
geldanamycin HSP90 cytotoxic; nsc-95397 cytotoxic phosphatase inhibitor;
withaferin-a in the shared top-500; genistein broad promiscuous kinase inhibitor).
The most "real" names (trifluoperazine #7, nicardipine, naltrexone, calcitriol,
estradiol) are old CNS/steroid agents with huge pleiotropic footprints — no
mechanistic ALS link beyond generic membrane/calmodulin effects.
Additionally, `top_candidates_annotated.csv` is misleadingly empty: ChEMBL/PubChem
lookups timed out during the run and the errors were **permanently cached as
empty results** (`annotation_cache.json`: 26/30 lookups contain error strings,
not data). Checklist item 5 ("ChEMBL annotations present") fails for reasons that
are part bug, part candidate quality.

## What was good

- Reproducibility contract is real: config.yaml + run.py regenerated everything I
  checked; TE-feature removal before normalization is thoughtful and correct.
- Honest self-flagged caveats in README.md (nsig bias, no covariates) match what
  I found — the experiment did not oversell itself.
- Streaming scorer and DE implementation are correct (verified by independent
  reimplementation producing consistent scores).

## Required next steps (pre-register before any re-run)

1. **Null-calibrated ranking**: empirical p-values for each drug via
   signature-label permutation or a null score distribution matched on n_sigs
   (e.g., rank-based statistic such as rank-product/median of within-cell-line
   percentile ranks, tested against exposure-stratified nulls). Report FDR.
   Drop `min`-over-signatures as a headline metric.
2. **Query QC gate**: deconvolve the GSE124439 signature against an atlas marker
   matrix (per H-001/GSE221692); require the neuronal or disease-relevant
   component to carry a defined fraction of the query norm before any reversal
   claim. Model subregion at minimum; drop samples without RIN metadata if
   obtainable from the supplement of the source paper.
3. **Positive-control gate with teeth**: define pass/fail BEFORE re-run
   (e.g., >=2 of {riluzole-or-analog, EAAT2-inducing class, dexamethasone} in the
   top decile exposure-matched). If the gate fails again, kill the CMap-reversal
   stream for ALS rather than re-tuning.
4. **Fix the annotation cache**: do not cache transport/API errors as negative
   results; key lookups by normalized name and retry failures.
5. Route any surviving shortlist through a functional-assay feasibility filter
   per research/02 sec 5 before it enters H-001/H-002 refinement.

## Disposition

- exp001 outputs: **KILLED** (no candidate survives; ranking metric invalid).
- H-001 (deconvolved queries) remains open and is now the mandatory prerequisite
  for any retry; H-002 inherits the same gate.
- The reviewer-built GSE255602 re-scoring path can be contributed to
  `pipelines/signature_reversal/` if maintainers want it as the replication arm.
