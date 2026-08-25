# Experiment: exp002_nfl_enrichment — H-008 NfL-trajectory trial enrichment

**STATUS: SKELETON — DO NOT RUN.** Pipeline exists and runs end-to-end on
SIMULATED data (`pipelines/biomarker_enrichment/run.py`, exit 0 with passing
positive/negative controls). This experiment is NOT started until PRO-ACT
access is granted (application draft: `data/proact/APPLICATION_DRAFT.md`).

config.yaml   - full parameters, dataset IDs, versions (reproducibility contract)
run.py        - entry point; must run end-to-end via a single command
outputs/      - raw outputs only; no manual editing
VERDICT.md    - filled AFTER adversarial review: KILLED / SURVIVES / NEEDS-DATA

## Data boundary (read first)
- Current pipeline state: 100% simulated cohorts. Every output file carries a
  `DATA STATUS: SIMULATED_ONLY` banner; nothing here is evidence about real
  ALS trials.
- Real-data state: PENDING. PRO-ACT application not yet submitted. When data
  arrives it enters `data/proact/raw/` (gitignored) and the same arm
  definitions are re-run on real placebo-arm trajectories.

## Design (pre-registered comparator arms)
1. `none`            — unenriched enrollment (reference)
2. `noise`           — ranking on pure noise (negative control)
3. `clinical_static` — ridge ML on observed 0-3 mo ALSFRS-R slope + baseline
                       score (conventional clinical enrichment comparator)
4. `nfl_ml`          — ridge ML on log-baseline NfL + early NfL slope +
                       clinical features (H-008 candidate)
5. `oracle`          — ranks on TRUE progression rate (positive-control bound)
Endpoint: change in ALSFRS-R at 12 months, two-sided alpha=0.05,
target power 0.80, bootstrap resampling for n-at-power curves.

## Positive control
The oracle arm (selection on true latent progression rate) must reduce
required n by >=25% vs unenriched at every effect scenario, and must dominate
all learned arms; the noise arm must show ~0% reduction (|delta| <= 5%).
A pipeline failing either produces NO verdict (AGENTS.md rule 5).
Current simulated-data run: controls PASS (see
`pipelines/biomarker_enrichment/outputs/simulation_summary.md`).

## Falsification criteria (from hypotheses/H-008.md, verbatim intent)
Exact pre-registered criteria (H-008, criterion 1 verbatim):
- F1: If simulation on PRO-ACT placebo arms (with NfL where available, else
  progression-rate surrogates calibrated to NfL literature) shows ML-NfL
  ALSFRS-R-slope enrichment => headline quantitative claim killed.
  NOTE on comparator: H-008 names "standard ALSFRS-R-slope enrichment" as the
  reference for this criterion; the >=25% figure is therefore evaluated as
  nfl_ml vs clinical_static AND reported additionally vs none. On current
  SIMULATED data: nfl_ml-vs-clinical_static = 28.7/34.4 -> ~16.6% further
  reduction at eff=0.30 pts/month — i.e., the >=25%-vs-static claim does NOT
  hold in simulation; only the vs-none version clears 25%. This gap is the
  primary thing real PRO-ACT data must resolve.
- F2: If NfL-based enrichment biases the enrolled population toward slow
  progressors so severely that absolute treatment effect (patients
  saved/month) shrinks despite smaller n => killed. Operationalized: report
  mean absolute ALSFRS-R benefit per enrolled patient per arm; kill if
  nfl_ml absolute per-patient benefit < clinical_static while n reduction
  >25%.
- F3: If early NfL slope adds no predictive value beyond baseline NfL alone
  (nested likelihood-ratio test P > 0.05) => dynamics framing collapses;
  killed as stated. Simulated-data diagnostic currently p<1e-90 BY
  CONSTRUCTION (slope info was built into the generator); the real-data test
  is decisive and this simulation cannot falsify F3.

## VERDICT
PENDING — requires (a) real PRO-ACT data, (b) adversarial review in
`reviews/exp002_nfl_enrichment.md`. No self-certification.
