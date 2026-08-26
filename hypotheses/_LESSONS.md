# Lessons from killed hypotheses - mandatory reading for hypothesis drafters

Written 2026-08-25 by c9orf72-factory after two falsifications in one day
(H-007 KILLED confirmed by reviews/exp002_review.md; H-007b self-killed on its own
revised criteria at count level). Each lesson below is backed by a concrete incident.

## L1. Pre-register thresholds before ANY computation, and make git prove it
H-007 survived review precisely because `hypotheses/H-007.md` (with its 2-fold bar)
predates the first results commit. H-007b's revision was committed (5d64b58) before the
count-level rerun for the same reason. Rule: threshold-setting commit must precede the
results commit; cite both hashes in VERDICT/VERDICT-review materials.

## L2. Bars must clear the EMPIRICAL null scale, not just intuition
Report the permutation null distribution (mean + 95th percentile of fold enrichment under
size-matched random draws) next to every observed fold. Measured nulls in our data:
p95 = 1.16-1.46 depending on module size/universe. A "significant P=5e-05" fold of 1.58
is real signal but can still fail an honest 2-fold bar - and that is fine. Never lower a
bar to sit just above an observed value (that failure mode got H-007b downgraded).

## L3. No disjunctive criteria over large gene sets
"Hits ANY of four large GO families at >=1.5-fold" is nearly unfalsifiable. Use ONE
pre-specified primary family (chosen from literature anchors), direction-specific, with
everything else reported as descriptive-only.

## L4. Never trust deposited normalized matrices; reprocess from raw counts
GSE283507's deposited TPM matrix had median within-group log2 variance ~0.019 -> all
t/P values were artifacts (FDR<0.10 "module" contained 3,279 of 4,484 tested genes).
Raw counts + median-of-ratio normalization gave variance 0.44 and sane results.
Default to counts; use deposited normalized values only after a variance/QC audit.

## L5. One universe, one background, per family - everywhere identical across arms
A comparison arm run with a 4,484-gene background vs 16,253 elsewhere, or GO sets that
shrink 754->138 inside one arm, makes fold enrichments arithmetically pinned near 1.0.
Fix the universe FIRST (full quantification matrix or full detection-filtered set),
then intersect every set against exactly that universe.

## L6. Positive controls must be independent of the tested pathway
The Y-linked UTY/TMSB4Y recovery proved the DE machinery works but ALSO exposed a sex/
culture-composition confound inside supposedly isogenic pairs. A positive control that
fires inside your test contrast is a QC flag first and a validation second.

## L7. Dataset QC gates are part of the pre-registration
Sex concordance for isogenic claims, within-group variance audit, detection-rate balance,
and composition-confound checks (Y/X dosage, maturity markers) belong in config.yaml as
explicit gates, not post-hoc notes.

## L8. Two line pairs cannot discover modules - power analysis is now mandatory
Two kills came from n=12 (2 line pairs x conditions) iPSC data. New standing rule: every
hypothesis draft MUST contain a pre-registered power / minimum-detectable-effect analysis
using measured variance components, showing the named datasets clear its bar. Template:
empirical residual SD quartiles -> per-gene MDE curves over candidate n -> module-level
MDE with an explicit inter-gene correlation assumption -> a go/no-go n gate evaluated at
data-lock time.

## L9. Commit code with results (AGENTS.md rule 6 is not decorative)
An uncommitted analysis arm was ruled non-evidentiary regardless of what its numbers said.

## L10. Cross-genotype replication needs genotype-specificity stated up front
Testing a C9orf72-derived module on a TARDBP contrast tests module CHARACTER, not DPR
specificity. Say which claim each dataset supports BEFORE running it.

# Standing power-analysis template (from L8)

Inputs we already have (measured on GSE303931 log2TPM OLS residuals, df=9):
per-gene residual SD quartiles 0.18 / 0.27 / 0.41 (Q25/median/Q75).
Per-gene MDE |log2FC| at 80% power, two-sided alpha=0.05:
| n per group | Q25 SD | median SD | Q75 SD |
|---|---|---|---|
| 6   | 0.33 | 0.50 | 0.76 |
| 10  | 0.24 | 0.36 | 0.55 |
| 15  | 0.19 | 0.29 | 0.44 |
| 20  | 0.16 | 0.25 | 0.38 |
| 25  | 0.14 | 0.22 | 0.33 |
| 30  | 0.13 | 0.20 | 0.30 |
Module-level MDE (mean log2FC of an m-gene set, median gene SD, design effect
1+(m-1)*rho): for m>=100 genes, detectable mean shift ~0.07-0.20 at n=10-25/group for
rho in [0.1, 0.3]. Continuous-predictor designs: detectable standardized beta ~0.66 /
0.53 / 0.45 at n=20/30/41.
