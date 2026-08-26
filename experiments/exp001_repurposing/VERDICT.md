# VERDICT: KILLED (adversarial review 2026-08-25)

Reviewed by: review-exp001-1 (`reviews/exp001_review.md`, branch `reviews/exp001`).

**The candidate ranking from this run is invalid and no output advances.**

Key reasons (details and evidence in the review file):
1. Positive-control gate FAILED: exposure-matched, riluzole ranks at the 54th
   percentile (raw "13.9th percentile" was a min-of-N signature-count artifact,
   Spearman(n_sigs, best_score) = -0.63). No mechanism-plausible compound beats
   chance under fair aggregation.
2. Top-100 candidates do not replicate against an independently recomputed
   GSE255602 (iPSC motor neuron) screen (overlap 3 vs ~1 expected by chance),
   despite moderate global rank concordance (Spearman rho = 0.70).
3. The disease query carries cell-composition/tissue-quality signal
   (oligodendrocyte up, endothelial/interferon down; no covariates modeled);
   deconvolution check per hypotheses/H-001 is mandatory before any retry.
4. `top_candidates_annotated.csv` is empty for 26/30 drugs because failed API
   calls were permanently cached as negative results.

Disposition:
- This run is killed; do not cite its drug list.
- A pre-registered revision is permitted (null-calibrated exposure-matched
  ranking + deconvolved query + pre-committed positive-control gate), tracked as
  prerequisite work in H-001. If the gate fails again after one revision, the
  CMap-reversal stream for ALS is retired.
- Historical note: the pipeline infrastructure itself is sound and reproducible;
  the failure is in statistics and query construction, not engineering.
