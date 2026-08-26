# Adversarial review: exp001b_deconvolved

Reviewer: review-exp001b-2 (independent adversarial reviewer; predecessor died pre-output)
Branch: `reviews/exp001b` | Date: 2026-08-25
Subject: `experiments/exp001b_deconvolved` (worktree ~/ALS-worktrees/exp001b, branch exp/exp001b-retry)

## VERDICT: PENDING — review in progress (this file updated incrementally)

Pre-committed gate claimed by runner: FAIL, 0/3 primary controls qualified
(rule: >=2 of riluzole/baclofen/gabapentin with exposure_matched_percentile<=10
AND FDR<=0.25; on_fail_action: kill the CMap-reversal stream for ALS).

Verification checklist:
- [ ] Gate integrity: config thresholds provenance vs outputs (git --follow)
- [ ] Statistical correctness: independent recompute of exposure-matched
      percentiles for riluzole + tacrolimus from sig_scores.npz + drug_ranking.csv
- [ ] Deconvolution sanity: attribution/composition-shift isolation of neuronal component
- [ ] Final call per lessons L2/L3: CONFIRM STREAM-KILL / OVERTURN / NARROW-EXCEPTION
