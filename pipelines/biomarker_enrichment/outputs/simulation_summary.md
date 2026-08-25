# H-008 NfL-enrichment simulation summary

> **DATA STATUS: 100% SIMULATED.** No PRO-ACT or other real patient
> records were used anywhere in this run. Every number below is a
> property of the generative model in `config.yaml`, calibrated to
> published summaries (PMIDs 25362243, 29598923, 35585374, 30014505,
> 31432691, 34690913, 38674431, 31280619). Real-data validation is
> **PENDING**: see `data/proact/APPLICATION_DRAFT.md`.

## Effect-size assumption
Treatment slows ALSFRS-R decline by an absolute [0.15, 0.3, 0.5] pts/month scenarios;
the headline claim uses 0.30 pts/month (H-008 realistic range 0.3-0.5).

- Pool calibration: {'corr_log_nfl_baseline_vs_log_lambda': 0.593, 'median_progression_pts_per_month': 0.797}
- Eligible subjects per arm (fraction 0.5): {'none': 10000, 'noise': 10000, 'oracle': 10000, 'clinical_static': 10000, 'nfl_ml': 10000}
- Nested-model test, NfL slope beyond baseline: p=1.45e-96

| effect (pts/mo) | arm | n/arm @80% power | reduction vs none | screened/enrolled |
|---|---|---|---|---|
| 0.15 | none | 186.5 | 0.0 | 2.0 |
| 0.15 | noise | 186.9 | -0.002 | 2.0 |
| 0.15 | clinical_static | 153.1 | 0.179 | 2.0 |
| 0.15 | nfl_ml | 126.7 | 0.321 | 2.0 |
| 0.15 | oracle | 56.1 | 0.699 | 2.0 |
| 0.3 | none | 40.2 | 0.0 | 2.0 |
| 0.3 | noise | 42.2 | -0.05 | 2.0 |
| 0.3 | clinical_static | 34.4 | 0.145 | 2.0 |
| 0.3 | nfl_ml | 28.7 | 0.285 | 2.0 |
| 0.3 | oracle | 14.0 | 0.653 | 2.0 |
| 0.5 | none | 12.0 | 0.0 | 2.0 |
| 0.5 | noise | 12.3 | -0.023 | 2.0 |
| 0.5 | clinical_static | 11.1 | 0.08 | 2.0 |
| 0.5 | nfl_ml | 9.4 | 0.221 | 2.0 |
| 0.5 | oracle | 5.6 | 0.535 | 2.0 |

## Pending real data
- Re-run with PRO-ACT placebo-arm trajectories once access is granted
  (application draft: `data/proact/APPLICATION_DRAFT.md`).
- Replace simulated NfL blocks with observed PRO-ACT/NfL values where
  available, else keep progression-rate surrogates calibrated to the
  NfL literature (per H-008 falsification criterion 1).
- ENCALS comparator parameters requested in APPLICATION_DRAFT.md.
