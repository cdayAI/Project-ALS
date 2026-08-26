# Adversarial review: exp001b_deconvolved

Reviewer: review-exp001b-2 (independent adversarial reviewer; predecessor died before writing anything)
Branch: `reviews/exp001b` | Date: 2026-08-25
Subject: `experiments/exp001b_deconvolved` (branch `exp/exp001b-retry`, worktree ~/ALS-worktrees/exp001b)
Method: artifacts judged directly from outputs/ + code + git history; README/VERDICT polish of the
terminated runner was incomplete and was NOT weighted.

## VERDICT: **CONFIRM STREAM-KILL** — confidence: MEDIUM-HIGH

But NOT for the reason claimed. The claimed evidence (gate_result.json FAIL,
0/3 controls) is **VOID**: the scoring stage silently crashed (query mapped to
0/15,000 L1000 genes -> every score exactly 0 -> every empirical p exactly 1.0).
The gate "FAIL" carries zero evidential weight and is hereby overturned as
evidence. The stream still dies on an independent, pre-registered ground the
runner recorded honestly: the H-001 variance floor FAILED (overall R2 = 0.0885
< 0.10, reproduced independently), i.e., the bulk ALS signature has essentially
no recoverable cell-type-deconvolved structure, so a "deconvolved neuronal
reversal query" has no premise. Per AGENTS.md rule 5, this produces no verdict
in favor of any candidate. Any future revival of the reversal idea must be a NEW
pre-registered hypothesis (different query source, e.g., cell-type-pure data);
fixing run.py alone cannot rescue it, because the fatal negative is upstream of
the drug scoring entirely.

## Finding 0 (process): gate pre-registration integrity — PASS

- `git log --follow` on config.yaml shows exactly two commits:
  70553c5 "pre-register config with positive-control gate BEFORE scoring"
  (thresholds: >=2 of riluzole/baclofen/gabapentin, exposure_matched_percentile<=10,
  FDR<=0.25, on_fail kill-without-re-tuning) and 5367391 (pipeline code + stage A/B
  snapshots). The only config edit after pre-registration is a cosmetic key-value
  rename `neuronal_component` -> `neuron_component` to match code, committed
  21:32, BEFORE the first output file mtime (21:27-21:35 window; scoring pass
  logged 21:34:57-21:35:20). Thresholds never changed.
- Caveat: outputs/ is gitignored (`experiments/*/outputs/`), so gate_result.json
  is untracked and its provenance rests on mtimes + internal consistency with the
  run log, not git history. Internal consistency checks out (log line, npz
  contents, and ranking table agree exactly).

## Finding 1: the gate FAIL is a mechanical artifact, not a result — FATAL TO CLAIMED OUTCOME

The run log's first line says everything:

    [21:34:57] [primary] query mapped to 0/15000 L1000 columns

Root cause, confirmed by direct inspection:
- LINCS Level5 GCTX row IDs are Entrez IDs ('5720', '466', ...); the deconvolved
  query uses HGNC symbols ('A1BG', ...). Overlap: **0/15,000**.
- `run.py::score_with_nulls` builds the query vector over GCTX rows and never
  loads `lincs.gene_info`, despite config.yaml listing it (dead config key).
- There is NO guard: hit==0 is logged, then execution continues. q = zero vector
  => every cosine = 0/-norm = -0.0.
Consequences verified independently from artifacts:
- sig_scores.npz: 205,034 signature scores, **one unique value (-0.0)**.
- drug_ranking.csv: all 19,811 drugs have median_score = 0.0, empirical_p = 1.0
  (= (1+1000)/(1000+1), the ceiling of the permutation p under total ties),
  fdr = 1.0, exposure_matched_percentile = 100.0.
- Independent recompute from raw npz (groupby perturbagen, median, exposure-match
  |log2 n_sigs ratio| <= 1, same rule as config): identical degenerate table.
  riluzole (n_sigs=45, 1,435 peers), tacrolimus (n_sigs=110, 585 peers),
  baclofen (n=3), gabapentin (n=4): ALL percentile 100.0 by pure tie-breaking;
  reported "ranks" (riluzole 18,684 vs baclofen 331) are arbitrary stable-sort
  order among identical values. Null calibration did not fail statistically —
  it never ran on any signal.
This is precisely the failure class review-exp001-1 flagged in the annotation
stage of exp001 (silent caching of failures); here it hit the primary metric.
A gate evaluated on a crashed scorer cannot trigger "kill" OR "pass"; it can
only say INVALID RUN.

## Finding 2: the stream still dies — on the pre-registered H-001 variance floor — CONFIRMED

Independent of the scoring crash, stage B produced a genuine, correctly-computed
negative:
- deconvolution_summary.json: overall_R2 = 0.08848, h001_variance_floor_met =
  false against the PRE-REGISTERED floor min_variance_explained: 0.10
  (comment in config: "H-001 falsification: <10% => premise fails, record
  honestly"). Recorded honestly it was.
- Reproduced independently: OLS of disease moderated-t (15,000 shared genes) on
  the six cell-type components from signature_components.csv.gz gives
  R2 = 0.0885. Every single cell type's partial R2 (0.069 neuron down to 0.083)
  is below the floor too; univariate neuron-vs-t R2 is 0.021.
- composition_shifts.csv: NO ALS-vs-control neuronal composition difference
  (neuron delta = +0.036, MW p = 0.86; nothing survives except endothelial at
  uncorrected p = 0.027). Consistent with the original exp001 review's suspicion
  that this bulk contrast carries little cell-type-resolved biology.
Interpretation: the deconvolution machinery is sane enough to demonstrate its own
premise is false on this dataset (embryonic spinal-cord reference GW7/GW9 vs
adult post-mortem bulk may contribute, but that reference was part of the
pre-registration). A reversal screen whose query explains 8.8% of signature
variance cannot support a positive-control argument even if run.py were fixed.

## Finding 3: deconvolution sanity — QUALIFIED PASS on mechanics

- atlas assignment: neuron is the largest assigned cluster (4,858 cells, marker
  mean-z 1.61; best-cluster z 1.65) — neurons were isolated plausibly.
- Weaknesses: 38.9% of cells UNASSIGNED; endothelial cluster missing PECAM1 from
  its markers; progenitor gets the highest panel loading (2.56) — fetal-reference
  composition dominates. None of these change Finding 2.

## Process lessons (repo-level, cheap fixes)

1. Hard-fail guards: abort (exit nonzero, write FAILED status) if query-to-matrix
   gene mapping hits < some fraction (e.g., 50%) of query genes, and if the score
   distribution is degenerate (zero variance) before any gate adjudication.
   Config lists gene_info; either use it or delete the key.
2. Commit small summary outputs (gate_result.json, positive_controls.csv,
   deconvolution_summary.json) or their hashes, so provenance doesn't rest on
   mtimes of untracked files.
3. The runner's instinct was right in shape (report FAIL honestly, no re-tuning)
   but wrong in kind: a crashed pipeline must be labeled INVALID/ERROR, never
   fed into the pre-committed decision rule. Pre-committed gates need a
   validity precondition clause.

## Final call per lessons L2/L3

**CONFIRM STREAM-KILL** for the CMap-reversal stream on ALS, with corrected
rationale: gate_result.json is void (overturned as evidence; the run is INVALID,
not FAILED), and the kill rests on the independent, pre-registered H-001 variance-
floor falsification plus absent neuronal composition signal. This is not post-hoc
rigging: nothing is being re-tuned to force a pass; the surviving rationale uses
only thresholds committed before any output existed, and explicitly forecloses
re-running the same design after a run.py patch (which cannot change stage B).
