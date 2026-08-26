# VERDICT: PENDING

Awaiting adversarial review per AGENTS.md rule 4. Pre-registered criteria live in
config.yaml (restating hypotheses/H-007b.md, committed at 472e32c BEFORE scoring).

## Stage status

UPDATE (post count-level FC-B rerun): FC-B FAILS on raw counts too (rescue indices
-0.014 / -0.024, 0/2 disease genotypes). H-007b marked killed per its own pre-registered
FC-A (fold 1.578 < 2.0 on clean counts). LINCS scoring stages REMAIN PAUSED - the
reversal target module no longer meets its own replication bar. Pipeline code
(pipelines/perturbation_signatures/) is dataset-agnostic and reusable if any future
hypothesis survives module validation.

- FC-A module enrichment: PASSED at pre-registration (both datasets; see exp002).
- FC-B ropinirole treatment response: RUN - FAIL (0/2 disease genotypes rescue-positive;
  rescue indices negative: TDP43M337V -0.112, TDP43M337V+DRD2KO -0.101). Ropinirole does
  NOT restore the cytoskeletal/ECM module toward control within 6-84h ex-vivo exposure.
  Reviewer note: literature rescue efficacy (PMID 40826812) was for cell death/splicing/
  mitochondria, not necessarily this module's transcriptional state.
- FC-C translational bar: pending annotation of top-50 after scoring completes.
- FC-D cross-dataset consistency positive control: pending scoring completion.

## Provenance

LINCS GCTX read read-only from data/lincs/ (owner: sprint1-repurposing agent); scorer
vendored from exp001-sprint1-handoff@fcfff64 with attribution in lincs_score.py.

---

## FINAL CLOSURE (2026-08-26): KILLED-SUPERSEDED

review-h007b-1 confirmed the H-007b kill via independent raw-data reproduction
(reviews/h007b_review.md). With the reversal target module dead:

- Scoring stages were never run (paused at gate failure); no Level5 scores exist.
- This experiment is closed KILLED-SUPERSEDED. Parent claim released.
- PRESERVED for future hypotheses: pipelines/perturbation_signatures/lincs_score.py
  (streaming multi-query Level5 -cosine scorer, vendored from
  exp001-sprint1-handoff@fcfff64 with attribution) and the pre-registration pattern in
  config.yaml (incl. FC-D cross-dataset consistency positive control). Any future
  hypothesis must survive module validation BEFORE this pipeline runs again.
