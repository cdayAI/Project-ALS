# VERDICT: exp001b_deconvolved

**STATUS: STREAM KILLED (adversarial review 2026-08-25).**

Reviewed by: review-exp001b-2 (`reviews/exp001b_review.md`, branch `reviews/exp001b`).

CONFIRM STREAM-KILL — with corrected rationale:

1. The claimed gate FAIL is VOID as evidence. The scoring stage silently
   crashed: the query mapped to 0/15,000 L1000 genes (Entrez IDs vs HGNC
   symbols; `lincs.gene_info` never loaded, no abort guard), so every score was
   exactly -0.0 and every empirical p exactly 1.0. gate_result.json's FAIL is
   overturned; this run is INVALID, not FAILED.
2. The stream still dies on an independent PRE-REGISTERED ground: the H-001
   variance floor failed (overall R2 = 0.0885 < min_variance_explained 0.10;
   independently reproduced), and no ALS-vs-control neuronal composition shift
   exists (neuron delta +0.036, p = 0.86). The deconvolved-neuronal reversal
   query has no premise.
3. No re-tuning and no same-design retry after a run.py patch can change this:
   the fatal negative is in stage B, upstream of drug scoring. Revival requires
   a NEW pre-registered hypothesis (e.g., cell-type-pure query source).

Per AGENTS.md rule 4/5: no candidate advances; CMap-reversal stream for ALS is
killed without re-tuning (on_fail clause honored, different justification).

## What this run does (pre-registered)
1. Deconvolves the GSE124439 bulk ALS-vs-control signature into cell-type
   components using GSE221692 human spinal cord scRNA markers (H-001), and
   builds the reversal query from the NEURONAL component only.
2. Ranks 19,811 LINCS L1000 trt_cp perturbagens by MEDIAN cosine reversal,
   calibrated against a 1000-permutation gene-label-permutation null of the
   query, stratified by exposure (n_signatures); BH-FDR across all drugs.
3. Positive-control gate was committed in config.yaml BEFORE scoring.
4. Replication vs GSE255602 (independent iPSC-motor-neuron signature).

## Outcome summary
See outputs/ and README.md. The gate result (pass/fail) is reported there;
per the pre-committed `on_fail` clause a failed gate kills the stream.

