# VERDICT: PENDING

Status per AGENTS.md rule 4: awaiting adversarial review (`reviews/exp002_c9orf72_module.md`).
No self-certification. Evidence below is what a reviewer must weigh.

## Result summary (run.py, single command, reproducible)

- Dataset: GSE303931 (isogenic C9orf72 iPSC-derived cortical neurons, PMID 42221822).
  2 isogenic pairs (C929, C952) x {Mutant, IsoControl} x 3 reps = 12 samples.
- DE (log2(TPM+1), OLS y ~ condition + line, BH FDR): 32 genes at FDR<0.05,
  69 at FDR<0.10, 246 at FDR<0.25 out of 18,351 tested expressed genes.
- Positive control (AGENTS.md rule 5): pipeline recovered Y-linked UTY
  (+3.06 log2FC, FDR 4.7e-05) and TMSB4Y (+2.32, FDR 1.3e-04) as top hits - known
  truth for a sex/Y-dosage contrast. CAVEAT: isogenic pairs should be sex-matched;
  this signal appearing INSIDE the mutant-vs-control contrast flags a culture-composition
  or annotation confound in the dataset. Reviewer must treat top hits accordingly.

## H-007 falsification criterion 1 - EVIDENCE TRIGGERED

Criterion: module needs >=2-fold enrichment of ribosome-biogenesis/nucleolar AND
speckle/splicing genes vs matched random sets (permutation P <= 0.05), else "the module
is not a faithful DPR correlate".

Observed (results/enrichment.csv, GO goa_human + obo descendant closure, 20k perms,
seed 20260825):
- nucleolar+ribo_biogenesis: fold 0.00 (FDR<0.10 module) / 0.69 (top-500|t| module)
- speckle+splicing:          fold 0.91 (FDR<0.10 module) / 0.50 (top-500|t| module)
- DPR curated literature set: fold 2.17 but P=0.11 (not significant)

Direction is DEPLETION, not enrichment, for both required gene-set families.
Consistent with the source paper itself: its headline finding is cytoskeletal/ECM/synaptic
(FLNB exon-30 skipping), not nucleolar stress (Series summary, GSE303931).

## What would change the verdict (for reviewer)

1. Replication attempt on GSE283507 (C9/FUS/TARDBP/SOD1 iPSC-MN panel) before killing:
   criterion 1 names iPSC-neuron data generally; one dataset may be underpowered (n=12,
   2 line pairs) or confounded by the Y-signal issue above.
2. Alternative module definitions (e.g., splicing-level rather than gene-level, since the
   series also reports splicing alterations) - not tested here.
3. If reviewer confirms criterion 1 as met across datasets: H-007 killed per its own
   pre-registered bar; the LINCS reversal arm (criterion 3) becomes moot.

## LINCS cross-ref status

Metadata-only cross-ref ran (results/lincs_crossref_metadata_only.csv): Level5 metadata IS
present locally. Module coverage of LINCS space is poor (69-gene module: 34/12328 measured,
only 2/978 landmark). Full connectivity/reversal scoring requires reading the Level5 GCTX
matrix owned by another agent via pipelines/perturbation_signatures/ - dependency noted,
NOT run here (no re-download performed).
