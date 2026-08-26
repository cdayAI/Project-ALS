# Adversarial review: exp002_c9orf72_module (H-007)

Reviewer: review-exp002-1 (independent)
Date: 2026-08-25
Branch under review: exp/h007-c9orf72-module (merged to main at 4b63a9a)
Claimed outcome under test: H-007 KILLED; H-007b drafted.

## 1. Was the kill pre-registered and is it justified?

**Pre-registration: YES (weak-form but genuine).**
- `hypotheses/H-007.md` containing the falsification criterion ">=2-fold enrichment of
  ribosome biogenesis/nucleolar AND speckle/splicing genes vs matched random gene sets,
  permutation P <= 0.05, else killed" was committed at 19:39 (4797131, hypothesis-drafter
  agent). The first exp002 results commit (12fdebd, config + enrichment.csv + PENDING
  verdict) is 20:06 the same day. The criterion text in config.yaml matches H-007.md
  verbatim; `git diff 12fdebd..472e32c -- config.yaml` shows only the replication-arm block
  was added later - the bar itself was never moved.
- Caveat: this is same-session pre-registration (~27 min lead), not an external registry.
  However, the rigging incentive runs the wrong way: a kill criterion written to be passed
  would not have been set at 2-fold when the data deliver ~1.0. Accepted as genuine.
- The criterion as triggered: every nucleolar/ribo and speckle/splicing row in
  results/enrichment.csv shows fold <= 0.91 with non-significant (or depletion-consistent)
  permutation P. FC1's kill condition (< 2-fold AND P > 0.05) is met on all rows.

## 2. Statistical validity of the GSE303931 arm (n=12, two line pairs)

- Design `y ~ condition + line` (OLS, df=9) is estimable but minimal: no line x condition
  interaction, no surrogate-variable/composition covariates. With 2 pairs, any line-level
  confound (differentiation batch, culture composition) is only partially controlled.
- The Y-linked positive control (UTY log2FC +3.06 FDR 4.7e-05; TMSB4Y +2.32 FDR 1.3e-04,
  top hits in BOTH pairs' contrasts) proves the mutant cultures carry more male-cell signal
  than their own isogenic controls - a sex/culture-composition imbalance that should not
  exist between truly isogenic matched cultures. This contaminates parts of the DE list
  (Y/X-dosage genes, maturity-composition correlates).
- Does it undermine the kill? **No.** A composition confound injects spurious sex/maturation
  genes into the DE list; it dilutes real signal but cannot manufacture systematic DEPLETION
  of nucleolar/ribo and speckle genes. Depletion is observed across every module definition,
  including top500|t| (fold 0.69 / 0.50), which does not depend on BH thresholds. Removing
  confounded genes could not rescue >=2-fold nucleolar enrichment from fold 0.0-0.9.
  Residual uncertainty: n=12 means the DE list is noisy overall; the kill is valid for the
  module AS COMPUTED from this dataset (which is what FC1 tests), not as proof that C9orf72
  neurons lack nucleolar-stress biology.

## 3. GSE283507 replication arm - materially flawed; treat as directional-only

Three independent problems:
1. **Near-zero within-group variance** (agent-disclosed, median within-group log2TPM var
   ~0.019): all t/P are inflated artifacts; the FDR<0.10 "module" contains 3,279 genes -
   essentially everything tested. All permutation P values in this arm are meaningless.
2. **The protocol was NOT identical**, contradicting VERDICT.md. Implied background from
   the committed CSVs (N = m_set * k / expected_overlap) is ~4,484 genes here vs 16,253 in
   GSE303931. GO sets also collapsed: nucleolar+ribo 754 -> 138 genes, speckle+splicing
   516 -> 89, DPR-curated 30 -> 5 (only 5 of 30 curated genes present!). With m=138 inside
   a k=3,279-gene module, fold enrichment is arithmetically pinned near 1.0 regardless of
   biology.
3. **No code committed for this arm.** Commit 472e32c adds only the results CSV; no script
   reproduces gse283507_replication_enrichment.csv (AGENTS.md rule 6 violation).
4. (Minor) The contrast is TARDBP M337V vs TDP43+/+ - cross-genotype by construction, so it
   never tested C9orf72/DPR specificity anyway.

Should it be excluded entirely? Functionally yes: its numbers are uninterpretable, and the
kill does not need it. Excluding it leaves a **one-dataset kill on the designated discovery
dataset GSE303931, which is decisive by itself** (see section 2). 'Directional evidence
only' understates the problem; I treat the arm as non-evidentiary for the verdict, while
noting the direction was consistent with failure.

## 4. H-007b quality

- FLNB anchoring PMIDs verified via NCBI E-utilities: PMID 42221822 exists ("Global
  transcriptional changes across multiple isogenic C9orf72 patient iPSC-derived neurons",
  iScience 2026) and is the ONLY PubMed hit for FLNB+C9orf72 / FLNB exon-30 skipping;
  PMID 40826812 (ropinirole ALS, J Neurochem 2025) checks out. Anchoring is legitimate.
- Falsifiability concerns:
  a. The bar drops from 2-fold (H-007) to 1.5-fold, AND becomes disjunctive over four large
     gene sets ("cytoskeleton ... ECM ... synapse OR cell adhesion"). Large sets (731 genes)
     plus a lowered disjunctive bar = easy to pass. Post-hoc threshold-setting after seeing
     fold 1.09/1.66 data. This part reads as rigged toward survival.
  b. FC2 (ropinirole module-reduction positive control) rests on the SAME artifact-laden
     GSE283507 TPM matrix flagged above. Module-score changes there are not interpretable
     until counts are reprocessed or a clean matrix is obtained.
  c. Evidence line cites fold 1.09 P=0.0001 as "replication" - a 9% excess overlap from a
     broken pipeline is not evidence of anything.
- Salvage path: require the out-of-sample third dataset FIRST, recompute everything from
  raw counts with the full 16k background, keep the >=2-fold standard (or justify 1.5-fold
  prospectively against null calibration), and make FC1 conjunctive or effect-size-based
  rather than disjunctive. As drafted, H-007b is NEEDS-DATA quality, not open-and-testable.

## 5. Process notes

- VERDICT.md was moved from PENDING to KILLED before formal review (parent-agent
  authorization documented). Technically an AGENTS.md rule-4 breach, but in the safe
  direction (kill, not self-certification of survival) with provenance noted. Accepted.
- Positive control handling is honest: Y-gene recovery is reported as a dataset confound
  flag, not spun as validation.

## VERDICTS

**[H-007]: KILLED - CONFIRMED.** Pre-registered criterion genuinely triggered on the
designated discovery dataset (GSE303931, depleted across all module definitions); the flawed
replication arm is unnecessary for the kill.

**[exp002 experiment]: KILL ACCEPTED, EXECUTION DEFECTS LOGGED.** Result stands; required
remediation before reuse of either arm: (1) commit the GSE283507 analysis code, (2) resolve
the 4.5k-vs-16k background and GO-set-collapse discrepancy, (3) reprocess GSE283507 from raw
counts before any downstream use (including H-007b FC2).

**[H-007b]: NOT APPROVED AS DRAFTED - downgrade to NEEDS-DATA** pending out-of-sample
dataset, count-level reprocessing, and de-rigged (conjunctive/effect-size) criteria.
