# Adversarial review: H-007b remediation + kill confirmation (exp002 addendum)

Reviewer: review-h007b-1 (independent)
Date: 2026-08-26
Under review:
- commits 5d64b58 (H-007b revision + remediation code, 2026-08-25 21:16:25 EDT)
  and b218517 (count-level rerun results + VERDICT addendum, 21:18:35 EDT) on
  `exp/h007-c9orf72-module` (NOT yet merged to main at time of review)
- prior review forcing this remediation: reviews/exp002_review.md

## SCOPE 1 - Is the H-007b KILL valid?

### 1a. Was the remediation sequence genuinely pre-analytic?

**Commit order verified: YES.** 5d64b58 contains the revised falsification criteria
(conjunctive single-primary-family: DOWN-regulated DE x cytoskeleton organization
GO:0007010 only, >=2-fold AND perm P<=0.05, full gene universe, out-of-sample leg (b),
null calibration reported) plus the committed analysis code replication_gse283507.py;
b218517 adds only results + status flip. Criteria text and code gating logic match
(`fc_a = fold >= 2.0 and perm_p <= 0.05` evaluated solely on the primary row).

**Caveat, logged:** the two commits are 2m10s apart. Parsing goa_human.gaf.gz,
running the OLS DE pass, and 20k permutations over ~20 module x family combos in under
two minutes is implausible; the run most likely executed before 5d64b58 was committed.
The "pre-analysis" is therefore weak-form in letter. It is nonetheless conservative in
direction: of the three significantly enriched secondary families, synapse (2.037)
would have PASSED the 2.0 bar while the chosen primary (cytoskeleton organization,
1.578) FAILS. A post-hoc selection made after seeing results picked the failing family
as sole gate. Rigging toward survival is absent; and section 1c shows the numbers
reproduce independently, so the verdict does not depend on resolving this ambiguity.

### 1b. Was the conjunctive criterion applied as written?

YES. Gate row in results/gse283507_countlevel_enrichment.csv:
down_FDR10 x primary_cytoskeleton_organization: k=2069, m=305, overlap 121,
expected 76.686, **fold 1.578**, perm P=5e-05 (permutation floor at 20k draws).
1.578 < 2.0 -> criterion leg (a) FAILS; conjunctive design makes leg (b)
(out-of-sample) moot. Internal arithmetic checks out exactly
(2069*305/8229 = 76.69; 121/76.686 = 1.578).

### 1c. Independent reproduction (fresh GEO download, different GO source)

I re-downloaded GSE283507_raw_Count_FPKM_TPM.csv.gz from GEO and re-implemented the
committed pipeline (median-of-ratios size factors, log2, genotype+time OLS, BH FDR):
- universe after symbol collapse: 8229 (matches); tested genes 6130 (matches de.csv);
- down_FDR10 = 2069, up_FDR10 = 2265 (both match);
- median within-group log2 normalized-count variance = 0.439 (claimed 0.44);
  deposited-TPM equivalent by my definition = 0.055 vs their reported ~0.019 -
  exact figure differs by computation detail but the conclusion holds: the artifact
  TPM matrix carries near-zero within-group variance and the count-level rebuild
  restores it (~8x higher here). Variance-artifact claim CONFIRMED.
- Enrichment re-checked against QuickGO human annotations with is_a-only descendant
  closure (independent of their gaf snapshot): m=311, observed overlap 120,
  **fold 1.535** vs their 1.578 - same conclusion, small delta attributable to
  annotation snapshot date/symbol mapping. Hypergeometric z ~= 5.2, so the P<=0.05
  leg passes easily; it is the FOLD leg that fails. Fold 1.55-1.6 sits far above the
  empirical null (p95 ~1.16) but clearly below the pre-registered 2.0 standard that
  killed H-007. Applying any lower bar post hoc would be the exact error the prior
  review prohibited.

### SCOPE 1 VERDICT: KILL CONFIRMED.
Criterion genuinely fails at leg (a); variance artifact resolved; result reproduced
independently end-to-end.

## SCOPE 2 - The secondary observation (ECM 1.94 / synapse 2.04 / adhesion 1.85, all P=5e-05)

Is the agent's refusal to rescue H-007b with these justified? **YES.**
Switching the primary family to synapse AFTER seeing that the pre-chosen primary fails
is precisely the disjunctive/moving-bar rigging failure mode flagged in
reviews/exp002_review.md sec 4a. Four correlated families were tested; taking the max
requires max-statistic calibration, and a nominal 2.037 from a post-hoc selection is
not admissible evidence. Additionally, FC-B failed independently (rescue index
R = -0.014 / -0.024 in BOTH disease genotypes at count level): the
treatment-reversal premise - the actual deliverable of H-007b - is dead for the
cytoskeleton module, so even an FC-A pass would not have yielded a live hypothesis.

**However: a narrow follow-up IS defensible, hence DRAFT-NOTE.** The synapse signal is
not noise: fold 2.037 vs null p95 1.134 (P at permutation floor), coherent
direction-specificity (up-module synapse DEPLETED at 0.728; top500|t| synapse 2.254),
and an existing literature anchor (synaptic signaling alterations in PMID 42221822).
A draft H-007c may be written IF it: (i) declares its own post-hoc origin explicitly
(family selected from GSE283507 count-level screen); (ii) treats ALL of GSE283507 as
discovery-only and gates survival SOLELY on >=2-fold + P<=0.05 for synapse
(GO:0045202 descendants) among DOWN-regulated DE genes in an untouched out-of-sample
dataset (GSE284339 or comparable), one shot, no interim peeking; (iii) specifies a FRESH
treatment-response readout (ropinirole response must be re-tested on a synapse module -
FC-B failure does not automatically transfer across modules, but neither can success be
assumed); (iv) keeps the 2.0-fold standard unchanged. That design cannot repeat the
post-hoc rigging mistake because the confirmatory test is fully out-of-sample.

## Process defects to fix (non-blocking for the verdict)

1. **Remediation is unmerged.** 5d64b58/b218517 exist only on `exp/h007-c9orf72-module`;
   origin/main still shows H-007b status "open", a VERDICT.md without the addendum, and
   none of the count-level result CSVs. Merge the branch so the kill record exists where
   the ledger points.
2. Commit-message provenance claim ("criteria ... BEFORE count-level recompute") is
   overstated given the 2-minute window; keep such claims tied to what commit history
   can actually demonstrate.
3. Minor: reconcile the reported artifact-TPM variance (~0.019) with the computation
   used (I get 0.055 under a plain all-genes median); state the definition next to the
   number.

## LEDGER STATUS CONFIRMATION

AGENTS.md ledger entry updated: H-007b KILLED per own pre-registered FC-A/FC-B
(count-level GSE283507 leg failed; reversal premise failed 0/2), kill CONFIRMED by this
review; DRAFT-NOTE recorded for optional H-007c synapse-primary out-of-sample proposal
under the conditions above.

## VERDICT

**CONFIRM KILL + DRAFT-NOTE for follow-up hypothesis.**
