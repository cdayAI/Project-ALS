# MDE / Power Audit - H-001..H-010 (retroactive, per _LESSONS.md L8)

Audit agent: c9orf72-factory. Branch: hypotheses/mde-audit.
THIS FILE IS COMMITTED IN TWO STAGES: criteria first (this commit), then computed
per-hypothesis results (next commit). Git history proves the ordering.

## Fixed method (set before seeing any per-hypothesis numbers)

1. Variance priors (measured, not assumed where possible):
   - iPSC-neuron bulk RNA-seq residual SD (log2 TPM, OLS): Q25/median/Q75 = 0.18/0.27/0.41
     - measured on GSE303931 (18,351 genes, df=9 residuals)
   - postmortem cortex residual SD: NOT locally measured; audited at a documented planning
     value of 0.60 log2 units (typical bulk-tissue range 0.4-0.8). Any hypothesis relying on
     this prior must re-measure SD at analysis time and re-check its gate.
2. MDE formulas (80% power, two-sided alpha=0.05):
   - two-group per-gene: (t_{0.975,df}+t_{0.80,df}) * SD * sqrt(1/n1+1/n2), df=n1+n2-4
   - module aggregate (m genes): same but SD*sqrt((1+(m-1)*rho)/m), rho in [0.1, 0.3]
   - continuous predictor (repeat length, progression slope): detectable standardized beta
     = (t_{0.975,n-2}+t_{0.80,n-2})/sqrt(n)
   - enrichment/rank tests: empirical null scale reported alongside (measured null fold
     p95 = 1.16-1.46 for gene-set overlap tests; compound-rank nulls to be calibrated at run)
3. Dataset sizes used (actual, verified):
   - GSE124439 postmortem cortex: 176 samples = 145 ALS-spectrum vs 17 controls (GEO characteristics)
   - GSE255602 iPSC neurons: 71 samples = 27 ALS vs 44 control (local inventory metadata)
   - GSE261875 iPSC neurons: 32 samples = 16 TDP-43 vs 16 control (local inventory metadata)
   - GSE303931: 12 (6v6); GSE283507: 54 total / 18 in the DMSO isogenic contrast
   - LINCS Level5: 473,647 signatures; 205,034 trt_cp signatures covering 19,811 unique
     compounds across 76 cell lines; 216,081 trt_sh signatures targeting 8,178 genes
   - AnswerALS portal: ACCESS-BLOCKED (see data/inventory/answerals_access.md);
     population facts only (>=41 C9 carriers among 830 WGS'd)
   - PRO-ACT: open after registration; >10,000 subjects clinically annotated; NfL lab
     subset size UNVERIFIED until download (published analyses suggest low hundreds)
4. Verdict rules (applied mechanically):
   - TESTABLE-NOW: every dataset needed for the PRIMARY endpoint is available or open with
     routine registration, AND the primary endpoint's MDE at available n <= the effect size
     the hypothesis needs to matter (or its bar is expressed in units power reaches).
   - NEEDS-BIGGER-DATA: pipeline feasible but a required dataset is access-blocked, too
     small for the endpoint's MDE, or an unverified subset size gates the claim.
   - NEEDS-REDESIGN: as written the hypothesis has no quantifiable falsifiable endpoint,
     depends on artifact-prone inputs without a counts-level fallback, or its bar cannot be
     cleared by any realistic n - rewriting criteria (with pre-registration) is required.
   - Contingent qualifiers appended when only SECONDARY arms are testable now.

## Per-hypothesis results

(computed in next commit - intentionally absent here)


## Per-hypothesis results (computed AFTER the criteria commit above)

### H-001 - cell-type-deconvolved signatures re-rank exp001 hits
- Data: GSE124439 (145 ALS vs 17 controls, postmortem cortex; locally present), open marker
  matrices, LINCS stack local.
- MDE: per-gene |log2FC| = 0.43 at planning SD 0.60 with unbalanced n (145/17);
  300-gene module aggregate MDE = 0.195. Deconvolution is a linear-unmixing step, not an
  extra data need. Re-ranking endpoint is rank-stability over 19,811 compounds - power comes
  from compound count, ample.
- Caveat: control n=17 limits subtype resolution within controls; postmortem covariates
  (age/PMI/RIN) must be modeled and will inflate the 0.60 prior.
- VERDICT: TESTABLE-NOW.

### H-002 - immune-dominant reversal -> YKL-40/CRP PD footprint convergence
- Data: inherits H-001 outputs + sig_info + MSigDB/Reactome sets. All local/open.
- Power: enrichment/convergence test over ~19,811 compounds; measured gene-set null p95
  1.16-1.46 means a pre-registered convergence bar must sit >=2-fold on set-overlap or use
  rank-based statistics calibrated on the compound-label permutation null.
- Contingent: blocked behind H-001 completion (pipeline dependency, not data).
- VERDICT: TESTABLE-NOW (contingent on H-001).

### H-003 - TBK1 causal colocalization + cell-type direction-of-effect
- Data: GWAS Catalog sumstats (open), GTEx v8 brain eQTL (open), MinE summary stats (DAC
  PENDING), AnswerALS expression (PORTAL-BLOCKED).
- Power: colocalization itself needs only fine-mapped credible sets (available for large ALS
  GWAS); the direction-of-effect claim needs cell-type-resolved expression/QTL - AnswerALS-
  gated. GTEx cortex n~100-190 gives eQTL power for common variants; adequate for coloc arm.
- VERDICT: NEEDS-BIGGER-DATA (primary endpoint blocked on MinE DAC + portal).
  Preliminary arm (GWAS x GTEx coloc, no cell-type split): TESTABLE-NOW.

### H-004 - UNC13A cryptic-exon sQTL colocalization + stratification marker
- Data: GWAS Catalog (open), brain sQTL resources (open), AnswerALS RNA-seq w/ genotype
  (BLOCKED) for the cryptic-exon patient-stratification endpoint.
- Power: sQTL coloc feasible on open data; stratification-marker claim requires donor-level
  splicing outliers - portal-gated. Note: splicing-load readout design now exists in H-012;
  consider merging endpoints there to avoid duplication.
- VERDICT: NEEDS-BIGGER-DATA (stratification arm). Coloc arm: TESTABLE-NOW.

### H-005 - pTDP-43 kinase/phosphatase modifier nomination via perturbation-signature reversal
- Data: module derivation needs AnswerALS RNA-seq/proteomics (PORTAL-BLOCKED). Reversal layer
  ready: LINCS trt_sh covers 8,178 targeted genes / 216k signatures; trt_cp 19,811 compounds.
- Power of reversal layer: adequate by construction (thousands of perturbagens); but the
  DERIVATION side is fully blocked, and kinase-reversal inference from knockdown signatures
  needs its own validation arm before results mean anything.
- VERDICT: NEEDS-BIGGER-DATA (derivation dataset).

### H-006 - KIF5A/DCTN1 transport module predicts fast progression
- Data: AnswerALS multi-omics + progression labels (BLOCKED); GWAS overlay open.
- Power: "predicts fast vs slow progression" needs donor-level omics-clinical pairing -
  exists only on the portal. With ~1,000 clinical subjects but far fewer multi-omics donors,
  expect low hundreds at best; detectable AUC increment >=0.05-0.10 requires n>=150-200
  donors with outcome labels (standard ROC comparison bounds). Unverifiable until access.
- VERDICT: NEEDS-BIGGER-DATA.

### H-008 - NfL-trajectory ML trial enrichment (PRO-ACT)
- Data: PRO-ACT open after registration (>10k clinically annotated subjects); NfL lab subset
  UNVERIFIED - published PRO-ACT NfL analyses operate in the low hundreds.
- Power: trial-enrichment simulations with n in the thousands are well-powered to detect the
  claimed >=25% sample-size reduction as a difference in simulated required-n distributions;
  binding constraint is NfL coverage per subject (baseline + >=1 follow-up needed for slope).
  Gate: if NfL-annotated subjects <200 at download, downgrade to simulation-on-published-
  kinetics only and mark the empirical arm NEEDS-BIGGER-DATA.
- VERDICT: TESTABLE-NOW (with the stated NfL-coverage gate at download).

### H-009 - ML molecular subtypes aligning with clinical progression strata
- Data: subtype training needs AnswerALS expression (BLOCKED); clinical strata via PRO-ACT
  (open); blood-marker panel from PMID 38129934 (literature).
- Power: alignment testing needs both sides per donor - portal-gated. Clinical-only
  clustering on PRO-ACT is possible but does not test the transcriptomic-alignment claim.
- VERDICT: NEEDS-BIGGER-DATA.

### H-010 - aggregated neural-weighted LINCS reversal scores pass exp001 positive-control recalibration
- Data: ALL LOCAL (Level5 gctx 473,647 sigs, sig_info, exp001 positive_controls list).
- Power: stabilization claim measurable directly - rank-correlation/variance reduction across
  contexts with thousands of compounds; positive control defined already in exp001 config.
  No access blockers, no sample-size limits beyond compute (~23 GB matrix streaming, proven
  pattern in lincs_score.py).
- VERDICT: TESTABLE-NOW.

## Triage summary

| Hypothesis | Verdict | Binding constraint |
|---|---|---|
| H-001 | TESTABLE-NOW | none (145v17 bulk imbalance noted) |
| H-002 | TESTABLE-NOW* | contingent on H-001 outputs |
| H-003 | NEEDS-BIGGER-DATA | MinE DAC pending + portal blocked (coloc arm NOW) |
| H-004 | NEEDS-BIGGER-DATA | portal blocked for stratification arm (coloc arm NOW) |
| H-005 | NEEDS-BIGGER-DATA | portal blocked for module derivation |
| H-006 | NEEDS-BIGGER-DATA | portal blocked (omics+progression pairing) |
| H-007 | KILLED | excluded from audit |
| H-007b | KILLED | excluded from audit |
| H-008 | TESTABLE-NOW | NfL-coverage gate (>=200 subjects) at download |
| H-009 | NEEDS-BIGGER-DATA | portal blocked (subtype training) |
| H-010 | TESTABLE-NOW | none (all local) |

Sprint-capacity implication: immediate compute value sits in H-010, H-001 (+H-002 chained),
and the preliminary coloc arms of H-003/H-004. Five hypotheses unblock on exactly two keys:
the AnswerALS portal credentials and the MinE DAC application - start both early.
