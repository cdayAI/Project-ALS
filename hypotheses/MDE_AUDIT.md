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
