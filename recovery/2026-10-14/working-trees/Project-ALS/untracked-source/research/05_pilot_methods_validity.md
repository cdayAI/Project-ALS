# Pilot Sprint #1 — Statistical and Validation Contract

**Status:** retrospective implementation audit plus prospective retry contract, current through 2026-08-25.

**Repository baseline reviewed:** `origin/main` at `84f0867`.

**Latest remote delta reviewed:** `origin/main` at `e734295` (formal exp001, H-007, and H-007b kills; H-008 simulation; reviewed trial matcher and causal-target scaffold; unreviewed data-QC, structure-readiness, and mocked AnswerALS gate-check merges; H-007c/H-011/H-012 stubs; exp003's appended `KILLED-SUPERSEDED` closure; and two later heartbeat-only status commits). Active exp001 handoff reviewed at `origin/exp001-sprint1-handoff@fcfff64`.

**Interpretation status:** **HOLD**. No compound rank from the reviewed implementation is evidence of ALS efficacy, and this document does not report a drug candidate or an experiment result.

This contract converts the broad strategy in `research/02_ai_methods.md` into falsifiable rules for the transcriptome-reversal pilot. It is deliberately stricter than “run differential expression, find an anti-correlated drug, and inspect the top names.” Transcriptome reversal is allowed here only as a low-cost hypothesis generator. Candidate interpretation begins only after the unit of analysis, null model, controls, replication, context, and exposure gates below pass.

---

## 1. Executive decision

A computational screen is technically feasible in principle, but the implementation reviewed at `84f0867` is **not ready for candidate interpretation**. The primary blockers are:

1. The code cannot complete as configured: the GEO parser calls `startswith` on a list of sample-characteristic values; `h5py` is asked to open a gzip-compressed `.gctx.gz` file directly; and the GCTX metadata paths and matrix orientation do not match the official schema.
2. Public records list many RNA libraries but substantially fewer independent donors. Treating libraries or tissues as independent people creates pseudoreplication and anticonservative uncertainty estimates.
3. Plain log-CPM ordinary least squares with a custom variance shrinker is not equivalent to limma-voom and does not account for repeated donors, region, or other measured covariates.
4. The reviewed ranking direction is inverted. `score_reversal = -cosine`, so larger positive values represent stronger reversal, but signatures are ranked ascending and compounds are ordered by the minimum score. The current “top” compounds therefore favor disease mimics, not reversers.
5. Cosine similarity across all 12,328 L1000 rows is not the published CMap bidirectional WTCS/NCS/tau query and makes measured landmark genes indistinguishable from inferred genes.
6. Selecting one best signature per compound rewards compounds tested in more cells, doses, and times and makes an outlier the primary result, even after the score direction is corrected.
7. The proposed “positive-control” panel is not a set of known ALS-active drugs. It includes drugs with negative ALS trials and minocycline, which produced significantly faster ALSFRS-R decline than placebo in a Phase III trial.
8. Molecular weight and hydrogen-bond-donor counts are coarse chemistry descriptors, not evidence of blood-brain-barrier penetration or achievable CNS exposure.
9. Independent-cohort replication and end-to-end empirical nulls are planned in prose but are not part of the reviewed rank-generating path.

These are stop gates, not cosmetic improvements. A technically completed 473,647-signature screen that fails them remains an exploratory code output.

The score-direction defect is visible directly in `experiments/exp001_repurposing/run.py`: line 188 negates cosine and labels it reversal strength, while lines 362–369 rank ascending and aggregate the minimum. Fixing only the sort direction would still leave the noncanonical score, row-alignment assumption, context bias, and multiplicity problems below.

### Current exp001 disposition

Main now marks exp001 **KILLED**, and that disposition is correct: no drug from the run advances. The merged review identifies exposure/min-of-N bias, failed GSE255602 top-hit replication, query confounding, and annotation failures. It does not, however, detect the more basic score-order defect above and incorrectly says the scorer is correct. Consequently, the review's percentiles, exposure matching, and candidate comparisons all operate on the wrong biological direction: with `score_reversal = -cosine`, an exact disease mimic scores `-1`, an exact inverse scores `+1`, and ascending/minimum ranking selects the mimic.

Main still contains the original non-running implementation. Calling its committed GEO parser on the committed `series_matrix.txt` reproduces `AttributeError: 'list' object has no attribute 'startswith'` at `run.py:51`; its GCTX path also retains the compressed-file, schema, orientation, and undefined-`h5py` defects. The handoff branch repairs parsing, decompression, metadata paths, orientation, and gene-ID joining, but retains ascending/minimum ranking at lines 428–434. All 30 rows in its annotated top-candidate snapshot have negative best scores (`-0.388` to `-0.348`), so they are the strongest disease-direction similarities under the documented formula, not reversers.

The default branch contains no exp001 outputs. The handoff snapshot omits the full differential-expression table, all signature scores, the primary full ranking, sample sheet, and annotation cache; its committed code and annotation-table schema also disagree. The review's GSE255602 reanalysis code and outputs are not committed. A retry must first fix and unit-test direction and aggregation, then use donor-aware disease modeling, frozen null-calibrated context aggregation, valid controls, and committed replication evidence. Re-labeling the same scores cannot salvage a candidate.

### Live hypothesis-ledger preflight

The H-001–H-010 batch added to `origin/main` during this audit does not clear these gates:

- H-001 and H-010 inherit the invalid exp001 control doctrine, including minocycline as a required positive despite its harmful ALS trial result.
- H-006 assumes a thresholded ChEMBL “BBB score.” ChEMBL does not supply a validated BBB score that can establish human CNS exposure.
- H-010 proposes judging aggregation by whether it preserves the current control ranks. Controls must validate the assay and scoring implementation independently; preserving an invalid rank is not evidence of improvement.

Those hypotheses may remain open research questions only after their inputs and falsification criteria are corrected. They must not be treated as preregistered validation of the present pipeline.

The later dataset inventory improves file-level provenance by adding checksums and accession manifests, but it does not increase the number of independent biological units. It also calls the GSE255602 controls “isogenic”; the accession manifest instead contains five separate non-diseased control donor lines alongside four ALS donor lines. That label must be corrected before the inventory is used as a design specification.

Later heartbeats and the merged verdict correctly say the run is killed, but status claims still exceed the evidence preserved on main. The current default branch has no exp001 outputs, the executable/snapshot handoff is divergent and unmerged, and the review replication is uncommitted. Heartbeat claims of clean Git hygiene are also contradicted by the 181 ignored-but-tracked raw experimental blobs already in history. Status prose is not a substitute for code, immutable inputs, outputs, manifests, and review artifacts in one reproducible commit graph.

### Immediate code gates

These failures occur before any biological interpretation:

| Gate | Reviewed behavior | Required disposition |
|---|---|---|
| GEO metadata | `parse_series_matrix` iterates over characteristic-row lists and calls `startswith` on each list. Running it on the committed matrix raises `AttributeError: 'list' object has no attribute 'startswith'`. It also omits the subject ID needed for inference. | Parse the value at each sample index, retain subject ID, and assert a reconciled subject/sample/region manifest before DE. |
| GCTX transport and schema | The configured file is gzip-compressed, while `h5py.File` expects an HDF5 file. The code asks for `metadata/rows/id` and `metadata/cols/id`; the official GCTX nodes are `/0/DATA/0/matrix`, `/0/META/ROW/id`, and `/0/META/COL/id`. The official parser transposes the stored matrix. | Decompress, verify the official checksum, use [`cmapPy`](https://github.com/cmap/cmapPy/blob/master/cmapPy/pandasGEXpress/parse_gctx.py) or the exact schema, and assert dimensions and orientation before streaming. |
| Gene alignment | The code uses the row number of `gene_info` as the GCTX row number, then reads GCTX IDs only after scoring. | Join GCTX row IDs to `pr_gene_id`, require a one-to-one mapping, and assert landmark/BING counts and order before constructing a query. |
| Score sign | An exact disease mimic receives `score_reversal = -1`; an exact inverse receives `+1`; ascending/minimum ranking promotes the mimic. | Add a synthetic sign unit test before any real data and fail closed on a mismatch. |
| Stop behavior | `sanity_ok` is not used, positive controls do not gate execution, and a candidate table is written regardless. | Every technical and statistical gate in section 10 must terminate or explicitly downgrade the run before candidate output. |

The annotation path also requires preflight: the configured PubChem URL already contains `/name/properties/JSON`, then the code appends another name/property route. Identity lookups must pass fixed test cases before they annotate a ranking.

---

## 2. What the datasets actually represent

| Dataset | Publicly listed material | Independent biological unit | Valid role in the pilot | Critical limitation |
|---|---|---|---|---|
| [GSE124439](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124439) | 176 postmortem cortex transcriptomes from 95 subject IDs in the current series matrix | Donor, with region/tissue repeated within some donors | Human postmortem discovery and region-specific replication | The current matrix contains 145 ALS libraries from 76 subjects, 17 non-neurological-control libraries from 12 subjects, and 14 other-neurological-disease libraries from 7 subjects; 48 subjects span regions and 7 subject-region pairs have repeated libraries. The paper-level categories differ slightly, reporting 148 ALS transcriptomes from 77 ALS donors and 28 controls from 18 people (12 non-neurological and 6 neurological), so the manifest must reconcile the discrepancy ([Tam et al., 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6866666/), PMID 31665631, DOI 10.1016/j.celrep.2019.09.066). |
| [GSE255602](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE255602) | 71 bulk RNA-seq libraries spanning chip/well, flow/static, experiments, and replicates | iPSC donor line; pure technical replicates are subsamples | Limited human iPSC orthogonal support and cross-model sensitivity analysis | Although the GEO summary describes SC-chips generated from 8 sporadic ALS and 10 non-diseased donors in the broader study, the 71 samples listed in this bulk accession contain only 4 ALS and 5 control lines, repeatedly assayed across experiments, chips/wells, flow/static conditions, and replicates. The accession—not the paper-level cohort headline—defines the usable inferential `n`. GEO describes the processed gene matrix as Salmon/tximeta-derived pseudocounts; its scale, transcript-length/offset handling, and suitability for the chosen likelihood must be verified or the data reprocessed before a count model is fit ([Lall et al., 2025](https://pubmed.ncbi.nlm.nih.gov/40562033/), PMID 40562033, DOI 10.1016/j.stem.2025.05.015). |
| [GSE261875](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261875) | 32 bulk RNA-seq libraries: four replicates for each of eight TDP-43/ATXN2 conditions | Replicate within an isogenic motor-neuron system | Mechanistic and protective-intervention calibration | It is not an independent patient cohort. It tests EV/YFP/TDP-43/TDP-43ΔNLS with ATXN2 intact or knocked out in one isogenic system; GEO lists no linked peer-reviewed paper as of the review date. |
| [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) | 1,319,138 Phase I L1000 profiles; Level 5 is 473,647 × 12,328 | Aggregated perturbational signature, stratified by cell, dose, time, and perturbagen | Perturbation lookup and reversal hypothesis generation | Level 5 aggregates replicates. L1000 directly measures 978 landmark genes and infers the remainder; most contexts are non-neural and many are cancer-derived. GSE70138 is the related production-phase release. |

The sample sheet committed for GSE124439 uses group counts that differ slightly from the publication's paper-level categories. The run must preserve a versioned sample manifest and reconcile every inclusion/exclusion to the paper and GEO metadata before modeling. The inferential `n` reported in results must be the number of independent donors contributing to each contrast, never just the number of columns in a count matrix.

The same principle applies to single-cell follow-up: cells from one donor are not independent people. Donor-level pseudobulk or an appropriate mixed model is required; methods that treat cells as replicates can generate widespread false discoveries ([Zimmerman et al., 2021](https://www.nature.com/articles/s41467-021-21038-1); [Squair et al., 2021](https://www.nature.com/articles/s41467-021-25960-2), PMID 34584091).

---

## 3. Locked biological estimands

Do not merge raw counts across postmortem cortex, spinal-cord chips, isogenic perturbations, and single-cell assays. They answer different questions.

Predeclare these estimands before viewing any drug rankings:

1. **Primary human-tissue estimand:** ALS versus non-neurological control in frontal cortex (65 ALS libraries/63 subjects versus 9/9 controls), using one contribution per donor-region or an explicit repeated-donor model. Lateral motor cortex (37/35 versus 4/4) and medial motor cortex (38/37 versus 4/4) are underpowered directional sensitivity analyses, not standalone primary CMap-query generators, unless a score-blind minimum-information/stability simulation specified before DE shows adequate calibration. The generic “Motor Cortex” label has 5 ALS libraries from 4 subjects and no non-neurological controls, so that contrast is not estimable.
2. **Cross-region estimand:** compare the three region-specific coefficients and uncertainties or fit a prespecified repeated-donor region-interaction model. Do not pool the library columns and call the result a larger cohort; the very small control strata require correspondingly wide uncertainty and sensitivity analyses.
3. **iPSC estimand:** sporadic ALS versus control within a prespecified homogeneous chip/flow/experiment design, with donor line as the biological unit and clone/differentiation/technical replication modeled or collapsed. GSE255602 experiment 1.3 contains controls only, and the same 4 ALS plus 5 control lines recur across later experiments; those runs are repeated contexts, not independent donor replications.
4. **Mechanistic estimand:** TDP-43 perturbation and the TDP-43 × ATXN2-knockout interaction in GSE261875.
5. **Composition estimands:** unadjusted total-tissue disease effect and a separately labeled composition-adjusted sensitivity analysis. Gliosis and neuronal loss may be consequences of ALS, so composition adjustment changes the question and may remove real disease biology.
6. **Subtype estimands:** molecular or genetic subtypes only if defined independently of the CMap result and supported by enough independent donors. GSE124439 itself contains heterogeneous molecular subtypes, so a global average may erase or reverse subgroup biology.

For this protocol, an **independent replication cohort** means a separately sampled cohort with no overlapping donors or donor lines and a biologically comparable estimand. Different regions from the same GSE124439 donors and different experiments reusing the same GSE255602 lines are not independent replications. A postmortem ecosystem signature and an iPSC motor-neuron-intrinsic signature must remain separate even if some genes overlap; their concordance is useful orthogonal support, not replication of one common effect. An “ALS core signature” requires directionally consistent evidence in at least two independent donor cohorts measuring a comparable estimand.

---

## 4. Within-study differential-expression contract

### Required metadata and QC

Before fitting a model, produce a sample-level table containing, where available:

- donor ID and all repeated samples derived from that donor;
- diagnosis, molecular/genetic subtype, tissue and region;
- sex, age, ancestry information, disease duration, age at death;
- postmortem interval, RIN, brain bank, library chemistry, sequencing batch;
- for iPSC work: donor, clone, differentiation and assay batch, day/maturity, format, flow, and estimated cell composition;
- raw/processed file checksum, source URL, retrieval date, and annotation version.

Report missingness and case-control balance for every proposed covariate. No batch method can recover a disease effect if diagnosis is perfectly confounded with batch, region, or donor source. In that situation the affected contrast is non-identifiable and must be stopped or relabeled descriptive.

Minimum QC includes library size, expressed-gene count, count-distribution plots, sample correlation, PCA before and after justified covariate handling, sample-identity checks, and predeclared outlier rules. Do not remove a sample because its inclusion weakens a preferred drug rank.

### Primary statistical models

- For independent-donor count data, use a maintained implementation such as DESeq2 or edgeR, or limma-voom when its linear-model flexibility is needed. DESeq2 explicitly models count dispersion and stabilizes fold-change estimates ([Love et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4302049/), PMID 25516281). Voom estimates the mean-variance relation and supplies observation-level precision weights before limma empirical Bayes analysis ([Law et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4053721/), PMID 24485249).
- Before using a count-likelihood method on GSE255602, verify that the Salmon/tximeta gene-level pseudocounts retain a valid estimated-count scale and the required effective-length/offset information. If they are arbitrary normalized or noninteger summaries without compatible offsets, use an appropriate abundance-aware method or reprocess transcript quantifications/raw reads; do not force them into a count likelihood.
- For repeated tissues, clones, differentiations, or time points, use donor aggregation within a homogeneous stratum or a repeated-measures method such as voom plus `dream` ([Hoffman and Roussos, 2021](https://pubmed.ncbi.nlm.nih.gov/32730587/), PMID 32730587). Pure technical replicates may be collapsed after verifying their identity.
- Include measured nuisance variables in the design when they are estimable. Do not run ComBat across unrelated studies and then treat the pooled adjusted matrix as one experiment. Analyze studies separately and combine evidence afterward.
- Independently filter low-expression genes before testing. Preserve the unshrunk coefficient and standard error for meta-analysis; use a shrunken effect or moderated statistic for stable ranking.
- Control gene-level false discovery with Benjamini-Hochberg and report effect size, uncertainty, and FDR. A small p-value without a stable direction and meaningful effect is not a CMap query gene.

Run at least one prespecified alternate DE method as a sensitivity analysis. The primary method must be chosen before any drug score is seen; do not select whichever DE pipeline yields the more attractive compound list.

---

## 5. Cross-study replication and meta-analysis

1. Fit every study and biological stratum independently.
2. Lock at least one independent, nonoverlapping donor cohort measuring a comparable estimand as confirmation. Do not use it to select covariates, query size, or CMap aggregation rules. A distinct model system may provide orthogonal support but cannot be relabeled cohort replication.
3. Use a random-effects meta-analysis of log2 fold changes and standard errors only when at least three genuinely comparable independent cohorts are available; then report tau², I², Cochran Q, a prediction interval, and leave-one-cohort-out analysis. With only two cohorts, report the two effects and uncertainties separately and call same-sign evidence directional replication; tau², I², and a pooled common-effect claim would be too unstable here.
4. Where effect scales are incompatible, use a directional z-score or robust-rank aggregation only as sensitivity analysis. Fisher p-value combination alone discards direction. Robust Rank Aggregation supplies an explicit null model for ordered lists but does not solve biological incompatibility ([Kolde et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3278763/), PMID 22247279).
5. A replicated query gene must show the same disease direction in at least two independent nonoverlapping donor cohorts measuring a comparable estimand, as defined above. Cross-model concordance is separately labeled orthogonal support. Report discordant genes rather than forcing them into a consensus.

Do not call the intersection of significant-gene lists a meta-analysis, and do not batch-correct tissues into apparent agreement. Methods that combine within-study evidence are designed to retain study-specific biological and technical variation ([Rau et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4021464/)).

---

## 6. Query construction and CMap scoring

### Query construction

Freeze the exact gene annotation, symbol-to-Entrez mapping, CMap release, and gene universe. The primary query is:

- 100 replicated upregulated and 100 replicated downregulated BING genes, ranked by a locked combination of effect direction, magnitude, and uncertainty;
- sensitivity queries of 25, 50, and 150 genes per direction;
- BING-space WTCS as the primary published-CMap-compatible analysis, with a landmark-only query as a prespecified sensitivity analysis for dependence on inferred genes;
- failure, not padding, if too few stable genes survive.

Official CLUE guidance allows 10–200 genes per direction for external queries ([How to construct CMap queries](https://clue.io/connectopedia/how_to_construct_cmap_queries)). Keep total-tissue, composition-adjusted, iPSC, mechanistic, and adequately powered subtype queries separate.

### Scoring

The primary score should reproduce the published bidirectional CMap logic:

1. calculate a weighted connectivity score (WTCS) from separate up and down gene sets;
2. normalize within cell line and perturbagen type to obtain NCS;
3. express specificity as tau only against a frozen, matching reference-query compendium;
4. interpret negative connectivity as reversal, not as an efficacy probability.

The L1000 paper defines WTCS, NCS, tau, the landmark/inferred split, and the reference-query comparison ([Subramanian et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5990023/), PMID 29195078, DOI 10.1016/j.cell.2017.10.049). The [official CMap algorithm definition](https://clue.io/connectopedia/cmap_algorithms) specifies that tau compares each NCS with exemplar queries from Touchstone matched by cell and perturbagen type. A large negative tau is a connectivity heuristic, not an ALS-specific p-value.

Download and hash the GSE92742 `sig_metrics` file in addition to expression and identity metadata. Filter or stratify on documented signature quality and retain cell, dose, time, perturbagen identity, and replicate information. `sig_metrics` is necessary for quality filtering but is not sufficient to calculate tau: the run also needs the frozen Touchstone signature universe and matching exemplar-query reference distribution. If that reference cannot be reproduced, report WTCS, NCS, and a separately named empirical percentile—never relabel that percentile as tau. Level 5 is already a replicate aggregate; do not describe its 473,647 columns as independent experiments.

### Compound aggregation

Canonicalize salts, aliases, vendor names, BRD IDs, parent structures, and ChEMBL IDs before aggregation. The primary compound statistic must not be the single most-extreme reversing signature, regardless of whether the selected score convention encodes reversal as positive or negative.

Coverage is materially unequal in the committed Phase I metadata: sirolimus has 470 `trt_cp` signatures, riluzole 45, and edaravone 11. A best-of-many statistic therefore gives some compounds far more opportunities to produce an extreme outlier even under the null.

For every active moiety:

- retain each eligible cell × dose × time score;
- before viewing any reversal scores, freeze a score-blind dose/time eligibility rule using external PK or a class-wide rule and separate that primary subset from exploratory high-dose results;
- summarize within cell with a robust median after excluding failed-QC signatures;
- compare compounds on a prespecified common eligible-cell panel or use a cell-stratified hierarchical aggregate; require at least two high-quality cells with adequate overlapping coverage for a cross-context primary claim;
- report sign consistency, number of contexts, heterogeneity, dose/time dependence, and whether the result is high-dose-only;
- require support in at least two distinct high-quality cell lines, not merely two doses or times in one line, or label the signal context-specific;
- show the best single signature only as diagnostic detail;
- match empirical nulls on the number of available signatures so heavily tested compounds do not gain a multiple-opportunity advantage;
- report CMap's maximum-quantile cell summary as a prespecified sensitivity analysis, not as post-hoc rescue of a preferred compound.

Do not treat an unmeasured cell as neutral evidence. Matching the empirical null to each compound's signature count corrects best-of-many opportunity bias, but it does not make a median over two favorable cells estimand-equivalent to a median over nine different cells. Compounds lacking the locked overlapping panel must be labeled context-specific and not interleaved in the primary cross-compound ranking.

Most L1000 cells are not motor neurons. A compound whose target/pathway is absent in the assayed cell is not a convincing negative; a compound active only in one cancer line is not a convincing positive. Cell-relevance weights, if used, must be defined from baseline target/pathway expression before observing reversal scores.

Generic stress, cytostasis, and cell death can masquerade as reversal. A systematic oncology validation found that apparent reversal utility disappeared after accounting for general antiproliferative effects. That result does not directly validate or invalidate ALS reversal, but it establishes a failure mode that this pilot must counter-screen ([Koudijs et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9851289/), PMID 36445193).

---

## 7. Controls that can actually falsify the method

Use three distinct control layers. Passing one does not substitute for the others.

### A. Implementation and assay controls

- Reproduce at least one documented CLUE query with frozen inputs, expected WTCS/NCS/tau outputs, and a numeric tolerance.
- Pass a synthetic sign test: an identical perturbation signature must be positive connectivity/not reversal, an exact inverse must be the strongest negative connectivity/reversal, and flipping the query must flip the score sign.
- In leave-one-context-out tests, retrieve the same perturbagen or a chemically distinct member of its perturbagen class from high-quality signatures.
- Show vehicle and negative-control perturbations are null.
- Verify that reversing the query sign reverses score direction.
- Confirm compound aggregation is not correlated with the number of signatures available per compound under matched null data.

Predeclare quantitative thresholds, for example same-perturbagen recovery in the top 1% and perturbagen-class enrichment at FDR < 0.05 across held-out contexts. These proposed cutoffs are pilot heuristics, not published universal validity thresholds. They validate implementation and internal assay behavior, not ALS efficacy.

### B. ALS-relevant biological calibration

- Use the TDP-43 × ATXN2-knockout factorial contrast in GSE261875 as a mechanistic protective-intervention calibration, with the limitations of a single isogenic system stated.
- Consider treated patient-derived ALS motor-neuron expression from [GSE209696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE209696) as a ropinirole rescue/calibration dataset. Ropinirole came from an ALS iPSC phenotypic screen ([Fujimori et al., 2018](https://pubmed.ncbi.nlm.nih.gov/30127392/), PMID 30127392), but its small Phase 1/2a trial did not show an ALSFRS-R difference in the double-blind period and requires further validation ([Morimoto et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37267913/), PMID 37267913).
- Riluzole and edaravone may be reported as descriptive benchmarks if adequately represented. Failure to rank them does not by itself invalidate CMap because modest clinical benefit need not arise through a strong transcriptional reversal in a non-neural cell.
- Tofersen is a genotype-specific antisense oligonucleotide, not a small-molecule CMap positive control.

### C. Negative controls and empirical nulls

The reviewed “positive-control” list contains compounds that must not be treated as known disease-modifying successes:

- **Minocycline** caused faster ALSFRS-R decline than placebo in a 412-participant Phase III trial ([Gordon et al., 2007](https://pubmed.ncbi.nlm.nih.gov/17980667/), PMID 17980667).
- **Memantine** showed no primary- or secondary-outcome efficacy in a randomized ALS trial ([de Carvalho et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20565333/), PMID 20565333).
- **Gabapentin** did not establish a significant disease-modifying effect in the cited placebo-controlled ALS study ([Miller et al., 1996](https://pubmed.ncbi.nlm.nih.gov/8960715/), PMID 8960715).
- **Baclofen** is used for symptom management; that does not make it a disease-reversal positive.

The current code's top-10%/top-60% labels are arbitrary and uncalibrated, and they inherit the reversed rank direction. The list also counts rapamycin and sirolimus separately before aliasing them to the same active moiety. Canonicalize identity before evaluating any control panel.

If present in the frozen LINCS release, failed or harmful ALS interventions may form a descriptive clinical non-success panel. They are imperfect negatives because exposure, timing, population, and trial design can fail independently of mechanism.

The primary nulls are computational:

1. permute disease labels at the **donor** level within valid exchangeability blocks, carry all regions from each donor together, preserve valid region/batch availability, and rerun the full DE → query → CMap → compound-aggregation pipeline;
2. use reduced-model residual permutations when nuisance covariates make direct label exchange invalid;
3. generate random up/down gene sets matched on size, expression, variance, BING/landmark status, and gene-correlation structure;
4. run at least 1,000 end-to-end permutations for initial calibration, recognizing that this gives only about 0.001 minimum p-value resolution and is generally inadequate for BH testing across thousands of compounds; extend surviving tests sequentially toward 10,000 or more permutations or use a prespecified max-statistic shortlist gate;
5. compute each compound's null using the same available cells/doses/times as its observed statistic, estimate p-values as `(1 + null scores at least as extreme) / (B + 1)`, then control FDR at the compound level.

Report the null top-`k` hit rate in addition to BH q-values. If the observed shortlist is indistinguishable from label-permuted or matched-gene queries, the screen is killed regardless of how plausible individual drug names sound.

---

## 8. Compound identity, exposure, and safety gates

ChEMBL and PubChem are identity/annotation sources, not efficacy validators.

For every surviving compound, record:

- canonical active moiety, parent structure, salt/formulation, InChIKey, PubChem CID, and ChEMBL ID;
- current approval, marketing, withdrawal, and indication status from a regulator or current label, not `max_phase` alone;
- concentration and duration producing the L1000 signal versus measured or credibly modeled unbound human exposure;
- plasma, CSF, and brain exposure where the proposed mechanism requires CNS access;
- protein binding, metabolites, route, dosing burden, boxed warnings, hepatic/QT/respiratory/sedation risks, and major ALS-care interactions;
- whether a peripheral mechanism makes BBB penetration unnecessary;
- target/pathway expression in human motor neurons and relevant glia.

`MW < 450` and `HBD <= 3` may be retained as coarse cheminformatics descriptors but must not be labeled “CNS penetration.” A high-dose-only L1000 result is a failure flag only when the in-vitro and human concentrations have a defensible unbound, compartment-aware comparison; otherwise exposure remains unresolved. `max_phase = 4` means a drug reached approval for some indication; it does not establish current availability, safety in ALS, or a usable formulation.

Screen for generic cytotoxicity, proliferation arrest, and assay stress. PRISM and DepMap contain cancer-cell data; they can help identify cytotoxic or target-dependency patterns but cannot demonstrate tolerability in non-transformed human motor neurons.

---

## 9. Required uncertainty and robustness outputs

Every candidate table must include, not hide:

- donor counts and study/model provenance for each contributing disease signature;
- effect direction, effect size, standard error/FDR, and cross-study heterogeneity for query genes;
- landmark versus inferred-gene contribution;
- all eligible cell/dose/time scores, the locked aggregate, and the best-score diagnostic;
- compound context count and a context-count-matched empirical null;
- donor-bootstrap selection probability and rank interval;
- leave-one-donor/cohort/context-out sensitivity;
- query-size, DE-method, composition-adjustment, and aggregation sensitivity;
- label-permutation p-value and compound-level q-value;
- identity, current regulatory status, exposure comparison, and safety flags;
- explicit evidence tier: computational only, orthogonal transcriptomic replication, human ALS iPSC rescue, or external wet-lab validation.

Use deterministic seeds where applicable and capture the operating system, package lockfile, command, config, checksums, timestamps, stdout/stderr, and code commit. Raw downloads belong outside Git; a manifest with URLs and cryptographic hashes belongs inside Git.

Repository preflight found 181 raw experimental-data blobs tracked on current `origin/main` (57,972,488 bytes, about 55.29 MiB), including both `GSE124439_RAW.tar` and extracted count files. With six inventory files and one PRO-ACT application document, `data/` contains 188 tracked blobs totaling 58,016,260 bytes. The raw subset contradicts `AGENTS.md` and `.gitignore`; adding ignore rules does not untrack existing blobs. This does not by itself invalidate an analysis, but it breaks the repository's declared data-governance contract and obscures which files are authoritative inputs versus duplicated downloads.

The newly merged `pipelines/data_qc/te_contamination_check.py` usefully reproduces the aggregate finding that colon-delimited repeat rows carry 25.51% of GSE124439 raw counts, and its synthetic self-test passes. It is not yet a reliable mandatory gate. The committed exact-demo whitelist is missing; its per-sample fractions divide each sample's class count by the **grand total across all 176 samples** instead of that sample's library size, producing reported TE fractions of 0.06%–0.35%. Direct per-library calculation gives 12.78%–42.05% (median 26.40%). Without a whitelist, valid symbol-shaped rows are labeled `gene_like_unverified` and then counted as non-gene; an all-gene A1BG/TP53 toy matrix therefore returns 100% non-gene and `NO-GO`. The current self-test supplies a whitelist and misses both defects. The exact demo whitelist is not committed/versioned, and the plain-text loader cannot directly consume the README-linked gzip-compressed NCBI file. The README also says the executed exp001 signature failed because TE rows were normalized in, while the handoff implementation and merged review explicitly remove them before normalization. Fix and test both verdict paths and the denominator, commit/freeze the whitelist, support or stop advertising compressed input, justify thresholds, obtain the required independent review, and distinguish the broken main implementation from the handoff run before making this gate binding.

---

## 10. Go/no-go and falsification rules

### Locked pilot thresholds

These are deliberately conservative pilot heuristics, not universal CMap or clinical-validity standards. Freeze them, the score convention, the eligible cell panel, and all seeds in a hashed config before any drug name or rank is viewed:

1. **Disease-query gene:** discovery BH q ≤ 0.05, absolute log2 fold change ≥ 0.5, and the same coefficient sign in at least 80% of 1,000 donor bootstraps. Confirmation requires the same direction and one-sided BH q ≤ 0.10 in an independent nonoverlapping donor cohort measuring the same estimand. If no such cohort exists, the query remains exploratory and cannot pass the ranked-table gate.
2. **Assay recovery:** the frozen same-perturbagen control must rank in the top 1% in its held-out context, and held-out perturbagen-class enrichment must have BH q < 0.05. Vehicle/negative controls must show no systematic score shift, and their fraction outside the central 95% matched-null interval must be compatible with 5% under a prespecified binomial or permutation-calibrated bound. Numeric reproduction tolerance for each frozen CLUE fixture must be locked from the reference precision before implementation testing.
3. **Compound significance and context:** compound-level empirical BH q ≤ 0.05; reversal direction in at least two eligible high-quality cells and at least two-thirds of the locked eligible cells; no promotion when one cell alone changes the aggregate sign. Cross-context candidates must use the same eligible-cell panel or a locked cell-stratified model.
4. **Stability:** top-1% selection in at least 70% of 1,000 donor bootstraps, a 95% bootstrap rank interval contained within the top 10%, and unchanged aggregate direction in every leave-one-cohort-out analysis. Report failure rather than relax these thresholds after ranks are known.
5. **End-to-end null:** every promoted compound must retain compound-specific empirical BH q ≤ 0.05 under a signature-count- and context-matched donor-label null, and the observed maximum standardized reversal must exceed the 95th percentile of the permutation maxima (Westfall–Young/max-statistic gate). If discovery-count calibration is also reported, calculate each permuted dataset's q-values leave-one-permutation-out; do not score a permutation against an ensemble containing itself. No recognizable drug name overrides these gates.
6. **Exposure label:** when both sides are defensibly comparable, define `exposure margin = measured unbound human steady-state concentration in the relevant compartment / measured or estimated unbound in-vitro concentration at the lowest replicated, direction-consistent active dose`. Require the activity direction at two adjacent doses or in an independent signature so one outlier cannot define the denominator. Margin ≥ 1 is exposure-plausible, 0.3 to <1 is unresolved, and <0.3 is high-dose-only/implausible for this pilot. Use a lower confidence bound when available. If in-vitro binding, relevant-compartment human PK, or another needed conversion is missing, label exposure preliminary/unresolved and do not infer BBB access from descriptors or compare nominal media concentration directly with unbound human exposure.

### The pipeline may produce a ranked hypothesis table only if

- donor identity and repeated measures are resolved;
- no primary contrast is perfectly confounded;
- the disease signature meets the locked gene thresholds in at least two independent nonoverlapping donor cohorts measuring a comparable estimand;
- implementation controls pass their locked thresholds;
- compound aggregation meets the locked common-cell, context-consistency, and signature-count-imbalance rules;
- results outperform donor-label and matched-gene nulls after compound-level multiple-testing control;
- top results meet the locked leave-one-cohort and donor-bootstrap thresholds;
- score-blind dose eligibility was frozen before ranking, and detailed human exposure plus current safety/regulatory status are documented before promotion.

### Kill or downgrade the screen if

- the GEO parser fails or the sample/donor/region manifest does not reconcile to accession metadata;
- the LINCS checksum, decompression, GCTX schema/orientation, or row/column identity assertions fail;
- the synthetic score-sign test or frozen CLUE reproduction test fails;
- sample columns are used as independent `n` when donors repeat;
- no reproducible disease signature exists;
- the implementation cannot recover assay controls;
- rank is driven by one cohort, one cell, one dose, one outlier signature, or the number of signatures tested;
- generic stress/cell-death programs dominate;
- the ranking is null-equivalent under permutation;
- relevant human ALS motor-neuron data show the opposite direction;
- a defensible unbound, compartment-aware comparison shows required exposure is implausible, exposure cannot be resolved before promotion, or safety is incompatible with the proposed use.

The strongest surviving label before functional validation is **computational repurposing hypothesis**. “ALS drug,” “therapeutic candidate,” “validated hit,” and clinical advice are prohibited labels for this stage.

---

## 11. Minimum next implementation sequence

1. Repair the GEO parser, freeze a donor-aware metadata manifest, and reconcile every GSE124439 group/region/donor identity.
2. Decompress and checksum LINCS, parse the official GCTX schema/orientation, and assert row/column identities before scoring.
3. Replace the custom two-group model with a maintained donor-aware DE method and locked covariates.
4. Analyze GSE124439 regions and GSE255602 donor-level strata separately. Treat GSE255602 as limited cross-model orthogonal support, not a confirmatory postmortem cohort, unless a different accession-level dataset supplies an independent, nonoverlapping donor cohort measuring a comparable estimand.
5. Produce replicated up/down gene sets and one-to-one landmark/BING mappings before touching compound ranks.
6. Download `sig_metrics` and the required frozen reference compendium; implement quality filtering and WTCS/NCS, and use tau only if its reference distribution is reproducible.
7. Pass synthetic sign and frozen CLUE-query tests.
8. Replace best-signature compound ranking with a robust context-aware aggregate and context-preserving null.
9. Replace the current positive-control list with the three control layers above.
10. Add donor-label permutations, matched-gene nulls, bootstrap stability, and compound-level FDR.
11. Apply the already frozen score-blind dose/time eligibility rule during aggregation; after statistical survival, perform detailed identity, current-status, measured-exposure, and safety review without changing the cutoff to rescue a rank.
12. Run an independent adversarial review. Keep `VERDICT.md` pending until that review examines actual outputs.

This order prevents a memorable drug name from influencing the statistical choices that generated it.
