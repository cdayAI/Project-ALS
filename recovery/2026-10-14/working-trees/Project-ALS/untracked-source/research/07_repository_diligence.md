# Project ALS — Repository Diligence and Evidence Ledger

**Audit status:** independent source/pipeline and remote-state read-only diligence, with local audit/task-board edits only; current through 2026-08-25 22:09 EDT.

**Repository:** [cdayAI/Project-ALS](https://github.com/cdayAI/Project-ALS)

**Default-branch snapshot:** `origin/main@e734295e8a9f9e160eab08ca36bb9b3a2576cb05`

**Material branch snapshots:**

- exp001 handoff: `origin/exp001-sprint1-handoff@fcfff6470f76f1ef38bed3f82c6c69a35ee77417`
- causal-target source: `origin/pipeline/causal-targets@9c291488c62c46ca15e7e715a9dc0d3e231a53cc` (content squash-merged to main as `d25deea`)
- C9/AnswerALS source: `origin/exp/h007-c9orf72-module@591971429ab91b7c1594a3e922a3b7973c0ccc7` (merged to main)
- structure-readiness source branch: `origin/pipeline/structure-readiness@1eaf9e0928aff3c9289ab5b4eb8aaf721a7895ed` (merged to main)

**Decision boundary:** this report distinguishes committed facts, independently reviewed conclusions, active-branch WIP, reproduced behavior, and proposed work. It does not treat a heartbeat, green syntax check, simulated result, model-confidence score, or unreviewed branch as scientific validation.

---

## 1. Executive verdict

**Overall disposition: HOLD for scientific or translational decisions.** Project ALS is a promising, unusually candid early research factory, but it is not yet a reproducible discovery platform. No drug, target, biomarker strategy, docking pose, or therapeutic hypothesis in the reviewed repository is validated for advancement.

The strongest outcomes are negative or bounded:

1. exp001 is correctly **KILLED** and yields no valid compound. Its reviewer found serious statistical failures, while this audit found an even more basic score-direction inversion the reviewer missed.
2. H-007 is defensibly **KILLED only as the preregistered module computed in one small discovery dataset**. It is not evidence that C9orf72 neurons lack nucleolar-stress biology.
3. H-007b is now **KILLED with a merged independent review**. The reviewer reports a fresh reanalysis and the conservative negative result is credible, although the review code/output bundle is absent and the merged front matter/config/verdict artifacts remain malformed and contradictory.
4. H-008 is useful simulation plumbing labeled `SIMULATED_ONLY`; its headline 25% improvement versus the specified clinical comparator is not met.
5. The trial matcher is a useful reviewed retrieval prototype with appropriate medical boundaries, but is not reproducible or an eligibility engine.
6. Causal-target and structure work are a reviewed association/tractability inventory and an unreviewed environment bring-up, respectively—not causal target validation or prospective docking proof.

The repository's scientific governance prose is stronger than its enforcement. There is no collected automated test suite, CI workflow, dependency lock, protected-branch rule, or complete canonical run package. The data-QC CLI has a synthetic self-test, but it misses both reproduced logic defects. Raw experimental data is tracked despite an explicit prohibition. These gaps make even honest negative results harder to reproduce than they should be.

### Evidence-status ledger

| Workstream | Current evidence status | Defensible claim | Claim that is not supported |
|---|---|---|---|
| exp001 transcriptome reversal | **KILLED; no candidate valid** | The run and current ranking method failed | Any named compound reverses ALS or deserves functional follow-up from this run |
| H-007 / exp002 | **Narrow one-dataset operational kill** | The specified module, as computed in GSE303931, missed its preregistered bar | C9orf72 neurons lack nucleolar/speckle biology; two independent datasets confirmed the biology |
| H-007b | **Review-confirmed kill; reported reanalysis** | Count-level cytoskeleton enrichment missed the 2.0 bar and ropinirole rescue failed in both disease genotypes | A replacement cytoskeletal module is rescued or ready for LINCS screening |
| H-007c/H-011/H-012 | **Hypothesis stubs; data access/gates unresolved** | Post-hoc follow-up conditions and power gates are written down | A new C9 module, repeat gradient, or pharmacodynamic readout has been detected |
| H-008 NfL enrichment | **Simulation plumbing only** | Code produces deterministic synthetic power tables and exposes design questions | NfL dynamics reduces ALS trial size by 25%, predicts outcomes, or is validated in PRO-ACT |
| H-009 subtype bridge | **Not implemented / non-identifiable as drafted** | The hypothesis is an idea ledger entry | Answer ALS molecular labels can be directly regressed onto separate PRO-ACT patients |
| Trial matcher | **Reviewed bounded prototype** | It retrieves and ranks public trial records and shows source text for human review | It determines eligibility, exhaustively retrieves under every option, or reproduces a historical run |
| Causal targets | **Reviewed main-branch association inventory** | It filters known placeholder labels and joins 79 reported GWAS gene labels with tiered tractability lookups | It identifies causal ALS genes, colocalizes signals, or supplies target direction of effect |
| Structure readiness | **Unreviewed environment smoke test** | Boltz was reportedly installed and a known complex was reportedly reconstructed | A prospective pose, affinity, ALS target, or therapeutic mechanism has been validated |

---

## 2. Stop findings

### P0 — No candidate or target may advance from current outputs

The repository contains no decision-grade positive result. The correct action is to preserve the three negative dispositions—exp001, H-007, and H-007b—repair their provenance, and prevent candidate language from leaking out of status notes, branch snapshots, or tables.

### P1 — exp001 ranking selects disease mimics as “reversers”

Both main and the handoff define `score_reversal = -cosine`. Therefore:

- exact disease-direction similarity: cosine `+1` → reversal score `-1`;
- exact inverse signature: cosine `-1` → reversal score `+1`.

The code then sorts ascending and aggregates each drug's minimum score. That selects the most disease-similar signature, not the strongest inverse. The handoff README correctly states that larger values mean stronger reversal, making the contradiction explicit. Its top 30 annotated rows are all negative (`-0.388` through `-0.348`), so the advertised top list comprises the strongest disease-direction similarities under the documented formula.

The merged review correctly kills the run for exposure/min-of-N bias, confounding, failed replication, and annotation failure, but incorrectly states that the scorer is correct. Its percentiles and exposure-matched controls therefore analyze the wrong direction. The kill remains safe, but the review cannot be used to validate the scoring implementation or design exp001b until direction is corrected and tested.

Additional exp001 blockers:

- Main's parser reproducibly raises `AttributeError: 'list' object has no attribute 'startswith'` at `run.py:51` on the committed matrix.
- Main opens a compressed `.gctx.gz` through `h5py`, uses wrong metadata paths/orientation, and has an undefined `h5py` reference found by Ruff.
- The handoff fixes transport/schema/alignment but retains the sign/order error.
- The disease model treats 145 ALS and 17 control libraries as independent, although they represent 76 ALS and 12 control subjects and 48 subjects span regions. The model has no donor blocking or region/covariate adjustment.
- Main has a final verdict and review but no exp001 outputs. The handoff omits the full DE table, all signature scores, primary ranking, sample sheet, and annotation cache.
- Handoff code and its committed annotation-table schema disagree, the large-file checksum is truncated, and dependencies are unpinned.
- The review's GSE255602 replication code and outputs are not committed.

**Required disposition:** do not cite any drug name from exp001. A single retry is reasonable only after a synthetic sign test, donor-aware query construction, context/exposure-matched empirical nulls, valid precommitted controls, and a complete canonical run package.

### P1 — default-branch verdicts are not reproducible from default-branch code

A clean checkout cannot reconstruct the central results:

- exp001's executing implementation and curated snapshot exist only on a divergent handoff branch; main retains the older broken code and zero outputs;
- exp002 hard-codes analysis choices rather than reading `config.yaml`;
- the original GSE283507 result has no producing code;
- count-level GSE283507 remediation code and outputs are now on main, but clean reconstruction still live-downloads unpinned count/GO inputs; the original TPM arm remains non-evidentiary, and GSE303931/GCTX plus other runtime inputs are not committed;
- mutable GO, BioMart, API, model, and weight inputs lack frozen versions/hashes;
- structure outputs and logs are explicitly ignored.

A verdict file is not a reproducibility bundle. Each result needs code, exact config actually consumed, input manifest/checksums, environment lock, output hashes, run log, and review tied to the same immutable commit.

### P1 — no enforceable software-quality gate

Measured at `e734295` (the two commits after `d25deea` modify only `digests/status.log`):

| Check | Result | Interpretation |
|---|---|---|
| Tracked files | 269 | 188 are under `data/`; code surface remains small |
| Python sources | 11 | All compile syntactically |
| Pytest tests | 0 | `pytest -q` collected none and exited 5; the separate data-QC CLI self-test misses both reproduced defects |
| CI workflows | 0 | No `.github/workflows/` |
| Dependency manifests/locks | 0 | No `pyproject.toml`, requirements file, environment file, or lock |
| Ruff | Failed: 43 findings | Includes real `F821 Undefined name h5py` in exp001 plus new gate-check hygiene findings |
| Branch protection | Absent | GitHub API reports main unprotected |
| Repository rulesets | 0 | Governance rules are advisory only |
| Merged pull requests | 2 | PRs #2 and #3 had no GitHub reviews and no checks; PR #3 did include an in-repo review artifact |
| First-parent main commits after `AGENTS.md` | 28 | Most governance happened outside an enforced PR/check path |

Positive checks should be preserved but not overstated: all 11 Python files compile; all five JSON and 13 CSV artifacts parse; 49 tracked Markdown files have no broken relative links; `git fsck` passed; and all 179 tracked gzip streams plus the tar archive passed integrity checks. One of five `.yaml` files and one of 14 hypothesis-record front-matter blocks fail YAML parsing.

### P1 — raw-data policy is already violated

`git ls-files -ci --exclude-standard` reports 181 ignored-but-tracked raw files occupying 57,972,488 Git-blob bytes (55.29 MiB): 178 GSE124439 files and three LINCS metadata files. `GSE124439_RAW.tar` duplicates much of the extracted material. The raw files and ignore rule were added together in commit `7869f60`; ignoring a tracked file does not remove it from history.

This contradicts `AGENTS.md`, `.gitignore`, inventory prose, and heartbeats claiming no Git-hygiene violation. It also creates unaddressed third-party redistribution/licensing questions because the MIT code license does not automatically license deposited datasets.

**Required disposition:** decide what the project is legally allowed to redistribute; preserve accession, checksums, download scripts, and immutable manifests; remove raw blobs from the current tree in a dedicated reviewed change; and make any history rewrite a separate, explicitly authorized operation.

---

## 3. Scientific workstream findings

### 3.1 H-007: the kill is narrower than the status language

The merged review correctly discards the original GSE283507 TPM arm as non-evidentiary: its tested universe shrank from about 16,253 to about 4,484 genes, GO sets collapsed, within-group variance was pathological, and no producing code was committed for that arm. Later count-level GSE283507 remediation belongs to H-007b and does not rehabilitate the original H-007 test. The only defensible H-007 conclusion is that the preregistered nucleolar/speckle module, as computed in the 12-library/two-isogenic-pair GSE303931 experiment, failed its bar. The [GSE303931 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303931) confirms the very small design.

That does not establish biological absence. Two line pairs cannot robustly separate genotype from differentiation, culture composition, or pair-specific effects. The recovered Y-linked genes are an unexpected confound in an isogenic contrast, not a known-truth positive control. The committed code calculates one pooled model and does not support the repeated claim that the Y signal was shown in both pair-specific contrasts.

There is also an enrichment-tail bug: when observed overlap is zero, code tests random overlap `>=1`, even though the correct upper-tail event is `>=0` and therefore has probability 1. This can make zero-overlap rows look spuriously unusual.

**Allowed wording:** “The specified module, as estimated in GSE303931, missed the preregistered operational bar.”

**Disallowed wording:** “Two datasets prove C9orf72 neurons lack this biology.”

### 3.2 H-007b: kill confirmed, record still inconsistent

The remediation and independent review are now merged. The reviewer reports independently downloading the raw data and reanalyzing the tested-gene/module counts, obtaining a primary cytoskeleton fold of `1.535` versus the committed `1.578`; both are below the locked `2.0` gate. It also reports confirming negative ropinirole rescue indices in both disease genotypes. Reviewer code and a review output bundle are not committed, so this is not a fully reconstructable independent reproduction. The two failed committed gates and conservative disposition nevertheless make the H-007b kill credible.

The merged record still contains defects:

- count-level primary cytoskeleton fold `1.578`, below the revised `2.0` bar;
- rescue indices `-0.014` and `-0.024` in the committed CSV, both failing;
- an experiment verdict that separately states `-0.112` and `-0.101`, which does not match that CSV;
- H-007b front matter is invalid YAML because its `provenance` quote is unterminated;
- `experiments/exp003_h007b_reversal/config.yaml` also fails YAML parsing and preserves the obsolete disjunctive `1.5`-fold criterion rather than the revised gate;
- exp003's `VERDICT.md` begins `PENDING`, calls FC-A passed, leaves FC-C/FC-D pending, and reports stale `-0.112/-0.101` rescue values before appending a final `KILLED-SUPERSEDED` closure; the current CSV/update instead reports `-0.014/-0.024`;
- exp003's advertised reusable LINCS path raises `ValueError` because it calls `np.save` on a dictionary with `allow_pickle=False`; even past that line, it returns score-only data while the next stage requires `pert_iname` and `pert_type` metadata, so aggregation would raise `KeyError`;
- the review's reanalysis code, command log, environment, and output bundle are not committed;
- the review flags the two-minute criterion-commit-to-results interval as implausibly short for the claimed pre-analysis chronology;
- the permutation universe includes deposited genes outside the DE-detectable set, so future competitive nulls should be detection matched.

These defects do not overturn the conservative, review-confirmed kill. They do prevent exp003 from serving as a clean reusable experiment bundle. Its LINCS arm should remain parked, and the appended `KILLED-SUPERSEDED` closure should become the sole canonical status by removing or reconciling the stale header, stages, config, and rescue values.

### 3.3 Data-QC checker: useful aggregate result, wrong per-sample metric

The newly merged checker correctly reproduces that colon-delimited repeat rows account for 25.51% of aggregate GSE124439 raw counts. Its synthetic self-test passes and a direct independent calculation confirms the aggregate number.

The per-sample implementation is wrong. `class_metrics` divides every sample's class count by the grand total across all 176 libraries, not by that sample's total library count. The committed report therefore shows TE fractions of 0.06%–0.35%. Direct per-library calculation gives 12.78%–42.05%, median 26.40%. The self-test checks only class counts, aggregate fractions, and verdicts, so it cannot detect the denominator error.

No-whitelist verdicts are also logically broken. Heuristic classification labels symbol-shaped rows `gene_like_unverified`, while `make_verdict` treats every class other than literal `gene` as non-gene. A clean two-gene toy matrix therefore returns 100% non-gene and `NO-GO`. The exact demo whitelist is not committed or version-hashed, and `load_whitelist` uses plain `open`, so it cannot directly consume the gzip-compressed NCBI file linked by the README. Thresholds are asserted rather than calibrated, and README incorrectly says the executed exp001 signature failed because it normalized with TE rows. The handoff code and its review explicitly removed TE rows before normalization; the broken main code never reaches a valid run. The checker also merged without the independent review required by the repository's own hard rule. Keep the aggregate contamination finding, fix and test both verdict paths and the per-sample metric, commit/hash the reference, and correct the provenance claim before making this a mandatory gate.

### 3.4 H-008: reproducible simulation does not validate the hypothesis

An isolated rerun on Python 3.13.7, NumPy 2.4.3, pandas 3.0.1, and SciPy 1.17.1 reproduced `power_curves.csv`, `sample_sizes.csv`, and `simulation_summary.md` byte-for-byte. `run_meta.json` differed only in the last floating-point digit of the F-test p-value. That is good evidence that the synthetic plumbing is nearly deterministic, and evidence that the environment must be pinned.

The scientific headline still fails in simulation: at effect `0.30`, `nfl_ml` reduces estimated `n` from 34.4 to 28.7 versus `clinical_static`, only 16.6%, not 25%. The 28.6% number uses the weaker random unenriched comparator.

Further design defects:

- half the subjects receive in-sample predictions and then compete with held-out subjects for enrollment;
- the written noise control is ±5%, but code accepts reductions down to -10%;
- written oracle dominance over learned arms is not enforced;
- F2 promises a per-patient-benefit check, but code does not report it and fixes mean treatment benefit by construction;
- the “none” arm still keeps a random 50% and reports two screened per enrolled;
- the clinical slope variable is negative for decline despite a comment saying positive;
- the comparator drifts between an early-slope model and the fuller ENCALS model named in the access draft;
- the very small NfL-slope p-value is generated by construction.

No PRO-ACT or linked longitudinal NfL records are used. No claim about real ALS patients, biomarker validity, prospective power, or trial enrichment is supported.

### 3.5 Trial matcher: good boundaries, incomplete reproducibility

This is the strongest governed implementation. It has merged code and a separate review, uses ClinicalTrials.gov's public API, checks structured locations, distinguishes fetch failure from zero results, quotes source text, and explicitly refuses to decide eligibility.

Live QA with `--country "United States" --mutation sod1 --limit 5 --max-pages 1 --json` fetched 71 records and exited successfully. It also showed why the source-quote boundary matters: a tofersen/SOD1 study received only a low score because the keyword appears in an exclusion-text clause with different scope.

Remaining defects:

- JSON ignores `--limit 5` and emitted all 71 results (204,441 bytes);
- `max_pages` can silently return a partial result while a continuation token exists;
- country confirmation is substring-based;
- text excerpts are truncated to 220 characters despite being described as the actual source sentence;
- tied scores preserve unstable live API order;
- no frozen response fixture, response hash, retrieval timestamp, provider `dataTimestamp`, or regression suite exists.

Keep it as advisory retrieval. The trial site and treating clinicians remain the eligibility authority.

### 3.6 Causal-target merge: reviewed scaffold, not causal validation

Main now contains Stream A plus `reviews/causal_targets.md`. The review performed useful adversarial work: it caught and fixed the pre-merge corruption that resolved `NR` and `intergenic` placeholders to NRAS and GAPLINC, separated approved-drug flags from weaker [Open Targets tractability tiers](https://platform-docs.opentargets.org/target/tractability), and improved top-level fetch failure messaging. The current 79-row outputs omit those placeholders and are internally parseable.

The `SURVIVES` verdict applies only to this extraction scaffold. Of the 79 reported gene labels, 20 meet the audit's `p <= 5×10^-8` tier; 59 are merely within the Catalog's [literature-curated top-association inclusion range](https://www.ebi.ac.uk/gwas/docs/faq) (`p < 1×10^-5` in all stages), which is not a causal-gene designation. This is a different extraction from `data/inventory/gwas_catalog_als_loci.json`, which has 56 labels across 97 gene-association rows and 37 significant rows. They overlap on 55 real labels, the pipeline adds 24, and label-level versus row-level significance counts are not interchangeable; neither path freezes the raw response, so there is no canonical inventory. Author-reported gene presence is not causal-gene assignment. The table collapses evidence to minimum p-value and mention count, drops rsID/position/effect allele/beta or odds ratio/SE/MAF/ancestry/sample-size/LD context, and turns multiple labels attached to one locus into separate apparent targets. It has no fine-mapping, causal-gene probability, eQTL/sQTL colocalization, Mendelian randomization, perturbational direction, or target-safety logic.

Partial Open Targets failure remains misrepresented. Per-gene network exceptions are converted to empty/false tractability, and the constructed `error`, `resolved`, and returned approved-symbol fields are not serialized. HTTP-200 GraphQL `errors` are not checked, so API failure remains indistinguishable from genuine lack of tractability. The positive control is also only a four-name presence check: a synthetic table containing just `SOD1`, `C9ORF72`, `UNC13A`, and `TBK1`, with no identifiers, p-values, mapping, or tractability, returns `passed: true`. The review's two byte-identical live runs establish short-term determinism, not reproducibility across mutable providers; outputs still lack retrieval times, provider releases, raw responses/hashes, environment lock, and a run manifest.

`COLOCALIZATION_PLAN.md` describes common/rare gene-list convergence, not statistical colocalization of a shared causal variant. Its mandatory TARDBP/FUS burden control is not a valid universal invariant across masks, cohorts, and power. The [PRJNA573105 record](https://www.ncbi.nlm.nih.gov/bioproject/573105) is an umbrella spanning public transcriptomic/SRA subprojects and genome-sequencing projects; it is not itself an analysis-ready WGS cohort. NYGC individual-level WGS requires authorized access through current [dbGaP accession phs003067.v3.p3](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs003067.v3.p3), so the plan has not identified an open, no-DAC WGS replication cohort. The merge also does not implement H-003/H-004's stated QTL PP4, therapeutic-direction, and AnswerALS gates, despite saying those hypotheses are unblocked. No target is decision-ready.

### 3.7 New C9 hypotheses: gated ideas, not evidence

H-007c is explicitly a post-hoc candidate stub and appropriately forbids promotion until an untouched out-of-sample dataset and power/variance gates pass. It nevertheless encodes the wrong follow-up: the independent review authorized a post-hoc **synapse** (`GO:0045202`) family, while the stub's primary gate reverts to **cytoskeleton organization** (`GO:0007010`); the proposed organoid-microglia context is not a like-for-like neuronal synapse replication. The repository records an AnswerALS attempt stopping at email verification, and this audit did not complete Neuromine authentication/data-request approval. The [official access process](https://app-wp-prod-001.azurewebsites.net/data-information/) requires an online DUA, applicable PI/director co-signature, and an institutional business-official signature. The repository therefore has no documented approved access or frozen cohort metadata; exact genotype-specific counts and repeat-length fields remain unverified rather than known absent, and the hypotheses' `open` labels overstate executable status.

H-011 incorrectly calls GSE124439 patient iPSC neurons. The material in this repository is postmortem cortex, and the required C9 annotation has not been established. H-012 requires intron-level PSI/excision ratios, but the cited available files are gene-level counts; it needs junction-, transcript-, or BAM-level input before its splicing test is executable. These files are useful preregistration scaffolds only. No C9 repeat-gradient, splicing-load, or treatment-response result exists.

Main now appends `KILLED-SUPERSEDED` to exp003 and includes an AnswerALS gate checker, but the new scaffold was merged without adversarial review and is validated only on mocked metadata. It counts sample rows rather than unique `donor_id` values for the `n >= 15` rule, accepts the literal value `unknown` as a sex annotation, and tests repeat size only for non-nullness rather than `>26`. It does not implement H-012's preregistered pilot gate—PSI SD `<=0.15` for at least 30% of informative introns—and instead invents a `>=60` total-row proxy. Because all gates feed one verdict, that invented proxy can also fail an otherwise adequate H-011-only 15+15 cohort. This audit constructed 60 rows from only three unique donors, set every sex to `unknown`, and used C9 repeat size `1`; the checker returned `PASS` on all gates. The README also documents a nonexistent `--n-per-group` mock-generator option. It is a preliminary scaffold, not validated access or cohort readiness.

### 3.8 Structure readiness: known-complex bring-up only

Main reports a 57:52 CPU Boltz run on TBK1 kinase domain plus co-crystal ligand 1FV, using one recycle, 50 sampling steps, one diffusion sample, and no MSA. It reports confidence `0.843`, ipTM `0.963`, complex pLDDT `0.814`, and approximate ATP/hinge contacts.

The evidence artifact is missing. No predicted structure, confidence JSON, log, environment lock, package version, model revision, weight checksum, ligand RMSD, aligned 4IM0 comparison, PoseBusters/physics check, or independent review is committed. `outputs/` and `run.log` are ignored. Installation is bare `uv pip install boltz`; DiffDock was not installed.

The input is a known TBK1–1FV co-crystal, so success is compatible with benchmark familiarity or memorization and is not a prospective test. ipTM/pLDDT measure confidence, not pose accuracy or binding affinity. This supports only “the stack reportedly ran on one known complex.” It does not support “Stream D ready” for target or compound decisions.

---

## 4. Governance, security, and documentation

### 4.1 State is contradictory

At the audited main SHA:

- README says exp001 is running, although its verdict and review say killed.
- `AGENTS.md` still marks exp001 in progress.
- H-007 front matter remains `open` despite its confirmed kill.
- H-007b front matter says killed but is invalid YAML; exp003 starts with `PENDING` and stale passed/pending stages before an appended final closure.
- the hypothesis batch says it is not merged although it is on main;
- the same initial-batch entry still says all H-001–H-010 hypotheses are open even though H-007/H-007b have since been killed;
- Stream A is now marked complete and says it unblocks H-003/H-004, although their QTL PP4, direction, and AnswerALS gates remain unimplemented;
- structure is “DONE pending review,” but it was merged without a review and a heartbeat calls it ready;
- the power-analysis rule lives in `_LESSONS.md` and new drafts but not the numbered hard rules in `AGENTS.md`;
- `AGENTS.md` says H-007b confirmation commits remain unmerged even though they are now on main;
- data-QC is called verified/mandatory in heartbeats while its board entry says “to become mandatory,” no review exists, and its per-sample calculation is wrong;
- the latest heartbeat says there are no `SURVIVES` verdicts even though the immediately preceding causal-target merge includes `reviews/causal_targets.md` with `SURVIVES`;
- the same heartbeat calls repository hygiene clean despite the 181 ignored-but-tracked raw blobs, and claims exp001b has 17 outputs/in-progress scoring although no exp001b remote branch or default-branch output package exists;
- the two newest heartbeat lines self-date themselves `23:40 EDT` and `00:05 EDT`, but their Git commits are timestamped `22:07:38 EDT` and `22:08:05 EDT`, respectively—92 and 117 minutes earlier—so status-log chronology is not trustworthy;
- “Day-1 digest published” refers only to another status-log line; no digest artifact exists.

Status should be generated from one machine-readable ledger and checked in CI. A heartbeat may announce activity but must not override verdict, review, or branch state.

### 4.2 Medical README language is inaccurate

README says: “Four approved drugs extend life by roughly 2–6 months combined.” That combines unlike endpoints and is no longer current. At Amylyx's request, FDA withdrew RELYVRIO's approval effective August 29, 2025, after PHOENIX failed its prespecified primary and secondary endpoints ([published Federal Register notice](https://www.federalregister.gov/d/2025-16646)). FDA's review says analyses used for riluzole's approval suggested about a three-month survival benefit, while edaravone was approved on ALSFRS-R functional decline and has no demonstrated survival benefit ([FDA summary review](https://www.accessdata.fda.gov/drugsatfda_docs/nda/2022/216660Orig1s000SumR.pdf)). Tofersen is restricted to SOD1-ALS and received accelerated approval based on NfL reduction reasonably likely to predict benefit; FDA's current tracker lists the required ATLAS postmarketing trial as ongoing, with projected completion June 30, 2028 ([FDA Qalsody decision](https://www.fda.gov/drugs/news-events-human-drugs/fda-approves-treatment-amyotrophic-lateral-sclerosis-associated-mutation-sod1-gene); [accelerated-approval tracker](https://www.fda.gov/drugs/accelerated-approval-program/ongoing-non-malignant-hematological-neurological-and-other-disorder-indications-accelerated)).

A safer replacement is:

> ALS has no cure. FDA-approved treatment options include riluzole and edaravone for ALS and tofersen for adults with SOD1-ALS under accelerated approval. Their indications, endpoints, and evidence differ; none should be described by combining their effects into one survival number.

### 4.3 Security has a good base but controlled-data planning is incomplete

Positive findings:

- repository is public and MIT-licensed;
- GitHub secret scanning and push protection are enabled;
- GitHub reports zero secret alerts;
- a high-confidence credential-pattern scan found no current-tree or reachable-history credential hit;
- the trial matcher carries an appropriate no-medical-advice boundary.

Gaps:

- Dependabot security updates/alerts are disabled; no code-scanning workflow or `SECURITY.md` exists;
- exp002 uses BioMart over plain HTTP and mutable provider inputs are not frozen;
- no `CONTRIBUTING.md`, `CITATION.cff`, issue template, PR template, or dataset `NOTICE` exists;
- GitHub community-profile health is 28%;
- the PRO-ACT draft lacks explicit access control, encryption, retention/deletion, DUA/IRB handling, audit logging, and incident-response rules before patient-level data is accepted.

---

## 5. Required next sequence

### Immediate — preserve truth and stop claim drift

1. State explicitly in README and the ledger that exp001 has no valid candidate; narrow H-007 to its one-dataset operational conclusion; close H-007b/exp003 consistently as reviewed and killed.
2. Correct the ALS treatment paragraph with drug-specific, sourced language.
3. Make `AGENTS.md`, hypothesis front matter, verdicts, reviews, README, and status log agree. Add the power-analysis rule to the binding policy only after defining its required inputs and exceptions.
4. Keep H-007c/H-011/H-012 at stub/draft status until their data-access, power, variance, and replication preconditions are actually met. Replace the AnswerALS row-count proxy with hypothesis-specific donor and PSI gates and adversarial tests. Treat the merged causal-target output as reviewed scaffolding only until its error provenance, reproducibility bundle, locus-level evidence, and H-003/H-004 gates are implemented.

### Engineering gate — make every result reconstructable

5. Add a pinned environment, CI, Ruff, tests, artifact-schema validation, Markdown-link checks, and a status-ledger consistency check. Protect main and require a reviewing approval plus passing checks.
6. Define one run-bundle schema: code SHA, consumed config hash, input accessions/checksums, retrieval timestamps/provider versions, environment lock, seed, output hashes, logs, and reviewer SHA.
7. Fix the raw-data policy in a dedicated reviewed change. Keep downloadable manifests and checksums; do not rewrite history without separate explicit authorization.

### Scientific gate — test the method before another real rank

8. exp001b: add exact mimic/inverse sign tests, donor-aware DE, region/covariate handling, measured-versus-inferred L1000 distinctions, null-calibrated context aggregation, valid controls, and committed independent replication. Retire the stream if the preregistered retry fails.
9. H-007/H-007b: fix the zero-overlap permutation tail, use detectable-gene matched nulls, resolve rescue-value/config/verdict inconsistencies, and preserve the confirmed kill. Require adequate biological replication before another module claim.
10. H-008: cross-fit or use a separate training cohort, enforce written controls exactly, implement F2, define one comparator, represent unenriched enrollment correctly, and make no performance claim until linked real data and external validation exist.
11. Trial matcher: apply `--limit` consistently, fail or label partial pagination, use exact structured location matching, store governed response manifests/fixtures, and add regression cases from its review.
12. Causal targets: preserve raw associations and missing/error states; add significance tiers, ancestry/effect metadata, fine-mapping, real QTL colocalization, therapeutic direction, and orthogonal replication before ranking.
13. Data QC: fix the per-sample denominator and no-whitelist verdict, support or stop advertising compressed NCBI input, add regression tests with unequal library sizes, freeze/version the whitelist, calibrate thresholds, obtain independent review, and correct the exp001 provenance statement.
14. Structure: pin code/models/weights, commit compact outputs and hashes, add ligand RMSD and physical-validity checks, benchmark prospectively against held-out complexes and classical docking, and obtain independent review before calling the stream ready.

---

## 6. Completion and provenance statement

Completed in this audit:

- fetched and pinned the moving remote repeatedly through `e734295`;
- inspected main plus the material branch snapshots listed above;
- traced architecture, code, data, results, reviews, hypothesis state, Git history, and GitHub controls;
- reproduced the exp001 parser crash and sign logic;
- ran syntax, test-discovery, lint, artifact-parse, archive-integrity, repository-integrity, biomarker-simulation, and live trial-matcher checks;
- verified current FDA evidence needed to correct the README treatment claim;
- updated the local methods/evidence ledgers in `research/05` and `research/06`;
- created this local repository-diligence ledger in `research/07` and updated the local task-board entry.

Not done or not authorized:

- no scientific pipeline code was changed;
- no raw data was removed;
- no GitHub setting, issue, PR, review, or default branch was changed;
- no commit or push was made;
- no wet-lab, clinical, external-cohort, or partner validation was performed;
- `research/05`, `research/06`, `research/07`, and the task-board edit remain uncommitted local audit changes.

The correct negative result is valuable: **the repository currently has no validated edge or candidate, and its most urgent work is to make future negative results reproducible and future positive claims much harder to earn.**
