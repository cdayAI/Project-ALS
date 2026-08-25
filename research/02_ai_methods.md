# Where AI/ML Genuinely Moves the Needle in Drug Discovery for Neurodegeneration and ALS
**Research brief — evidence-checked against NCBI PubMed E-utilities (all PMIDs resolved programmatically, Feb 2025-era literature through 2025).**

Scope: 2024–2025 state of practice, with landmark primary literature. Every bracketed number is a verified PMID unless marked otherwise (preprints carry their Europe PMC/bioRxiv IDs). The tone is deliberately skeptical: each section ends with what the method has *actually* produced for neurodegeneration, and where it fails.

---

## 1. Target Identification from Omics

### Network propagation and network medicine
The core idea: map disease genes onto a protein–protein interaction network, then diffuse ("propagate") signal to find proximal, druggable neighbors. Menche et al. established the topological scaffold — disease modules occupy distinct network neighborhoods, and their separation predicts comorbidity [PMID 25700523]. Guney et al. operationalized this for drugs: compounds whose known protein targets sit close to a disease module are predicted efficacious, validated across many indications [PMID 26831545].

**Applied to neurodegeneration, the honest scorecard:** these methods are now standard *prioritization* layers in academic pipelines and biotech triage, but they have not by themselves produced a first-in-class approved neurodegeneration drug. Their real value is (a) rationalizing polygenic risk into biology (Alzheimer's disease GWAS now implicates ~80 loci spanning immunity, lipid metabolism, and endosomal pathways — Bellenguez et al. [PMID 35379992]; earlier Kunkle/Lambert waves [PMID 30820047, 24162737]), and (b) ranking repurposing candidates for wet-lab follow-up. Treat them as hypothesis generators, not oracles.

### Causal inference: genetics-first prioritization
The most defensible causal method is Mendelian randomization (MR) over the *druggable* genome — asking whether genetically perturbed expression of a target changes disease risk, which mimics the drug mechanism. This genuinely moves the needle because genetic evidence doubles clinical success rates when present at the chosen target (the classic pharma success-factor finding underpinning human-genetics-first strategies).

Concrete neurodegeneration outputs:
- Storm et al. screened the entire druggable genome against Parkinson's disease risk and surfaced prioritized, genetically supported targets now feeding industry programs [PMID 34930919].
- Zhu et al. applied transcriptome-wide MR to ALS and FTD, identifying putatively causal, druggable loci [PMID 38019415]; Duan et al. independently prioritized **TBK1** as a drug-repurposing entry point for ALS via druggable-genome MR [PMID 38443977].
- Integrative brain-proteomic + genetic prioritization extends the same logic to protein levels rather than transcripts [PMID 36759259].

Caveats: MR inherits pleiotropy and eQTL/pQTL confounding; results vary with tissue and cis-restriction choices. Multiple independent ALS MR papers disagree on specific hits — treat any single "causal target" claim as provisional.

### Multi-modal integration: genomics × proteomics × single-cell
Single-cell atlases converted bulk "disease vs control" signatures into cell-type-resolved states. Mathys et al.'s multiregion single-cell dissection of Alzheimer's cortex mapped progression-associated glial and neuronal states at scale [PMID 39048816], and an accompanying atlas identified correlates of cognitive resilience despite pathology [PMID 37774677]. A concrete mechanistic payoff: Blanchard et al. used single-cell analysis to show APOE4 impairs myelination via cholesterol dysregulation in oligodendrocytes — a cell-type-specific, druggable axis invisible in bulk tissue [PMID 36385529].

In ALS, unsupervised ML clustering of post-mortem motor cortex and blood expression identifies molecular subtypes beyond the sporadic/familial split [PMID 38129934] — the substrate for patient stratification (Section 4). Fluid proteomics plus ML yields diagnostically valuable CSF/plasma protein signatures in ALS [PMID 30397248].

A different modality worth naming: human **iPSC phenotypic screening with computational hit-calling** produced ropinirole as an ALS candidate from patient-neuron screening, advanced into a Phase 1/2a trial with suggested slowing in a subset [PMID 37267913]. This is among the very few omics/AI-era ALS target-to-clinic stories that reached patients.

**Bottom line:** multi-modal integration is where target ID is genuinely improving, mainly through (i) causal genetics filtering and (ii) cell-type resolution. Pure network-propagation papers remain largely non-translated.

---

## 2. Drug Repurposing Screens (LINCS L1000/CMap, PRISM, DepMap)

### How signature matching works
- **LINCS L1000/CMap**: ~1M+ gene-expression profiles across cell lines treated with thousands of compounds at multiple doses/times [PMID 29195078]. The query logic ("signature reversal"): compute a disease signature, rank all perturbagen signatures by similarity/dissimilarity (e.g., weighted connectivity score); compounds whose profiles anti-correlate with the disease state are candidate therapeutics.
- **PRISM**: pooled barcoded viability profiling of non-oncology drugs across hundreds of cancer cell lines, revealing unexpected mechanisms (e.g., disulfiram-like and other repositioning signals) [PMID 32613204].
- **DepMap**: genome-scale CRISPR/RNAi dependency mapping connecting genes to cancer vulnerabilities [PMID 28753430] — primarily oncology, but its infrastructure (cell-line omics + effect matrices) seeded generalizable repurposing analytics.

### What it has actually produced for neuro diseases
- **Clemastine**: a high-throughput, expression/gene-signature-guided micropillar-array screen identified the antihistamine muscarinic antagonist as a pro-remyelinating agent in MS [PMID 24997607]; it subsequently showed delayed conduction improvement in chronic demyelination in the randomized ReBUILD trial [PMID 29029896]. This is the canonical small-molecule repurposing win adjacent to neurodegeneration — though it came from a functional screen with signature support, not pure CMap reversal.
- **Ropinirole** for ALS (iPSC-derived neuron phenotypic matching, Section 1) [PMID 37267913].
- Many published CMap-for-Alzheimer's/Huntington's/Parkinson's papers exist; almost none progressed past retrospective validation.

### Limitations — why signature matching disappoints
1. **Assay artifacts and context dependence**: signatures are cell-line-, dose-, and time-point-specific. A compound's L1000 profile at one dose in one line may invert at another.
2. **The bexarotene cautionary tale**: a signature-flavored systems argument ("APOE up → amyloid clearance") looked spectacular in mice [PMID 22323736], then failed independent replication within a year — no cognitive benefit, no consistent amyloid reduction across labs [PMID 23764200, 23704556, 26025659]. Systems-biology plausibility ≠ efficacy.
3. **Bulk signatures dilute cell-type signal**: an astrocyte disease program buried in whole-cortex RNA is often unrecoverable by simple reversal.
4. **Neuro-specific gaps**: L1000 cells are peripheral/cancer lines; BBB penetration and glia/neuron-specific pharmacology are absent from the prior entirely.

**Bottom line:** signature-matching is a cheap idea generator. Its genuine wins in neuro (clemastine, ropinirole) came when signature logic was fused to a *functional* assay with translational pharmacology — not when used standalone.

---

## 3. Structure-Based AI: AlphaFold3/Boltz-2, IDRs, Diffusion Models, and Oligonucleotide Design

### What AF2/AF3/Boltz-2 do well
AlphaFold2 solved monomer structure prediction [PMID 34265844]; AlphaFold3 extended to complexes including proteins, nucleic acids, ligands, and ions [PMID 38718835, addendum PMID 39604737]. Boltz-2 (MIT/Recursion lineage) is the leading *open* co-folding model and adds a binding-affinity head usable at virtual-screen throughput — its preprint (Passaro, Corso, Wohlwend et al., bioRxiv 2025, EPMC PPR1039145) reports near-docking-free affinity ranking; independent benchmarks already probe its reliability for structure and affinity prediction (e.g., PMID 42579383) and for classifying docking hits (e.g., PMID 41592323).

### The key caveat: TDP-43 and FUS are intrinsically disordered proteins (IDPs)
This is the single most important technical caveat for applying structural AI to ALS:
- Low pLDDT regions of AlphaFold predictions indicate **disorder, not unresolved structure** — they should never be interpreted as folds [PMID 40937679, 40131945].
- AlphaFold's amyloid-region behavior is systematically misleading relative to real aggregation landscapes [PMID 34023402].
- AlphaFold-Multimer can partially capture IDP-mediated interactions and dynamics in favorable cases [PMID 39446390], but confidence must be read carefully.
- TDP-43's C-terminal low-complexity domain and FUS's prion-like domain form liquid droplets via weak multivalent interactions tuned by post-translational marks: ALS mutations disrupt α-helical-structure-mediated phase separation of TDP-43 [PMID 27545621], methylation tunes TDP-43 LLPS [PMID 32132204], phosphorylation disperses FUS condensates [PMID 28790177], and pathological FUS mutations accelerate a liquid-to-solid transition into fibrillar aggregates [PMID 26317470]; hnRNPA2 behaves analogously [PMID 29358076].

**Implication:** a static AF3 "structure" of the TDP-43 LCD is close to meaningless. The relevant object is an *ensemble* governed by condensate physics, RNA stoichiometry, and modification state. Structural AI helps here only indirectly — e.g., designing binders to the folded RRM domains or to transiently ordered motifs, or as priors inside MD/coarse-grained LLPS simulations. Any pitch of "AI solved TDP-43" is hype.

### Protein and molecule generation
- **RFdiffusion** enables atomically plausible de novo protein scaffolds and binder design [PMID 37433327], with antibody extensions emerging. Genuine capability — but binders must target structured epitopes; diffusing binders against IDP targets remains immature.
- **DiffDock** (Corso et al., ICLR 2023 — arXiv:2203.01765, no PMID): generative diffusion over ligand poses; strong blind-docking recall versus classical tools, but pose accuracy still lags on hard pockets. Follow-up work shows hybrid physics/data refinement materially improves results [PMID 42023722], and systematic studies warn that docking against raw AlphaFold models — side chains unrelaxed, no induced fit — produces high false-positive rates [PMID 37546760, 41075093].
- **Generative small-molecule design**: reinforcement-learning generation (REINVENT) dates to 2017 [PMID 29086083]; reviews catalog the field's expansion [PMID 39106790, 38384298]. The flagship end-to-end clinical proof point remains Insilico's rentosertib (ISM001-055, TNIK inhibitor designed with generative chemistry + target discovery engine), which cleared Phase 2a in pulmonary fibrosis [PMID 40461817, 38459338] — notable, but fibrotic lung, not brain, and no CNS generative-design asset has replicated this.

### Antisense oligo and siRNA design ML
ALS-relevant because tofersen (SOD1 ASO), and the wider pipeline of ASOs for C9orf72, ATXN2, FUS, etc., makes oligo design a core industrial activity:
- Published ML aids are narrow but real: thermodynamic (Tm) prediction models for modified gapmer duplexes [PMID 39176173]; deep-learning siRNA efficacy predictors [PMID 30255786], graph neural networks such as siRNADiscovery [PMID 39503523], multimodal geometric predictors like ENsiRNA [PMID 40194620], and pairwise siRNA-mRNA attention architectures [PMID 42621926].
- Honest framing: modern oligo programs (Ionis-style) are dominated by medicinal chemistry rules (motifs, phosphorothioate placement, GalNAc conjugation) refined over decades, plus proprietary internal data. Public ML is a growing *adjunct*, not the design engine. Tofersen itself — including its presymptomatic deployment strategy [PMID 35585374] — was built on rational design, not deep learning.

---

## 4. Clinical Trial AI: Stratification, Digital Endpoints, Prognostic Modeling, Biomarkers

This is arguably where ML most clearly pays off in ALS today — because the disease's extreme heterogeneity inflates sample sizes and masks true drug effects.

### Data foundation: PRO-ACT
The Pooled Resource Open-Access ALS Clinical Trials database aggregated >8,000 placebo-arm patients from failed trials; the DREAM crowdsourced challenge on it demonstrated ML could predict progression better than individual teams' bespoke models [PMID 25362243]. PRO-ACT remains the field's public benchmark and the substrate for most published prognostic ML.

### Prognostic modeling and stratification
- The **ENCALS personalized survival model** (clinical + genetic features) is prospectively used in clinics and trials [PMID 29598923].
- Deep learning on MRI predicts survival beyond clinical features [PMID 28070484].
- Deployable variants exist for populations lacking full spirometry (VC-Free model) [PMID 34348539]; mortality-prediction ML validated across diverse databases [PMID 40952318]; deep sequence-style models predict ALSFRS-R trajectories [PMID 35962027]; dynamic updating of prognostic factors improves individual-level prediction [PMID 37014017].
- Methodologically important result: incorporating ML predictions as covariates in trial analysis **improves statistical power**, shrinking required sample sizes [PMID 32862509]; trial-innovation frameworks argue this is one of the highest-yield near-term fixes for ALS development [PMID 34315786], supported by formal validation of prognostic marker sets [PMID 39270623].
- **DIAAL-S**: a widely referenced web-based dynamic individual survival-prediction tool for ALS used in trial-design discussions. Note: our PubMed/Europe PMC/Crossref searches found **no indexed primary publication** under that name; cite it as a tool, and rely on ENCALS [PMID 29598923] and related dynamic models [PMID 37014017] for peer-reviewed anchors.
- Molecular subtyping (Section 1) feeds stratification: unsupervised ML-defined ALS subtypes offer a path to targeted enrollment [PMID 38129934].

### Biomarkers: NfL dynamics
Neurofilament light chain is the field's best-validated dynamic biomarker, and ML-on-NfL-trajectory is where biomarker science and trial AI converge:
- NfL rises before symptom onset in SOD1 carriers and predicts phenoconversion [PMID 30014505, 31432691].
- Meta-analyses confirm blood/CSF NfL as prognostic for survival across sporadic and familial ALS [PMID 34690913, 38674431]; broader biomarker panels extend this [PMID 41140053].
- The tofersen presymptomatic program operationalized NfL kinetics as a pharmacodynamic trigger — trial design explicitly powered around NfL change [PMID 35585374; context PMID 37382103]. That is a genuine AI-adjacent win: longitudinal biomarker modeling changed regulatory and clinical practice.

### Digital endpoints
- Wearable sensors objectively capture ALS symptoms and decline during trials at home [PMID 31859676].
- Multimodal **speech biomarkers** remotely track bulbar progression — validated as remote monitoring measures [PMID 38978682, 39126786] — with automatic speech-intelligibility scoring validated across languages and neurological diseases [PMID 39108340].
- Status: strong for enrichment and exploratory endpoints; FDA qualification for primary endpoints is still in progress. Do not oversell.

### ACT for ALS
The 2021 U.S. "ACT for ALS" law funds expanded-access investigation of experimental therapeutics for rapidly fatal diseases and supports research on outcomes in expanded-access populations. It is policy infrastructure (data generation from access cohorts), not an ML method per se; its value to trial AI depends on whether those datasets become analyzable. No robust PubMed-indexed evaluation exists yet — flagged as an open item.

---

## 5. Honest Assessment: Hype Versus Concrete Wins

### Context: Eroom's law persists
Pharmaceutical R&D efficiency has declined for decades — roughly nine-fold cost inflation per approved drug since the 1950s, attributed to the 'Better Than the Brill' effect, the cautious-regulator effect, throw-money-at-it tendency, and basic-research-brake traps [PMID 22378269]. Nothing in the AI era has yet bent that aggregate curve: approvals attributable to AI-designed molecules remain a handful globally.

### Quantified skepticism
- A quantitative review of AI-discovered small molecules entering the clinic found the count in the low dozens worldwide, with attrition patterns similar to conventional discovery — speed-to-candidate improved, not success rates [PMID 35132242].
- Medicinal chemists' assessments ask bluntly whether AI delivers beyond incremental automation [PMID 37738505].
- In Alzheimer's genetics specifically, ML adds classification polish more than new causal biology versus standard GWAS machinery [PMID 40691194].
- Structural-AI misuse is documented: docking into raw AlphaFold models inflates false positives [PMID 37546760, 41075093]; low-confidence regions encode disorder, not mystery folds [PMID 40937679].

### Concrete wins to defend
1. **Rentosertib** — first fully generative-AI target-plus-drug to clear Phase 2a [PMID 40461817, 38459338].
2. **Clemastine** — screen/signature-driven repurposing to positive Phase 2 remyelination signal [PMID 24997607, 29029896].
3. **Tofersen + NfL kinetics** — biomarker-model-driven presymptomatic intervention design [PMID 35585374, 30014505].
4. **PRO-ACT-derived prognostic ML reducing trial sample sizes** [PMID 25362243, 32862509, 29598923].
5. **Ropinirole** — iPSC-computational repurposing into ALS patients [PMID 37267913].

### Where hype has NOT delivered (neurodegeneration-specific)
- De novo generative design of CNS-penetrant drugs: zero approved assets; BBB pharmacology remains the bottleneck, not molecule ideation.
- AlphaFold-family models applied to TDP-43/FUS aggregation: category error against IDP/phase-separation biology [PMID 26317470, 27545621, 34023402].
- Standalone CMap reversal for AD/HD/PD: abundant papers, no clinic-validated winners; bexarotene shows how fast such stories collapse [PMID 22323736 vs 23764200, 26025659].
- ML-designed antisense oligos as headline technology: real but auxiliary [PMID 39176173, 39503523]; industrial oligo design remains rule-and-data driven.

### Net judgment (2024–2025)
For neurodegeneration/ALS specifically, AI/ML currently earns its keep in three places: (1) **causal, genetics-anchored target prioritization**; (2) **trial operations** — prognostic enrichment, stratification, digital and fluid biomarker dynamics; and (3) **engineering acceleration** (structure models for folded domains, co-folding screens like Boltz-2, oligo property models). Discovery-stage generative spectacle — new molecules, new folds — is real technology with as-yet thin neurodegeneration yield. Budget accordingly: fund the boring, measurable wins first.

---

*Verification note: all PMIDs above were retrieved and title-matched via NCBI E-utilities esearch/esummary during preparation of this brief. Preprint identifiers are given where no journal record exists (Boltz-2, DiffDock).*
