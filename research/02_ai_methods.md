# AI Methods for Computational ALS Drug Discovery
**A critical brief on the computational methods available to Project ALS's HYPOTHESIZE -> EXPERIMENT -> ATTACK -> TRIAGE loop (docs_plan.md)**

*Sources retrieved via PubMed (NCBI); PMIDs cited inline. Method papers are drawn from ALS literature where available and from adjacent neurodegenerative-disease / oncology literature where ALS-specific validation doesn't yet exist — those cases are flagged, since a method proven on Alzheimer's or cancer data carries real but unproven transfer risk to ALS.*

---

## 1. Transcriptome-based drug repurposing — the Pilot Sprint #1 core method

This is the method class `docs_plan.md` bets Pilot Sprint #1 on: derive an ALS disease signature, reverse it against a compound-perturbation database, rank candidates.

**Directly on ALS:**
- A 2026 consensus machine-learning study derived recurrent ALS transcriptional signatures from motor-cortex (E-MTAB-2325) and blood (E-TABM-940) datasets using four feature-selection methods across 100 repetitions of 4-fold cross-validation, then ran Connectivity Map drug-signature analysis on the recurrently selected probes. No single gene/probe was shared between the top motor-cortex and blood signatures, but pathway-level integration converged on glial/immune regulation, proteostasis, MAPK stress signaling, and RNA processing — the same mechanism list as `research/01_biology_and_therapeutics.md`. Deferoxamine and disulfiram showed reversal-compatible profiles in motor cortex; yohimbic acid, atovaquone, ciprofloxacin, and prochlorperazine in blood (PMID 42600917).
- A protein-inhibition-focused ML study trained per-target models for three ALS-associated kinases (Casein kinase 1, PTK2, EPHA4), screened FDA-approved drugs against all three, and validated hits by protein-ligand docking. Risperidone was the top multi-target hit (docking affinity -8.9 kcal/mol), with a strong ML-score/docking-affinity correlation supporting the pipeline's reliability (PMID 40745959).
- A network-biology + ML study meta-analyzed 4 ALS transcriptome datasets (motor neuron and muscle tissue separately), built differential co-expression modules per tissue, used KNN/SVM/Random Forest to confirm the modules discriminate disease from control, then ran drug-repurposing + text-mining on the module genes as targets — nominating nilotinib, trovafloxacin, apratoxin A, carboplatin, and clinafloxacin (PMID 40180646).
- A cross-disease (AD/ALS/MS/PD) study combined ML classification of mitochondrial gene expression with Mendelian Randomization for causal gene-disease links, then ran drug enrichment + docking; it flags a pyruvate-metabolism gene as a shared AD/ALS vulnerability and nominates celecoxib as a repurposing candidate (PMID 40864790) — mitochondrial dysfunction is independently flagged as "moderate-strong" evidence in `research/01_biology_and_therapeutics.md` §1.7.
- A TDP-43-focused repurposing pipeline (DRIAD-SP) updated to use cryptic-exon detection as a pTDP-43 proxy found the JAK inhibitors baricitinib and ruxolitinib protective in AD brain regions with multiple cryptic exons, and identified TYK2 as a CRISPR-screen hit for cdsRNA-mediated neural death — this was built and tested primarily on AD but explicitly targets the TDP-43/interferon mechanism shared with ALS (PMID 38895380).

**Methodological grounding (not ALS-specific, underpins the above):**
- A systematic reconciliation of connectivity/reversal scoring methods (ES, XSum, RGES, Tau, CSS, EMUDRA, etc.) — useful as a reference when choosing which CMap/L1000 scoring function the pipeline should use, since the field has not converged on one (PMID 34013329).
- A 2026 review of gene-expression-guided repurposing methodology (signature reversion, pathway-level analysis, validation studies), with the caveat — directly relevant to the pilot — that "transcriptomic responses are highly context-dependent" and off-target effects "complicate mechanism interpretation" (PMID 41751568).

**Cross-disease review:** A broader review of ML for drug combination optimization across neurodegenerative diseases (AD, PD, Huntington's, MSA, ALS) surveys SVM, CNN, RNN, and transformer approaches for virtual screening/repurposing/combination design (PMID 41502663).

---

## 2. Network-based target prioritization

Complementary to signature reversal: instead of (or alongside) reversing expression, propagate ALS genetic/omics evidence across a molecular interaction network to nominate targets not necessarily captured by differential expression alone.

- A co-essentiality network built from CRISPR screens across 769 cancer cell lines outperformed conventional PPI/molecular networks at prioritizing anticancer drug targets, including nominating 3 approved drugs repurposed to a new indication that other networks missed (PMID 41002174) — the general lesson (network choice materially changes repurposing output) applies directly to any ALS network-propagation step.
- A multi-layer network (PPI + gene regulation + metabolite interactions + multiple disease-signature types) combined with ML for target repositioning explicitly critiques over-reliance on PPI-only proximity as a source of bias, and demonstrates propagation-based, non-proximity features improve target discovery (prostate cancer case study) (PMID 38983753).
- Weighted gene co-expression network analysis (WGCNA) combined with Mendelian Randomization for causal gene prioritization and Connectivity Map for compound nomination, validated in a mouse model — this is close to a full closed loop (network -> causal filter -> CMap -> in vivo validation) and is a template worth adapting for ALS (Alzheimer's case study) (PMID 41409615).
- A hypergraph neural network approach to autism risk genes — modeling multi-gene complexes rather than pairwise PPI edges — found a 51% higher correlation with independent genetic evidence (TADA) than random-walk methods, and identified a druggability-filtered "super-hub" gene cluster (PMID 41595671). Relevant as a more expressive alternative to simple network propagation once the ALS GWAS/rare-variant gene set (`research/03_data_resources.md` §1.6) is available as propagation seeds.
- A network-and-genomics framework for the mTOR pathway (GeneCards/KEGG/STRING/UniProt/PathCards integration + personalized PageRank/Random Walk with Restart) is a directly reusable recipe for building an ALS-pathway-specific propagation network (PMID 41300706).

---

## 3. Structure-based / deep-learning methods for target and compound triage

Once candidate targets and compounds are nominated, these methods add mechanistic plausibility — the role AlphaFold DB/Boltz-2/DiffDock/RDKit play in `research/03_data_resources.md` §3.

- **Most directly relevant finding:** SKALE 2.0, a phase-resolved geometric deep-learning framework, was trained and validated explicitly on **SOD1, TDP-43, MAPT, and PRNP** — the core ALS aggregation-prone proteins from `research/01_biology_and_therapeutics.md` §1.1–1.3. It separates nucleation-phase from elongation-phase structural determinants of aggregation (something static AlphaFold-derived risk scores fail to resolve) and its predicted suppressor/enhancer/phase-switch mutations were validated in recombinant SOD1 experiments (PMID 42307331). This is a strong candidate tool for the design stream in `hypotheses/_TEMPLATE.md` (`streams: [design]`).
- A 2026 review of AI in protein misfolding and innate immunity across neurodegenerative disease surveys AlphaFold, I-TASSER, RoseTTAFold, Phyre2, and ESMFold applications to conformational modeling of misfolded proteins (Aβ, Tau, α-synuclein, **TDP-43**) and DAMP-mediated neuroinflammation — directly ties structure prediction to the neuroinflammation mechanism in `research/01_biology_and_therapeutics.md` §1.6 (PMID 42206050).
- AlphaFold-predicted structures were used to explain the biochemical effect of disease-linked mutations in TOE1 (a neurodegenerative-disease RNA-processing enzyme), demonstrating the "predict structure -> rationalize variant effect" workflow at protein-variant scale rather than proteome scale (PMID 42178111) — a template for interpreting ALS rare-variant hits from Project MinE/NYGC data (`research/03_data_resources.md` §1.1, §1.5).
- For compound-side triage: a deep-learning-rescoring virtual-screening pipeline (AutoDock Vina + GNINA 3D-CNN rescoring + MD simulation) screened ~180,000 compounds and reduced false-positive binding predictions relative to empirical scoring alone (PMID 42630437) — a concrete upgrade path beyond plain DiffDock/AutoDock docking.
- A review of geometric deep learning for drug design (graph neural nets, SE(3)-equivariant networks, geometric transformers) for binding-affinity prediction and de novo generation, useful as a technology-selection reference once the pilot moves from docking to learned scoring functions (PMID 42572360).
- One cautionary data point: a head-to-head comparison found LLM-assisted keyword literature/database search returned messier, less structurally-usable results than classical substructure search for compound discovery (PMID 42462355) — evidence against using an LLM as a substitute for structured chemical database queries in the pipeline.

---

## 4. Single-cell deconvolution and cell-type specificity

Needed to convert bulk ALS transcriptome signatures (used in §1) into cell-type-resolved evidence, and to project bulk signatures onto the atlases in `research/03_data_resources.md` §2.1.

No ALS-specific deconvolution methods paper surfaced in this search — the following are general-purpose methods with active development that a bulk-signature pipeline should benchmark against, not evidence about ALS biology itself:
- DeconX addresses a real gap in standard deconvolution: cell types entirely missing from the single-cell reference (e.g., due to poor dissociation) are invisible to standard methods; DeconX recovers their existence, proportion, and expression signature from deconvolution residuals (PMID 42635234) — relevant because motor neurons are notoriously fragile and under-represented in dissociation-based single-cell references, exactly the failure mode this method targets.
- SKIM, a fast dataset-sketching method with model dynamic-feedback, benchmarks explicitly against bulk RNA-seq deconvolution accuracy as one of its four validation tasks and reports >20x speedup at large scale (PMID 42574452) — relevant for keeping CELLxGENE-scale atlas work (`research/03_data_resources.md` §2.1) tractable on a single workstation.

---

## 5. LLM-based literature mining and hypothesis generation

Maps directly to Layer 1 of `docs_plan.md` ("Mine literature ... for contradictions, gaps, and cross-domain connections") and to the `evidence` / `Prior art checked` fields in `hypotheses/_TEMPLATE.md`.

- SKiM-GPT combines classical literature-based discovery (A-B-C term co-occurrence) with retrieval-augmented LLM evaluation: for each candidate hypothesis it retrieves the actual PubMed abstracts, filters for relevance, and has an LLM score agreement with a stated hypothesis, showing its work. On a 14-hypothesis disease-gene-drug benchmark it reached strong agreement with expert biologists (Cohen's κ = 0.84). Open source with a public web interface (PMID 41408154) — this is close to a drop-in tool for the "Prior art checked" step of the hypothesis ledger, and its transparency design (shows retrieved abstracts + justification) matches the project's falsifiability-first operating rules.
- eGoT combines automated knowledge-graph construction from a literature corpus with graph-of-thoughts multi-hop querying, outperforming other retrieval-augmented methods on multi-hop biomedical QA benchmarks; demonstrated finding a plausible disease-environment connection (Lupus and UV exposure) from a closed corpus (PMID 42412787) — a candidate approach if the project later wants multi-hop reasoning across mechanism + genetics + compound literature rather than single-hop retrieval.
- IID-KG built an open, ontology-aligned literature knowledge graph (from 30M+ PubMed abstracts and 1.4M+ PMC full-text articles) that ships with a drug-repurposing hypothesis-generation workflow (PMID 42244715) — not ALS-specific (infectious/immune-mediated disease), but the released methodology and pipeline are directly adaptable, and it's a working example of what a project-owned literature graph could look like.

**Caveat directly on point:** the SKiM-GPT paper is explicit that standalone LLMs (without retrieval grounding) are "hampered by hallucinations, lack of transparency in information sources, and the inability to reference data not included in the training corpus" (PMID 41408154) — an argument for building any Project ALS literature-mining tooling on retrieval-augmented, citation-transparent pipelines rather than raw LLM prompting, consistent with Operating Rule #1 ("No result counts until an adversarial review exists").

---

## 6. Honest constraints

1. **Cross-disease transfer risk.** Several of the strongest methodological papers here (WGCNA+CMap+MR in Alzheimer's, hypergraph propagation in autism, co-essentiality networks in cancer) were validated on non-ALS diseases. The methods generalize; the specific findings do not. Every method above needs re-validation on ALS-specific data before its output enters the hypothesis ledger as more than a starting point.
2. **No ALS-specific single-cell deconvolution literature was found** in this search — §4's tools are general-purpose and unvalidated on motor-neuron/spinal-cord data specifically. This is a real gap, not just an artifact of search terms (a second, broader query on "single cell RNA-seq deconvolution motor neuron disease" returned zero PubMed hits under MeSH-mapped terms).
3. **Positive controls remain mandatory** (Operating Rule #3): none of the ALS-specific repurposing papers in §1 report whether their pipeline rediscovers riluzole, edaravone, or tofersen as positive controls — this should be a required check before trusting any pipeline built from this brief on novel candidates.
4. **Repurposing hit lists disagree across studies.** PMID 42600917, PMID 40745959, and PMID 40180646 each used different data (motor cortex vs. blood vs. combined tissue meta-analysis vs. protein-target screening) and produced almost entirely non-overlapping drug candidate lists. This is expected given different signatures and methods, but it means no single paper's shortlist should be treated as consensus — convergent evidence across methods, not any one paper's ranked list, is the bar.

---

## 7. How this maps onto Pilot Sprint #1 (`docs_plan.md`)

| Pilot Sprint #1 step | Method class (section) | Best-supported reference |
|---|---|---|
| 1. Pull ALS-vs-control RNA-seq signature | — (data step, see `research/03_data_resources.md`) | Answer ALS / GEO accessions |
| 2. Compute disease signature; score compounds by reversal | §1 (transcriptome repurposing) + §1 methodology | PMID 42600917, PMID 34013329 |
| 3. Positive control check | §6 constraint #3 | *(gap — no paper found doing this for ALS; must self-implement)* |
| 4. Cross-check top candidates against safety/max-phase | §3 (structure/docking) + ChEMBL/DrugBank (`research/03_data_resources.md`) | PMID 42630437 |
| 5. Adversarial agent attacks the ranking | §5 (LLM literature grounding for "Prior art checked") + §6 constraints | PMID 41408154 |
| *(not yet in Pilot Sprint #1, but supported)* Cell-type specificity filter | §4 (deconvolution) | PMID 42635234 |
| *(not yet in Pilot Sprint #1, but supported)* Target-level network propagation ahead of/alongside signature reversal | §2 | PMID 41002174, PMID 41300706 |
