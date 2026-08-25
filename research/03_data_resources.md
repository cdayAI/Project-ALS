# Public Datasets & Computational Resources for Computational ALS Research
**Status date:** 2024–2025 · **Verification method:** every URL below was probed live (HTTP GET status codes) and/or queried via NCBI E-utilities during preparation of this brief. Status noted inline where useful.

---

## Executive orientation

ALS computational work today sits on three pillars: (1) large patient genomic cohorts (mostly controlled-access but obtainable in weeks-to-months), (2) an unusually generous open multi-omics program (Answer ALS), and (3) mature generic tooling (scRNA pipelines, perturbation-signature databases, structure prediction/docking). A single computational researcher can realistically do target-prioritization and drug-repurposing analyses **in months**, not years, by combining open tiers of these resources — provided they avoid getting blocked at dbGaP application queues by starting those applications early.

---

## 1. Patient / genomic data

### 1.1 Project MinE (whole-genome sequences)
- **What it is:** The largest ALS-specific WGS collection (~3,000+ genomes across releases; international consortium, Netherlands-led). Genotypes plus phenotypes (onset site, progression rate, C9orf72 status).
- **Access:** Controlled access through the European Genome-phenome Archive (EGA) under Data Access Committee agreements; apply via the Project MinE Data Sharing page. Historically managed via the "ALDFASS" system (legacy host no longer resolving as of this probe).
- **Verified URLs:** https://projectmine.com/ (HTTP 200) · https://projectmine.com/datasharing/ (HTTP 200)
- **Size:** ~30× WGS per sample; full release is tens of terabytes; summary-level and joint-called VCF subsets are far smaller.
- **Realistic use in months:** Rare variant burden testing in known ALS genes; replication of GWAS hits; ancestry-specific variant annotation; building a local gnomAD-style allele-frequency panel for ALS genes.
- **Caveat:** DAC approval typically takes weeks–months; budget that time up front.

### 1.2 Answer ALS (multi-omics + iPSC motor neurons) — the standout open resource
- **What it is:** ~1,000 well-phenotyped ALS and control subjects with longitudinal clinical data, from whom iPSC lines were made, differentiated into motor neurons/glia, and profiled with RNA-seq, whole-genome sequencing, DNA methylation, ATAC-seq, proteomics, and high-content imaging.
- **Access:** **Free and open** after a short registration on the data portal; bulk data hosted on Google Cloud Storage (`answer_als_open_data` bucket; console requires sign-in redirect, bucket exists as of probe).
- **Verified URLs:** https://answerals.org/ (HTTP 200) · https://dataportal.answerals.org/home (HTTP 200) · Google Cloud bucket `gs://answer_als_open_data`
- **Size:** Multi-petabyte raw tier; the curated expression/metadata tables are GB-scale and laptop-analyzable.
- **Realistic use in months:** Disease-vs-control differential expression in genetically stratified motor neurons; network/module analysis to nominate druggable nodes; cross-modality integration (methylation × transcriptome).

### 1.3 PACTALS (Pan-Asian Consortium for Treatment and Research in ALS)
- **What it is:** Asia-Pacific ALS registry/biorepository network harmonizing clinical phenotypes and biosamples across 10+ countries.
- **Access:** Consortium membership/collaboration agreements rather than a public download. Verified: https://pactals.org/ (HTTP 200); no public dataset endpoint found (`/pactals-dataset/` HTTP 404).
- **Use case:** Population diversity in ALS genetics — currently the weakest axis in existing cohorts; partner for East/Southeast Asian replication.

### 1.4 NEALS Biorepository
- **What it is:** Sample repository attached to the Northeast ALS Clinical Research Network (~100+ clinical trial sites): CSF, blood, tissue from trial participants, including negative-control placebo-arm samples.
- **Access:** Written proposal to NEALS (free; reviewed by committee); samples shipped or data shared per agreement.
- **Verified URL:** https://neals.org/als-researchers/neals-sample-repository/ (HTTP 200; older `/living-with-als/...biorepository` path returns 404)
- **Realistic use in months:** Biomarker validation on placebo-arm CSF (e.g., NfL, inflammatory panels) — wet-lab dependent; computationally, use NEALS trial data summaries for modeling progression rates.

### 1.5 NYGC ALS Consortium sequencing data
- **What it is:** New York Genome Center ALS Consortium (Phatnani lab): thousands of whole genomes from sporadic/familial ALS, plus matched RNA-seq from postmortem tissue; source of several landmark rare-variant papers (e.g., TBK1, KIF5A).
- **Access:** Data deposited publicly: NCBI BioProject **PRJNA573105** ("NYGC ALS Consortium data", verified HTTP 200 via https://www.ncbi.nlm.nih.gov/bioproject/PRJNA573105); newer releases flow through consortium collaboration agreements and AnVIL/Terra workspaces.
- **Note:** The landing page `nygenome.org/research/als-consortium/` returned HTTP 404 at probe time; use the news pages and BioProject instead.

### 1.6 GWAS Catalog — ALS hits
- **What it is:** Curated SNP-trait associations; ALS maps to EFO term `EFO_0000253`. Includes the big 2018–2022 ALS GWAS loci (C9orf72, UNC13A, TBK1, SOD1, etc. plus ~25 genome-wide significant loci).
- **Access:** Fully open REST API and downloads.
- **Verified URLs:** https://www.ebi.ac.uk/gwas/studies/EFO_0000253 (HTTP 200) · https://www.ebi.ac.uk/gwas/home (HTTP 200)
- **Realistic use in months:** Instant input for gene-set enrichment, chromatin mapping, colocalization with eQTL/sQTL (UNC13A cryptic exon biology), and Mendelian-randomization-style target triage.

### 1.7 ALS–FTD cohorts
- **GENFI (Genetic FTD Initiative):** GRN/C9orf72/MAPT carrier families; imaging + biomarkers. Verified: https://www.genfi.org/ (HTTP 200). Access via data-sharing request.
- **CReATe / RDCRN (Clinical Research in ALS and Related Disorders for Therapeutic Development):** US network, phenotype-genotype bank incl. ALS-FTD spectrum. Verified: https://create.rarediseasesnetwork.org/ (HTTP 200 after redirect).

### 1.8 UK Biobank
- **What it is:** ~500,000 participants, exome + full WGS release, linked health records; contains hundreds of incident ALS cases — small for ALS alone but uniquely deep for comorbidity, exposure, and generalizable method development.
- **Access:** Institutional academic application (modest cost recovery fee); analysis on the UK Biobank Research Analysis Platform (Terra/DNAnexus). Verified: https://ukbiobank.dnanexus.com/ (HTTP 200); main site https://www.ukbiobank.ac.uk/ returned HTTP 403 to our automated probe (bot-blocking, site is live in browsers).
- **Realistic use in months:** Phenome-wide scans for ALS-associated traits/prodromal features; polygenic risk score construction and calibration.

---

## 2. Molecular / cellular data

### 2.1 Single-cell atlases of human spinal cord / motor cortex
- **GSE221692** — "single cell RNA-seq and Visium data of human spinal cord" (n=49; verified via E-utilities esummary, accession resolves). Good scaffold for cell-type deconvolution of bulk ALS data.
- **Maniatis et al.–style postmortem spinal cord snRNA-seq** (2019, Nature Neuroscience) — foundational ALS spinal-cord atlas; find via GEO query `"spinal cord"[Title] AND "amyotrophic lateral sclerosis"[All Fields]`.
- **CELLxGENE Discover (CZI):** web-hosted, API-accessible census of >100M cells; multiple adult/fetal human spinal-cord collections. Verified: https://cellxgene.cziscience.com/ (HTTP 200) and `/collections?q=spinal+cord` (HTTP 200).
- **Allen Brain Map / BRAIN Initiative cell census:** human motor cortex taxonomy (human multiple cortices atlas, ~1M nuclei). Verified: https://brain-map.org/ (HTTP 200).
- **Human Cell Atlas:** https://www.humancellatlas.org/ (HTTP 200).
- **Realistic use in months:** Build cell-type marker matrices; project ALS bulk RNA-seq onto healthy atlas with scVI/deconvolution; identify motor-neuron-selective expression for target filtering.

### 2.2 iPSC-derived motor neuron datasets (verified GEO accessions)
All confirmed live via NCBI GDS esummary during this brief:

| Accession | Content | n samples |
|---|---|---|
| GSE299997 | Integrated profiling of iPSC-motor neurons carrying **C9orf72, FUS, TARDBP, SOD1** mutations | 156 |
| GSE283507 | **Ropinirole** improves TDP-43-mutant iPSC-motor-neuron phenotype (repurposing precedent) | 54 |
| GSE303931 | Transcriptomics across multiple **isogenic C9orf72** iPSC-derived neuronal lines | 12 |
| GSE313074 | C9orf72 repeat expansion rewires 3D chromatin in ALS | 39 |
| GSE229095 | Rapid induction of human spinal lower motor neurons; ALS cell screening | 15 |
| GSE284339 | Organoid-derived microglia; microglial dysfunction in C9orf72 ALS/FTD | 25 |

URL pattern: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ACCESSION>` (GEO root verified HTTP 200).

### 2.3 TDP-43 proteomics / interactomes
- **BioGRID** aggregates TDP-43 (TARDBP) physical/genetic interactions including BioID/PromID proximity datasets (Freibaum et al. 2015; Conlon et al. RNA-binding proximity screens). Verified: https://thebiogrid.org/ (HTTP 200).
- Use pattern: extract TARDBP neighborhood → intersect with ALS GWAS genes and LINCS reversal signatures → rank candidate nodes.

### 2.4 LINCS / Connectivity Map (CMap)
- **What it is:** ~1.5M+ perturbation-expression signatures (drugs × cell types, incl. neural lineages) built to be reversed against disease signatures.
- **Access:** Free account at clue.io; bulk downloads available. Verified: https://clue.io/ (HTTP 200); LINCS program data index https://lincsproject.org/LINCS/data (HTTP 200).
- **Realistic use in months:** Compute ALS disease signature (from Answer ALS / GSE299997), run CMap reversal, get a ranked list of candidate repurposing compounds in days.

### 2.5 PRISM Repurposing Dataset
- **What it is:** Broad Institute cancer-cell-line drug-sensitivity screen (~450 FDA-approved/in-clinical compounds × ~600 cell lines, 19Q4 secondary screen released openly).
- **Access:** Open download via DepMap repurposing hub. Verified: https://depmap.org/portal/page/repurposing (HTTP 200) and download listing `?repo=repurposing` (HTTP 200).
- **Caveat for ALS:** cancer-lineage sensitivity ≠ CNS efficacy; use PRISM mainly for compound-annotation and target-linkage metadata feeding into a motor-neuron-based assay plan.

### 2.6 DepMap
- **What it is:** Cancer dependency map: CRISPR/RNAi essentiality + omics across ~1,800 cell lines.
- **Access:** Open bulk downloads. Verified: https://depmap.org/portal/ (HTTP 200).
- **ALS angle:** validate that candidate targets' loss is tolerable in non-transformed contexts, mine co-dependency clusters containing ALS genes.

---

## 3. Models & tools

| Tool | Purpose | Verified URL | Notes |
|---|---|---|---|
| **AlphaFold DB** | Pre-computed structures, >214M proteins | https://alphafold.ebi.ac.uk/ (200) · /download (200) | Free bulk download; use for ALS protein variant impact & pocket inspection |
| **OpenFold** | Open reproduction/training code for AlphaFold | https://github.com/aqlaboratory/openfold (200) | Run your own predictions; training data nuances vs AF3 |
| **Boltz-2** | Co-folding + **affinity prediction** (open successor spirit to Boltz-1; strong 2025 baseline) | https://github.com/jwohlwend/boltz (200) | GPU needed for inference; good for compound-target docking triage |
| **RDKit** | Cheminformatics standard library | https://www.rdkit.org/ (200) | Free/BSD; fingerprinting, ADMET-ish descriptors for repurposing filters |
| **DiffDock** | Diffusion-based molecular docking | https://github.com/gcorso/DiffDock (200) | Blind docking without binding-pocket priors |
| **CellProfiler** | Image-analysis pipeline for high-content screens | https://cellprofiler.org/ (200) | Headless mode for cluster processing of iPSC-MN imaging |
| **scanpy** | Python scRNA-seq toolkit | https://scanpy.readthedocs.io/en/stable/ (200) | Standard pipeline backbone |
| **scVI-tools** | Probabilistic integration models (scVI/scANVI/totalVI) | https://scvi-tools.org/ (200) | Batch correction across donors/labs — critical for Answer ALS |

### ALS-focused computational resources
- **ALS Online Database (ALSoD)** — UCL-curated catalog of ALS genes, variants, and mutation frequencies. Verified: https://alsod.ac.uk/ (HTTP 200; some subpaths like `/genes`, `/about/` return 404 — navigate from homepage). Best single reference for the current ALS gene list.
- **mNDX:** we could **not verify** any established, live resource under this name (no DNS resolution for candidate hosts; no NCBI hits). If you encountered "mNDX," treat as unconfirmed or possibly superseded/renamed. Do not build plans around it without further evidence.
- Also useful: **GWAS Catalog ALS trait page** (above) doubles as the canonical ALS-gene shortlist source.

---

## 4. Access model & feasibility cheat sheet

| Resource | Access model | Typical size | Months-scale deliverable |
|---|---|---|---|
| Answer ALS | Free registration → open cloud data | GB (curated) – PB (raw) | Cross- genotype DE + network target nomination |
| Project MinE | EGA DAC application (controlled) | ~TBs | Rare-variant burden test in ALS genes |
| NYGC ALS WGS | Public BioProject deposit / consortium agreement | TBs | Meta-analysis replication |
| GWAS Catalog | Open | MBs | Colocalization & MR target triage in days |
| UK Biobank | Institutional application + platform fees | PB (platform-side) | PRS + prodrome phenome scan |
| LINCS/CMap | Free account | GBs | Compound-reversal ranking in days |
| PRISM | Open download | ~GBs | Compound/target annotation layer |
| DepMap | Open | GBs | Essentiality/tolerability filter |
| CELLxGENE atlases | Open (API) | GBs | Cell-type specificity scoring |
| AlphaFold DB / Boltz / DiffDock | Open | TBs (AF DB) / local GPU | Structure-informed variant & pocket triage |

Practical guidance:
1. Start controlled-access applications (Project MinE, UKB) **immediately**; do open-data work while waiting.
2. All open tiers above are analyzable on a single workstation; only raw WGS and AF DB bulk need object storage/compute.
3. The fastest defensible 3-month arc: Answer ALS signature → CMap/PRISM reversal → DepMap/ALSoD/GWAS Catalog filtering → Boltz-2/DiffDock structural check on top candidates.

---

## 5. Ranked shortlist — top 5 highest-leverage resources

**1. Answer ALS open multi-omics portal** — https://dataportal.answerals.org/home
The only place where clinically deep ALS cohorts meet iPSC-motor-neuron multi-omics with zero access friction. A disease signature derived here is the seed for everything downstream (targets, compounds, biomarkers).

**2. LINCS Connectivity Map (clue.io) + PRISM repurposing annotations** — https://clue.io/ · https://depmap.org/portal/page/repurposing
Reversing an ALS signature against ~1.5M perturbation signatures yields a concrete, testable repurposing hypothesis list within days. This is the shortest path from computation to a therapeutic claim.

**3. GWAS Catalog ALS associations (EFO_0000253) + Project MinE WGS** — https://www.ebi.ac.uk/gwas/studies/EFO_0000253 · https://projectmine.com/datasharing/
Human genetics remains the strongest causal-evidence engine for target selection. GWAS loci give directionality; MinE/NYGC genomes let you test rare-variant support before committing to a target.

**4. Human spinal-cord/motor-cortex single-cell atlases (CELLxGENE + GSE221692 + Allen Brain Map)** — https://cellxgene.cziscience.com/ · https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221692 · https://brain-map.org/
Every nominated target must pass a cell-type-specificity and druggability filter. Atlases provide the motor-neuron/glia selectivity map cheaply and immediately.

**5. Structure stack: AlphaFold DB + Boltz-2 + DiffDock (+ RDKit)** — https://alphafold.ebi.ac.uk/ · https://github.com/jwohlwend/boltz · https://github.com/gcorso/DiffDock
Converts genetic/network hits into mechanism-level hypotheses (variant destabilization, pocket compatibility, off-target profiles) without a crystallography pipeline.

*Justification logic:* items 1–2 generate hypotheses fast and openly; item 3 adds causal human evidence; item 4 enforces biological specificity; item 5 adds mechanistic plausibility. Together they cover the full target→compound→validation-prioritization loop inside a few months for one computational researcher.

---

## Verification log (probe snapshot)

Live (HTTP 200): projectmine.com · projectmine.com/datasharing/ · answerals.org · dataportal.answerals.org/home · pactals.org · neals.org/als-researchers/neals-sample-repository/ · ncbi.nlm.nih.gov/bioproject/PRJNA573105 · ebi.ac.uk/gwas/studies/EFO_0000253 · create.rarediseasesnetwork.org · genfi.org · ukbiobank.dnanexus.com · cellxgene.cziscience.com · brain-map.org · humancellatlas.org · thebiogrid.org · clue.io · lincsproject.org/LINCS/data · depmap.org/portal (+ repurposing page) · alphafold.ebi.ac.uk · github.com/aqlaboratory/openfold · github.com/jwohlwend/boltz · github.com/gcorso/DiffDock · rdkit.org · cellprofiler.org · scanpy.readthedocs.io · scvi-tools.org · alsod.ac.uk · geo.ncbi.nlm.nih.gov (root) · all cited GEO accessions resolved via E-utilities.

Not resolvable / moved / blocked: aldfass.projectmine.com (DNS failure) · nygenome.org/research/als-consortium/ (404) · ukbiobank.ac.uk (403 to bots; live) · ega study EGAD00010000904 (404 — cite Project MinE datasharing page for current EGA studies instead) · "mNDX" (unverifiable).
