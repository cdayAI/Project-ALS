# 00 — Synthesis: What Four Research Briefs Tell Us To Do
*Cross-brief synthesis of research/01-04. Every claim traceable to PMIDs in the source briefs.*

## The honest starting point
- No cure exists; approved drugs buy roughly 2-6 months combined (01).
- The 2024-25 pipeline is mostly negative: Relyvrio withdrawn (PHOENIX failed),
  reldesemtiv futile, CNM-Au8 missed, verdiperstat negative. Masitinib (Ph3,
  OS HR 0.53 in enriched cohort) is the one live oral candidate awaiting confirmation (01).
- Genetics-first knockdown works only where targets are crisp: SOD1/tofersen
  validates the paradigm for ~2% of patients; ATLAS shows presymptomatic NfL
  lowering alone did NOT delay onset - timing and target choice both matter (01).

## Where AI has EARNED trust (fund these)
1. **Causal target prioritization** - Mendelian randomization over the druggable
   genome + GWAS colocalization is the strongest validated method (02: PMIDs
   38019415, 38443977). Inputs exist and are open (03: GWAS Catalog EFO_0000253,
   Project MinE via EGA application).
2. **Trial stratification & enrichment** - ENCALS/PRO-ACT ML covariates boost
   power and cut sample sizes; NfL dynamics are the qualified biomarker story
   (02: 29598923, 25362243, 35585374). Also directly patient-relevant.
3. **Structure/co-folding engineering for FOLDED targets** - AlphaFold2/3,
   Boltz-2 with affinity heads are real tools (02). Constraint: TDP-43/FUS are
   intrinsically disordered - low pLDDT means disorder, not solved structure.
   "AI solved TDP-43" is a category error; LLPS physics governs them.

## Where AI is an idea generator ONLY (use cheaply, never trust directly)
- **L1000 signature reversal** (our Sprint #1) produced clemastine->ReBUILD Ph2
  (a win) but also bexarotene (failed replication within a year). The method's
  real wins fused signatures to FUNCTIONAL ASSAYS before advancing (02).
  => Treat exp001 output as a ranked hypothesis list requiring orthogonal validation.
- Generative small-molecule design: one non-CNS clinical asset globally
  (rentosertib); zero approved CNS assets (02). Not our near-term bet.

## The factory's hypothesis-stream ranking (informed by all briefs)
1. **Stream A - Causal targets**: MR/colocalization on ALS GWAS loci against
   druggable genome; cross-reference Answer ALS multi-omics + single-cell atlas
   cell-type expression. Highest evidence-backed value.
2. **Stream B - Repurposing screen** (RUNNING as exp001): ranked reversal list ->
   ChEMBL phase annotation -> flag CNS-penetrant, safety-known compounds ->
   propose functional-assay validation path for top hits.
3. **Stream C - Patient-facing tooling**: trial matcher (ClinicalTrials.gov API,
   already working) + progression modeling on PRO-ACT. Fastest human impact;
   also builds community credibility for the repo.
4. **Stream D - Structure-based design**: deferred until Streams A/B nominate a
   FOLDED, causally-supported target worth engineering against.

## Data arc (from 03+04, all access verified)
Answer ALS + GSE124439/GSE255602/GSE261875 -> LINCS GSE92742 open GCTX reversal
-> GWAS Catalog/MinE causal triage -> Boltz-2/DiffDock structural work on
surviving folded targets. Apply early to Project MinE DAC (weeks-months latency).

## What would change our mind
- Masitinib confirmatory results (validates enriched-cohort trial design).
- Any credible replication of signature-reversal hits WITHOUT assay fusion -
  would upgrade Stream B confidence.
- WVE-004 / C9orf72 ASO biomarker data (extends genetics-first paradigm beyond SOD1).
