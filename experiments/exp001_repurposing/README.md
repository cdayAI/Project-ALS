# Experiment 001: Transcriptome-based drug repurposing screen for ALS

Status: **PENDING ADVERSARIAL REVIEW** (see VERDICT.md - not self-certified).

## Question
Can we find perturbagen signatures in LINCS L1000 Phase-1 that reverse a
disease signature derived from human ALS post-mortem CNS tissue?

## Data (all open)
| Dataset | Role | Files |
|---|---|---|
| GSE124439 (GEO) | Disease signature: 145 ALS-spectrum MND vs 17 non-neurological control post-mortem CNS samples (frontal + motor cortex), gene-level STAR counts | `data/gse124439/` |
| GSE92742 LINCS L1000 Phase-1 Level5 COMPZ.MODZ | Reversal reference: 473,647 consensus signatures x 12,328 inferred genes, small-molecule perturbagens | `data/lincs/*.gctx(.gz)` + sig/pert/gene info |
| ChEMBL REST / PubChem REST | Candidate annotation (max_phase, approval year, MW, H-bond donors) | fetched at runtime, cached in `outputs/annotation_cache.json` |

## Method
1. Parse GEO series matrix -> sample groups. Drop transposable-element count rows
   (~1k features, ~26% of raw library size) BEFORE normalization.
2. CPM filter (>1 CPM in >=17 samples), log2-CPM, OLS with design `[1, ALS]`,
   limma-style empirical-Bayes moderated t (inverse-gamma prior on residual
   variance, method-of-moments), BH correction.
3. Map moderated t to the 12,328 L1000 genes via `pr_gene_id`; zero elsewhere.
4. Score every Level5 signature by `-cosine(disease_vector, signature)`
   (`score_reversal`, larger = stronger reversal). Single streaming pass over
   signature-row blocks of the HDF5 matrix (~14 s).
5. Collapse signatures to drugs by best score; keep n_sigs as exposure measure.
   Sensitivity ranking restricted to drugs with >=3 signatures is also written.
6. Positive control: ranks of known neuro/immune-active compounds.
7. Annotate top-30 candidates with ChEMBL max_phase + PubChem properties;
   flag CNS-relevant drugs (curated list OR MW<450 & HBD<=3).

## Reproduce
```bash
# from repo root: create env once
uv venv .venv && uv pip install --python .venv/bin/python numpy scipy pandas pyyaml h5py requests
cd experiments/exp001_repurposing && ../../.venv/bin/python run.py --config config.yaml
```
The 21 GB GCTX is NOT in git (`/.gitignore`). Download:
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz`
(sha512 `6a3115cf...08a2a`, see `data/lincs/`). run.py decompresses it automatically.

## Outputs (`outputs/`)
- `sample_sheet.csv` - parsed GSE124439 metadata (group, subregion)
- `disease_signature_full.csv` - DE result for all genes (log2FC, moderated t, p, padj)
- `als_signature_top.csv` - top 50 up / 50 down genes
- `sanity_check.md` - marker-gene + pathway-level sanity check
- `all_signature_scores_trt_cp.csv` - every small-molecule signature score
- `drug_ranking.csv` - drug-level ranking (primary: best signature score)
- `drug_ranking_nsig3_sensitivity.csv` - sensitivity: drugs with >=3 signatures only
- `positive_controls.csv` - positive-control compound ranks
- `top_candidates_annotated.csv` - top-30 with ChEMBL/PubChem annotation
- `summary.json`, `run_log.txt`

## Key caveats (honest)
- Individual canonical marker genes are mostly ns in this contrast; the sanity
  signal lives at the pathway level (ER protein processing / UPR / autophagy).
- Best-signature-per-drug ranking has a multiple-testing advantage for drugs
  with many signatures; use the n_sigs>=3 sensitivity file when interpreting.
- No covariate adjustment for RIN/PMI/age/sex (not in the series matrix).
