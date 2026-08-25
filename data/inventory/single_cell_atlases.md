# Single-cell atlases — verified access paths

Verified live (HTTP status in parentheses) by data-steward on the pipeline/data-steward branch.
No raw atlas files are stored locally yet; these are access paths for future streams.

## GSE221692 — human spinal cord scRNA-seq + Visium (embryonic development)

- Title: "single cell RNA-seq and Visium data of human spinal cord" (PMID 38177242)
- Design: 10x scRNA-seq, 10 time points x 3 segments (C/T/L), plus 10x Visium spatial at 4 time points
- GEO page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221692 (200)
- Bulk archive: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221692/suppl/GSE221692_RAW.tar (200; ~2.2 GB per filelist.txt)
- Per-file listing: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221692/suppl/filelist.txt (200)
- Format: 10x MTX triplets (barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz) per sample
- Note: embryonic developmental atlas — useful as reference for motor-neuron / progenitor states,
  not an ALS-vs-control contrast.
- Serves: Stream A (target triage context: developmental cell-state references); optional
  deconvolution reference for bulk ALS signatures.

## CELLxGENE (CZ CELLxGENE Discover) — ALS-relevant collections

Browse: https://cellxgene.cziscience.com/ (200) | Census API docs:
https://chanzuckerberg.github.io/cellxgene-census/ (200)

Collections API used for verification:
`GET https://api.cellxgene.cziscience.com/curation/v1/collections` (200)

### 1. ALS motor cortex and spinal cord single-nucleus multiome dataset
- Collection page: https://cellxgene.cziscience.com/collections/0986e4cd-7a58-405d-9b91-4b199bb4124e (200)
- collection_id: `0986e4cd-7a58-405d-9b91-4b199bb4124e`
- snRNA-seq + snATAC-seq, ALS vs control, cell-type annotated
  - Primary motor cortex: ~75,583 cells (dataset e2a00644...)
  - Lumbar spinal cord: ~62,711 cells (dataset b57462e3...)
- Highest-priority human ALS single-cell contrast for replication/validation.

### 2. Cortical brain samples from C9-ALS, C9-ALS/FTD, C9-FTD patients and age matched controls
- Collection page: https://cellxgene.cziscience.com/collections/aee9c366-f2fb-470b-8937-577d5d87d3fc (200)
- collection_id: `aee9c366-f2fb-470b-8937-577d5d87d3fc`
- Frontal + occipital cortex nuclei (~463k total across 4 datasets); C9orf72 repeat expansion
  (genetic ALS) context; useful for TDP-43/C9 mechanism checks.

### 3. Population-scale cross-disorder atlas of the human prefrontal cortex
- Collection page: https://cellxgene.cziscience.com/collections/84ce6837-548d-4a1f-919f-0bc0d9a3952f (200)
- collection_id: `84ce6837-548d-4a1f-919f-0bc0d9a3952f`
- ~3.9M nuclei dorsolateral prefrontal cortex, cross-disorder incl. neurodegenerative cohorts;
  background/reference atlas for cell-composition shifts.

## Download notes

- CELLxGENE datasets download as H5AD/CSV via the Census API (`cellxgene-census` Python package)
  or direct dataset URLs under https://cellxgene.cziscience.com/datasets/<dataset_id>.h5ad
- GSE221692 RAW.tar is >50 MB: download to `data/gse221692/` (gitignored); never commit it.
