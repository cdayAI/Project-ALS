# Data inventory — Project ALS

Local dataset registry. Raw data lives in `data/` (gitignored); only this `inventory/`
folder is committed. Checksums are SHA-256 (`shasum -a 256`).

## Local datasets

### data/gse124439/ — GSE124439 (owned by exp001 / sprint1-repurposing agent)
- Human post-mortem CNS bulk RNA-seq: 145 ALS-spectrum MND vs 17 non-neurological control
  (+14 other neuro), frontal + motor cortex.
- Status: downloaded and stewarded by the exp001 agent. Listed here for completeness only —
  do not modify without coordinating with the exp001 owner.
- Serves: exp001 primary disease signature.

### data/gse255602/ — GSE255602 (stewarded by pipeline/data-steward)
- Sporadic ALS iPSC-derived spinal-cord-chip motor neurons vs isogenic controls, bulk RNA-seq
  salmon gene counts (ENSG).
- Files:
  | File | Bytes | SHA-256 |
  |---|---|---|
  | GSE255602_Bulk_RNAseq_all_samples_salmon_counts.csv.gz | 5,972,463 | `1920716ded3780130c6708814cf9dbe8e56a13449f1cb2d25bb4c765355f4f2c` |
  | GSE255602_series_matrix.txt.gz | 6,501 | `4abc1dbd8d9e3142fb53b980fcf92dc49e3166f9d62eb0c820d6f63975e03e0b` |
- Metadata: `inventory/gse255602_metadata.csv` — 71 samples, all counts columns mapped 1:1 to
  GSM accessions. Labels: 27 ALS / 44 control across 9 cell lines, experiments 1.3-5.4,
  Chip(Flow/Static) + 96-well formats.
- Source URLs:
  - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255602/suppl/GSE255602_Bulk_RNAseq_all_samples_salmon_counts.csv.gz
  - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE255nnn/GSE255602/matrix/GSE255602_series_matrix.txt.gz
- On disk size: 5.7M
- Serves: replication cohort for exp001 signature reversal (iPSC human-cell context).

### data/gse261875/ — GSE261875 (stewarded by pipeline/data-steward)
- TDP-43 perturbation in human iPSC-derived motor neurons, bulk RNA-seq unnormalized counts.
- Design: 2x4 factorial x 4 bio-reps = 32 samples, Day 12.
  Constructs: empty vector / YFP only / TDP-43:YFP / TDP-43-deltaNLS:YFP.
  Backgrounds: wild type / ATXN2 knockout.
- Files:
  | File | Bytes | SHA-256 |
  |---|---|---|
  | GSE261875_counts_unnorm.txt.gz | 2,010,639 | `12e2c5ac5f08d905d864c7d4a1e54512ce4ee6ee1cc6c08ab90f7de3e09f48f6` |
  | GSE261875_series_matrix.txt.gz | 5,231 | `f145510c6aefa277b055886105ebb5fc52ddfd96168bd8f728414b044477dff5` |
- Metadata: `inventory/gse261875_metadata.csv` — 32 samples; `counts_column` maps each GSM to
  its column tag in the counts file (verified as a multiset match).
- Source URLs:
  - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261875/suppl/GSE261875_counts_unnorm.txt.gz
  - https://ftp.ncbi.nlm.nih.gov/geo/series/GSE261nnn/GSE261875/matrix/GSE261875_series_matrix.txt.gz
- On disk size: 1.9M
- Serves: orthogonal validation of exp001 hits via TDP-43 mechanism axis (AGENTS.md review
  checklist explicitly requires replication vs GSE255602/GSE261875).

### data/lincs/ — LINCS L1000 (owned by exp001 / sprint1-repurposing agent)
- GSE92742 Broad L1000 Phase 1 assets per research/04_pilot_access_report.md.
- Status: owned by exp001 agent; listed for completeness only.
- Serves: exp001 signature-reversal query space (~473k Level5 signatures).

## Reference-only resources (not downloaded)

See [single_cell_atlases.md](single_cell_atlases.md) for verified access paths to
GSE221692 (human spinal cord scRNA+Visium) and CELLxGENE ALS collections.

## Experiment-stream map

| Dataset | Stream / consumer |
|---|---|
| GSE124439 | exp001 primary disease signature |
| GSE255602 | exp001 replication (iPSC MN) |
| GSE261875 | exp001 orthogonal validation (TDP-43 axis) |
| LINCS GSE92742 | exp001 reversal library |
| GSE221692, CELLxGENE collections (access paths only) | Stream A scaffolding / future single-cell streams |

## Housekeeping

- `.gitignore` uses `data/*` with `!data/inventory/` so only inventory docs are committed.
- Never commit raw files >50 MB (AGENTS.md hard rule #1).
