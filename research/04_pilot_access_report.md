# ALS Drug-Repurposing Pilot — End-to-End Data Access Report

Tested live with Python `requests` from this machine (macOS). Status codes are exact HTTP responses observed on test date.

---

## 1. Disease data — NCBI GEO (READY TODAY, no auth)

### API
| What | Result |
|---|---|
| `GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=...` | **200** |
| `GET .../esummary.fcgi?db=gds&id=...` | **200** |
| Auth | None |

### Verified datasets (ALS vs control, human, bulk RNA-seq)
| GSE | Samples | Design | Supplementary verified |
|---|---|---|---|
| **GSE124439** | 176 (145 ALS-spectrum MND / 17 non-neurological control / 14 other neuro) | Human post-mortem CNS (frontal + motor cortex), Illumina HiSeq | `matrix/GSE124439_series_matrix.txt.gz` → **HTTP 200**, downloaded fully (11,491 B), parsed: sample groups confirmed from `!Sample_characteristics` |
| **GSE255602** | 71 (sporadic ALS iPSC-derived spinal-cord-chip motor neurons vs isogenic controls) | Bulk RNA-seq, salmon counts | `geo/download/?acc=GSE255602&format=file&file=GSE255602_Bulk_RNAseq_all_samples_salmon_counts.csv.gz` → **HTTP 200**, 5,972,463 bytes, valid gzip, gene-level ENSG count matrix parsed to completion |
| **GSE261875** | 32 | TDP-43 perturbation in human iPSC-derived motor neuron systems, bulk RNA-seq | Suppl listing **HTTP 200**: `GSE261875_counts_unnorm.txt.gz`; series matrix **HTTP 200** (HEAD, 5,231 B) |

Notes:
- Direct FTP-mirror path (`https://ftp.ncbi.nlm.nih.gov/geo/series/...`) works with plain GET.
- The `https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE...&format=file&file=<name>` endpoint returns **404 for files not in `suppl/`** (e.g., a series matrix name passed to it failed), but works correctly for actual supplementary files (verified: full 5.97 MB counts download). Use the ftp.ncbi.nlm.nih.gov mirror path for matrices.
- Also available if wanted: GSE330130 (human lumbar spinal cord + motor cortex, n=87; single-nucleus).

## 2. LINCS L1000 / CMap

| Resource | URL tested | Status | Auth | Verdict |
|---|---|---|---|---|
| clue.io app | https://clue.io/ | 200 | account needed to use | FREE REGISTRATION |
| clue.io REST API | https://api.clue.io/api/sigs | **401** `"User Key must be specified in the request header"` | free account → user key | FREE REGISTRATION |
| data.clue.io S3 | https://s3.amazonaws.com/data.clue.io/ | **403** | yes | BLOCKED without account |
| HMS LINCS DB API | https://lincs.hms.harvard.edu/db/api/ | **500** (server error) | — | BROKEN |
| **L1000 phase-1 data on GEO** | `ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/` | **200**, full listing incl. Level2–Level5 GCTX + all metadata (`sig_info`, `pert_info`, `cell_info`, SHA512SUMS) | none | **READY TODAY** |

Key finding: **GSE92742 (Broad LINCS L1000 Phase 1) is fully open via GEO** — no registration:
- Level5 COMPZ.MODZ signatures: `GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz` — HEAD **200**, 21.33 GB; Range GET returned **206** with valid gzip magic bytes.
- Metadata: `sig_info.txt.gz` (**200**, ~11 MB), `pert_info.txt.gz` (**200**, ~1.1 MB).
This means signature-reversal can be run locally (cosine/Kendall correlation of ALS signature vs 473k L1000 signatures) with zero accounts.

## 3. Alternative signature-reversal resources

| Resource | Test | Status | Verdict |
|---|---|---|---|
| L1000FWD web | https://maayanlab.cloud/l1000fwd/ | **200** (old l1000fwd.org does not resolve — DNS failure; use maayanlab.cloud) | UI READY |
| L1000FWD drug lookup API | `GET https://maayanlab.cloud/l1000fwd/search_all/riluzole` | **200**, JSON `[{"Phase":"Launched","MOA":"glutamate inhibitor","name":"RILUZOLE","id":"BRD-K21283037",...}]` | READY TODAY |
| L1000FWD reversal upload | `POST .../simulate_query` | **405** (POST rejected server-wide); all POST endpoints return 405 | BLOCKED programmatically right now |
| Enrichr addList | `POST https://maayanlab.cloud/Enrichr/addList` (multipart) | **200**, userListId 136322214 | READY TODAY |
| Enrichr enrich | `GET https://maayanlab.cloud/Enrichr/enrich?userListId=...&backgroundType=KEGG_2021_Human` | **200**, top hit "Protein processing in endoplasmic reticulum" adj-p 5.6e-4 | READY TODAY |
| SigCom LINCS | https://maayanlab.cloud/sigcom-lincs/ **200**; `data-api/api/v1/` root **200** but enrich returned empty body; exact contract undocumented in landing page | PARTIAL — usable via UI; API needs doc reading | FREE REGISTRATION-free but API unverified |
| signatureSearch (Bioconductor) | Data ships via ExperimentHub; host reachable; requires R/Bioconductor install | Not tested end-to-end here (R runtime) | READY TODAY (with R) |

## 4. Compound annotation APIs

| Source | Test | Status | Verdict |
|---|---|---|---|
| PubChem REST | `GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/riluzole/property/MolecularFormula,MolecularWeight,CanonicalSMILES/JSON` | **200** — CID 5070, C8H5F3N2OS, 234.20 g/mol | READY TODAY, no auth |
| ChEMBL web services | `GET https://www.ebi.ac.uk/chembl/api/data/molecule.json?pref_name__iexact=riluzole` | **200** — CHEMBL744, RILUZOLE, `max_phase: 4.0` (approved), first_approval 1995 | READY TODAY, no auth |
| DrugBank | `GET https://go.drugbank.com/drugs/DB00703` | **403** (bot-blocked; full data requires license) | BLOCKED for programmatic use |

## Recommended pilot stack

- **Disease signatures (ALS):**
  1. GSE124439 — human post-mortem CNS bulk RNA-seq, ALS (n=145) vs control (n=17). Primary human-tissue signature.
  2. GSE255602 — sporadic ALS iPSC-derived spinal-cord-chip MNs (n=71), processed salmon counts downloadable directly.
  3. GSE261875 — TDP-43 iPSC motor-neuron perturbation (n=32) as orthogonal validation context.
- **Signature reversal:** Download GSE92742 Level5 GCTX (21 GB, GEO, no auth) + sig/pert metadata and rank by correlation reversal locally. Use Enrichr (working API) for pathway sanity checks of the ALS query signature. L1000FWD search_all usable for drug lookup; its reversal-upload endpoint is currently broken.
- **Compound annotation:** ChEMBL API for max_phase/approval year (no auth) + PubChem REST for structure/properties (no auth). Skip DrugBank (403/licensed).
