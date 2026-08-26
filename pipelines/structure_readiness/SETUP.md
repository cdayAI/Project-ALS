# SETUP — Stream D structure stack (boltz / DiffDock)

*Verified on macOS arm64 (Apple Silicon), 16 GB RAM, ~143 GB free disk. Prepared by structure-factory.*

## TL;DR verdict
- **boltz (Boltz-2): installable and installed** in this directory's `venv/` (Python 3.11). CLI verified.
- **Boltz-2 weights: accessible** (HTTP 200, Hugging Face + model-gateway mirrors). Minimal viable set = `boltz2_conf.ckpt` (2.29 GB) + `mols.tar` (1.86 GB) ≈ **4.15 GB** downloads.
- Full stack including affinity head (`boltz2_aff.ckpt`, 2.06 GB), legacy Boltz-1 (`boltz1_conf.ckpt`, 3.60 GB), and DiffDock's ESM-2 weights (~2.5 GB) would total **>10 GB → deliberately NOT downloaded** per task rule (>10 GB stop threshold).
- **DiffDock: NOT installed.** Its pinned requirements conflict with modern Python and it needs ESM-2 weights; documented below for a future GPU machine.

## Exact commands used

```bash
cd pipelines/structure_readiness
uv venv venv --python 3.11          # system default is 3.14; boltz needs <=3.11-era stack
uv pip install --python venv/bin/python boltz    # resolves cleanly on macOS arm64, torch 2.13.0 arm64 wheel ~111 MB
./venv/bin/python -c "import boltz"              # OK
./venv/bin/boltz --help                          # OK

# weights (default cache ~/.boltz, override with BOLTZ_CACHE)
curl -L -C - -o ~/.boltz/mols.tar          https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar
curl -L -C - -o ~/.boltz/boltz2_conf.ckpt  https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt
```

## Smoke test
`smoke_test/tbk1_kd_1FV.yaml`: TBK1 kinase domain (Q9UHD2 residues 1–310, UniProt domain annotation "protein kinase") + ligand **1FV**, the co-crystallized inhibitor of PDB **4IM0** (human TBK1 KD–inhibitor complex). SMILES derived from the RCSB chemcomp InChI via RDKit.

Run:
```bash
./venv/bin/boltz predict smoke_test/tbk1_kd_1FV.yaml \
    --out_dir smoke_test/outputs \
    --accelerator cpu \
    --recycling_steps 1 \
    --sampling_steps 50 \
    --diffusion_samples 1 \
    --write_full_pae false
```
Gotchas found during bring-up:
- boltz refuses inputs without MSAs unless the protein block has `msa: empty` (single-sequence mode) or `--use_msa_server` is set. Our smoke input uses `msa: empty`.
- Even for structure-only runs boltz pre-fetches the affinity checkpoint (`boltz2_aff.ckpt`, ~2.06 GB) into the cache on first invocation.
- On Apple Silicon, Lightning sees MPS but boltz runs were launched with `--accelerator cpu`; MPS support is untested here.

Reduced recycling/sampling steps keep a CPU-only run tractable on 16 GB RAM; raise to defaults (`--recycling_steps 3 --sampling_steps 200`) once any CUDA/MPS-capable box is available.
Positive-control criterion for later full runs: predicted ligand pocket should coincide with the ATP-site hinge region seen in 4IM0 chain A (compare contact residues).

## Disk & RAM budget (16 GB-RAM machine)
| Item | Download | On-disk |
|---|---|---|
| venv (torch 2.13 arm64 + deps) | ~600 MB | ~2 GB |
| boltz2_conf.ckpt | 2.29 GB | 2.29 GB |
| mols.tar | 1.86 GB | ~5 GB extracted (ccd + ram generated on first run) |
| **Minimal working set** | **~4.8 GB** | **~9-10 GB** |
| Deferred: boltz2_aff.ckpt / boltz1_conf.ckpt / DiffDock+ESM | ~8 GB more | not downloaded |

RAM: Boltz-2 inference on one kinase-domain-sized complex fits in 16 GB with reduced sampling; full-length multi-chain targets may swap — chunk targets to domains (this shortlist's whole purpose).

## DiffDock (deferred, install recipe)
```bash
# DiffDock requires python <=3.10 in practice; separate venv:
uv venv diffdock_venv --python 3.10 && source diffdock_venv/bin/activate
git clone https://github.com/gcorso/DiffDock ../third_party/DiffDock   # path TBD if adopted
pip install -r ../third_party/DiffDock/requirements.txt
# weights: https://github.com/gcorso/DiffDock/releases/download/v1.0/original_weights.tar (~few hundred MB)
# plus facebook/esm esm2_t33_650M_UR50D sequence embeddings (~2.5 GB) on first run
```
Not attempted here: >10 GB cumulative threshold reached, CPU-only DiffDock inference is impractical, and briefs warn docking against raw AF models has high false-positive rates ([PMID 37546760](https://pubmed.ncbi.nlm.nih.gov/37546760/)) — prefer Boltz-2 co-folding as primary screen.

## Sources
- Boltz-2: Passaro, Corso, Wohlwend et al., bioRxiv 2025 (EPMC PPR1039145); code https://github.com/jwohlwend/boltz
- AlphaFold DB entries verified live (model_v6 URLs, HTTP 200) per accession in `target_shortlist.md`.
- PDBe best_structures API used for experimental-structure counts.
