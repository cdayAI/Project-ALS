# ALS Research Factory — Operating Plan
*How AI agents maximize research throughput and test ideas at machine speed*

## Core principle
Close the loop: HYPOTHESIZE -> EXPERIMENT -> ATTACK -> TRIAGE -> repeat.
Human attention goes only to judgment calls. Everything else is automated.

## The loop (5 layers)

### Layer 0 — Shared knowledge base (single source of truth)
A git repo (`als-factory/`) holding everything agents read/write:
- `hypotheses/` — hypothesis ledger. Each hypothesis is a structured file:
  mechanism, evidence links (PMIDs/GSE IDs), proposed experiment,
  falsification criteria, cost (compute only?), status (open/killed/validated).
- `data/` — curated dataset inventory + downloaded subsets
- `experiments/` — one folder per experiment: code, config, outputs, verdict
- `reviews/` — adversarial agent reports on each result
- `digests/` — weekly human-facing summaries

### Layer 1 — Hypothesis generation (orchestrator agents like me)
- Mine literature (PubMed E-utilities API, free) + omics data for contradictions,
  gaps, and cross-domain connections.
- Score each hypothesis: novelty, testability IN SILICO, expected information gain
  if killed, if confirmed.
- Only hypotheses cheaply falsifiable with public data enter the queue first.

### Layer 2 — Execution fleet (Claude Code / Codex / autoresearch loops)
- Each approved hypothesis -> an agent implements a standardized, reproducible
  pipeline in the repo (differential expression, signature-reversal scoring,
  network propagation, docking, antisense design).
- Standard templates so any agent can run any experiment identically.
- Autoresearch-style sweeps run unattended over parameter grids / compound sets;
  results append to the same results schema every time.

### Layer 3 — Adversarial review (a second agent per result)
Every result gets attacked before it counts:
- statistics (multiple testing, effect size, power),
- batch effects and confounders in the underlying data,
- positive controls (does the method rediscover known ALS drugs?),
- replication in an independent dataset.
Verdicts: KILLED / SURVIVES / NEEDS-DATA-WE-DON'T-HAVE.

### Layer 4 — Human checkpoints & external interface
- Weekly digest: top-N surviving hypotheses + next cheapest tests. You decide.
- What agents cannot do: wet lab, patients, clinical trials. Interfaces:
  * iPSC motor-neuron assay screens are commercially available (CROs) -
    a validated in-silico shortlist becomes a purchasable experiment.
  * Open-data consortia (Answer ALS, NEALS) for collaboration/validation.
  * Preprints to attract domain collaborators who CAN run wet-lab follow-ups.

## Tool division of labor
| Tool | Role |
|---|---|
| Prime Agent (me) | orchestrator, literature mining, synthesis, delegation, digests |
| Claude Code / Codex | heavy pipeline implementation in the repo, long refactors |
| autoresearch loops | unattended batch sweeps: compound screens, parameter grids |
| Shared git repo | the factory floor - all state lives here, not in chat history |

## Pilot Sprint #1 (starts today): transcriptome drug-repurposing loop
Goal: prove the full loop end-to-end in days, produce a REAL artifact.
1. Pull ALS-vs-control RNA-seq signatures from public GEO datasets
   (spinal cord and/or iPSC motor neuron cohorts).
2. Compute robust disease signature; score LINCS L1000/CMap compounds by
   signature REVERSAL.
3. Positive control: known ALS-active drugs should rank sensibly.
4. Cross-check top candidates against safety/max-phase (DrugBank/ChEMBL).
5. Adversarial agent attacks the ranking; survivors go in the digest.
Deliverable: ranked repurposing shortlist + fully reproducible repo +
a template every future experiment reuses.

## Honest constraints
- No wet-lab capability in-loop; validation requires partners or purchased assays.
- Public-data conclusions are hypothesis-generating until replicated externally.
- Speed comes from parallelism + falsification discipline, not from skipping steps.
