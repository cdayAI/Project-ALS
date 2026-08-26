# AGENTS.md — Rules for AI contributors (Claude Code, Codex, Prime Agent, sub-agents)

This repo is a multi-agent research factory. Humans review; agents do the work.
Follow these rules exactly — they exist because agents already broke them once.

## Hard rules
1. **Never commit raw data.** `data/` is gitignored. Datasets are downloaded by
   pipelines into `data/` and stay local. If `git status` shows a file >50 MB,
   STOP — something is wrong.
2. **Work on branches**, not main: `exp/<name>` for experiments,
   `pipeline/<name>` for shared pipelines. Merge to main only via PR after an
   adversarial review exists in `reviews/`.
3. **Claim tasks on the board below** by editing this file in your branch
   (`- [ ]` -> `- [~] your-name`). Never work on a task someone else claimed.
4. **Experiments follow `experiments/_template/EXPERIMENT.md`**: config.yaml +
   run.py + outputs/ + VERDICT.md. VERDICT stays PENDING until a reviewer agent
   writes `reviews/<experiment>.md`. No self-certification.
5. **Positive controls are mandatory.** A pipeline that cannot rediscover known
   truth (documented per-experiment) produces no verdict at all.
6. **Cite or it didn't happen**: research claims need PMIDs/URLs; code needs the
   dataset IDs and parameter values that reproduce it.

- [~] **Data QC pipeline** `pipeline/data-qc` `pipelines/data_qc/`
      te_contamination_check.py - reusable pre-flight QC detecting non-gene
      feature classes (TE/repeat rows etc.) vs a gene-reference whitelist,
      quantifying their library-size share, emitting GO/CAUTION/NO-GO.
      Worked example: GSE124439 (TE rows = 25.5% of raw reads -> NO-GO).
      CLAIMED BY: sprint1-repurposing agent. To become mandatory pre-flight
      step for all count-matrix experiments.

## Context to read first (in order)
1. `README.md` - mission and layout
2. `docs_plan.md` - factory design
3. `research/00_synthesis.md` - what we're doing and why (10 min)
4. `research/04_pilot_access_report.md` - verified data access

## Task board
- [x] Research briefs 01-04 (Prime Agent fleet) 
- [x] Cross-brief synthesis (Prime Agent)
- [~] exp001 repurposing screen (sprint1-repurposing agent - IN PROGRESS, owns experiments/exp001_repurposing/)
- [~] **Biomarker/enrichment factory (H-008/H-009)** `pipelines/biomarker_enrichment/`
      PRO-ACT simulation framework + real-data application draft.
      CLAIMED BY: biomarker-factory (Prime Agent fleet)
      STATUS: simulation framework runs end-to-end on SIMULATED data (controls
      pass); data/proact/APPLICATION_DRAFT.md ready for human sign-off;
      experiments/exp002_nfl_enrichment/EXPERIMENT.md skeleton pre-registered.
      Real-data validation blocked on PRO-ACT access.
- [~] **C9orf72 DPR module (H-007)** `exp/h007-c9orf72-module` + data/gse303931/
      Nucleolar-stress/speckle module from isogenic iPSC data.
      CLAIMED BY: c9orf72-factory (Prime Agent fleet)
- [~] **Structure readiness (Stream D prep)** `pipelines/structure_readiness/`
      Folded-target shortlist + Boltz-2/DiffDock tooling. IDPs out of scope.
      CLAIMED BY: structure-factory (Prime Agent fleet)
- [x] **Stream C: trial matcher** `tools/trial_matcher/` - ClinicalTrials.gov API v2 wrapper:
      input = country/region, mutation status, disease stage; output = ranked recruiting ALS trials.
      Working query pattern already proven (see Prime Agent session notes / research briefs).
      CLAIMED BY: claude (pipeline/trial-matcher branch) - DONE, independently reviewed
      (reviews/trial_matcher.md). Real matching bugs found in both self-testing and
      independent review; see README's Known limitations before extending this.
- [ ] **Stream A scaffolding: causal target triage** `pipelines/causal_targets/` -
      fetch GWAS Catalog EFO_0000253 loci, druggable-genome overlay, colocalization
      plan for Project MinE summary stats (application pending). CLAIMED BY: (available)
- [ ] **exp001 adversarial review** - BLOCKED until exp001 completes. Reviewer must
      check: multiple testing, batch effects (post-mortem tissue covariates!), 
      positive control outcome, replication vs GSE255602/GSE261875.

## Hypothesis ledger (initial batch)
- [x] Draft hypotheses/H-001..H-010 (hypothesis-drafter agent, branch `hypotheses/initial-batch`).
      Coverage: 2 signature-reversal refinements of exp001 classes (H-001, H-002),
      2 causal-target GWAS-x-druggable-genome (H-003 TBK1 axis, H-004 UNC13A/sQTL),
      1 TDP-43 PTM modifier w/ perturbation-signature test path (H-005),
      1 mitochondrial/axonal-transport node (H-006), 1 C9orf72 DPR countermeasure
      w/ signature correlate (H-007), NfL-dynamics trial enrichment (H-008),
      subtype-to-stratum enrichment bridge (H-009), reversal-context robustness
      calibration for exp001 (H-010). All status: open. Not merged to main;
      each awaits pipeline assignment and adversarial review per AGENTS.md rules.
- [~] H-001/H-002 blocked-on: exp001 outputs (sprint1-repurposing agent) — refinement
      queries are pre-registered against config.yaml positive controls.
- [~] H-003/H-004 depend on: Stream A scaffolding `pipelines/causal_targets/`
      (unclaimed; MinE DAC application should be started early per research/04).

## Working-tree isolation (MANDATORY - incident 2026-08-25)
Never share one checkout between agents. On task start run:
  git worktree add ~/ALS-worktrees/<branch-slug> <branch>
and work ONLY inside that directory. Merge to main via a separate temp
worktree or PR. If you find another agent's uncommitted changes in your tree,
do NOT checkout/reset over them - create your own worktree and report it.
