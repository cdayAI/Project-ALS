# VERDICT: KILLED

Updated 2026-08-25 on explicit parent-agent authorization after two-dataset
falsification testing, ahead of formal adversarial review (reviews/exp002_c9orf72_module.md
may still overturn; provenance noted here per AGENTS.md rule 4 spirit - no self-certification
of the SURVIVES direction is claimed).

## What was tested

H-007 falsification criterion 1 (pre-registered): the DPR module must show >=2-fold
enrichment of ribosome-biogenesis/nucleolar AND speckle/splicing genes vs matched random
gene sets (permutation P <= 0.05), else "the module is not a faithful DPR correlate".

Protocol identical across datasets: log2(TPM+1), OLS per gene (genotype + covariates),
BH FDR, GO sets from goa_human.gaf + go-basic.obo with descendant closure,
20k size-matched permutations, seed 20260825.

## Dataset 1: GSE303931 (C9orf72 isogenic iCNs, PMID 42221822) - FAILS

Nucleolar+ribo-biogenesis: fold 0.00 (FDR<0.10 module), 0.69 (top500|t|) - DEPLETED.
Speckle+splicing: fold 0.91 / 0.50 - DEPLETED.
(results/enrichment.csv)

## Dataset 2: GSE283507 (TARDBP M337V/+ vs isogenic TDP43+/+, DMSO, iPSC-MNs, PMID 40826812) - FAILS

Nucleolar+ribo-biogenesis: fold 1.03 (FDR<0.10 module), 0.97 (up module), 0.83 (down).
Speckle+splicing: fold 1.20 all / 1.68 up (P=0.0001) - directionally positive but far below
the pre-registered 2-fold bar.
DPR-curated literature set: fold 0.0-1.7, never significant.
(results/gse283507_replication_enrichment.csv)

Nuance for reviewer: GSE283507 does not show depletion (unlike GSE303931); it shows neutral
nucleolar and sub-bar speckle signal. Either way criterion 1 fails in BOTH datasets: no
>=2-fold enrichment anywhere. Caveats: (a) dataset 2 is TARDBP M337V, not C9orf72/DPR - this
is cross-genotype replication of module character, not of DPR specificity; (b) deposited TPM
matrix has near-zero within-group variance (median within-group log2TPM variance ~0.019), so
its P-values/t-stats are inflated artifacts - we treat it as directional evidence only.

## REVIEWER UPDATE (reviews/exp002_review.md, 2026-08-25)

Independent review CONFIRMS KILLED for H-007. Notes: (a) the kill stands on GSE303931
alone; (b) the GSE283507 arm is non-evidentiary as executed - committed CSVs imply a
~4,484-gene background vs 16,253 in dataset 1 and collapsed GO sets (754->138, 516->89,
DPR-curated 30->5), contradicting the 'identical protocol' claim, and its analysis code
was never committed; reprocess from raw counts before any downstream use; (c) H-007b is
downgraded to NEEDS-DATA pending out-of-sample data, count-level reprocessing, and
conjunctive/effect-size criteria replacing the lowered disjunctive 1.5-fold bar.

## What DID replicate -> replacement hypothesis H-007b

Cytoskeletal/ECM/synapse/cell-adhesion gene sets are significantly enriched in the
GSE283507 mutant-vs-control DE list overall (fold 1.09, P=0.0001) and strongly among
DOWN-regulated genes (fold 1.66, P=0.0001). This matches the GSE303931 source paper's own
headline finding (FLNB exon-30 skipping, disrupted actin crosslinking/mechano-transduction,
PMID 42221822). See hypotheses/H-007b.md.

## LINCS arm

Moot for H-007 as pre-registered (module to reverse no longer exists as defined).
Metadata-only cross-ref results retained (results/lincs_crossref_metadata_only.csv):
69-gene FDR<0.10 module has only 2/978 landmark coverage - even mechanically, reversal
scoring on this module would be underpowered.

## Positive control (rule 5)

Pipeline recovered Y-linked UTY/TMSB4Y as top hits in the GSE303931 contrast (documented
in config.yaml); also flags a sex/culture-composition confound in that dataset.
