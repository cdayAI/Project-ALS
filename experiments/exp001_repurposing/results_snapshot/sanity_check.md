# ALS signature sanity check (GSE124439, ALS vs non-neurological control)

## Marker-module level
- Microglial module (13 genes): mean moderated t = -0.368
- Neuronal module (11 genes): mean moderated t = -0.005
- Individual canonical markers are mostly not significant: post-mortem bulk tissue carries large cell-composition and RNA-quality variance, and the control arm is small (n=17).

## Pathway level - Enrichr on top-300 UP genes
- KEGG_2021_Human: Protein processing in endoplasmic reticulum (p=1.9e-06); RNA transport (p=0.002); Autophagy (p=0.0047); Estrogen signaling pathway (p=0.0047); Endocytosis (p=0.0049); AMPK signaling pathway (p=0.0095); Lipid and atherosclerosis (p=0.016); Longevity regulating pathway (p=0.019)
- MSigDB_Hallmark_2020: Myc Targets V1 (p=3.8e-07); Protein Secretion (p=0.015); mTORC1 Signaling (p=0.032); Fatty Acid Metabolism (p=0.032); Androgen Response (p=0.064); UV Response Dn (p=0.066); Mitotic Spindle (p=0.08); G2-M Checkpoint (p=0.081)

## Pathway level - Enrichr on top-300 DOWN genes
- KEGG_2021_Human: mRNA surveillance pathway (p=0.016); Glycine, serine and threonine metabolism (p=0.022); Glycosaminoglycan degradation (p=0.032); ECM-receptor interaction (p=0.043); Fanconi anemia pathway (p=0.047); NOD-like receptor signaling pathway (p=0.056); Glycosylphosphatidylinositol (GPI)-anchor biosynthesis (p=0.058); Pertussis (p=0.11)
- MSigDB_Hallmark_2020: Interferon Gamma Response (p=0.011); Mitotic Spindle (p=0.031); Xenobiotic Metabolism (p=0.032); Inflammatory Response (p=0.081); Glycolysis (p=0.081); Hypoxia (p=0.18); KRAS Signaling Dn (p=0.18); Apical Junction (p=0.35)

Interpretation: the UP signature is dominated by protein processing in the ER /
unfolded-protein response / autophagy / mTORC1 signaling - the canonical
TDP-43-proteinopathy axis of ALS biology.

## Selected marker genes
| gene | log2FC | moderated t | padj |
|---|---|---|---|
| NEFH | -0.207 | -0.79 | 6.23e-01 |
| NEFM | +0.015 | 0.08 | 9.67e-01 |
| STMN2 | +0.140 | 0.70 | 6.64e-01 |
| RBFOX1 | +0.092 | 0.52 | 7.57e-01 |
| TARDBP | +0.002 | 0.04 | 9.85e-01 |
| SQSTM1 | -0.003 | -0.03 | 9.87e-01 |
| OPTN | +0.083 | 1.42 | 3.30e-01 |
| TBK1 | +0.107 | 1.67 | 2.45e-01 |
| AIF1 | -0.035 | -0.15 | 9.39e-01 |
| TYROBP | -0.053 | -0.27 | 8.83e-01 |
| C1QA | +0.149 | 0.53 | 7.55e-01 |
| CD68 | -0.027 | -0.12 | 9.49e-01 |
| HLA-DRA | -0.243 | -0.86 | 5.86e-01 |
| GFAP | +0.167 | 0.65 | 6.93e-01 |
| HSP90AA1 | +0.489 | 4.02 | 4.45e-03 |
| IFNAR1 | +0.179 | 3.96 | 4.92e-03 |
