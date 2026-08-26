# Data QC: non-gene feature contamination check

- Generated: 2026-08-25 21:29:39
- Input: `/Users/christopherday/Project-ALS/data/gse124439`
- Features analysed: 28,953 x 176 samples
- Gene reference whitelist: **YES**

## Verdict: **NO-GO**

Non-gene features carry 25.7% of raw reads. Strip them BEFORE CPM/TPM/size-factor normalization; any analysis normalized on the raw matrix is globally diluted and should be re-run.

| class | n features | % features | % raw reads | per-sample read % (min/med/max) |
|---|---|---|---|---|
| gene | 26,994 | 93.23% | 74.29% | 0.18% / 0.40% / 0.89% |
| other_non_gene | 967 | 3.34% | 0.20% | 0.00% / 0.00% / 0.00% |
| transposable_element_repeat | 992 | 3.43% | 25.51% | 0.06% / 0.14% / 0.35% |

## Top non-gene features by raw reads

| feature | class | reads |
|---|---|---|
| MIRb:MIR:SINE | transposable_element_repeat | 64,491,011 |
| L2a:L2:LINE | transposable_element_repeat | 49,891,830 |
| MIR:MIR:SINE | transposable_element_repeat | 48,567,195 |
| L2c:L2:LINE | transposable_element_repeat | 42,458,567 |
| AluJb:Alu:SINE | transposable_element_repeat | 29,096,772 |
| MIRc:MIR:SINE | transposable_element_repeat | 28,496,992 |
| L2b:L2:LINE | transposable_element_repeat | 25,401,875 |
| MIR3:MIR:SINE | transposable_element_repeat | 23,956,408 |
| AluSx:Alu:SINE | transposable_element_repeat | 23,052,178 |
| 7SK:RNA:RNA | transposable_element_repeat | 21,172,145 |
| LOC100996724 | other_non_gene | 916,453 |
| MIR3687-2 | other_non_gene | 861,526 |
| LOC728392 | other_non_gene | 483,974 |
| LOC100506548 | other_non_gene | 293,870 |
| MIR3687-1 | other_non_gene | 279,384 |
| MIR6087 | other_non_gene | 198,179 |
| LOC642852 | other_non_gene | 172,367 |
| LOC100419583 | other_non_gene | 157,163 |
| LOC646762 | other_non_gene | 131,744 |
| LOC220729 | other_non_gene | 110,843 |
