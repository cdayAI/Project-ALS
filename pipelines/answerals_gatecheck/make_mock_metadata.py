#!/usr/bin/env python
"""Generate MOCK AnswerALS iPS-neuron sample metadata shaped per PMID 35115730 /
PMC8825283 supplementary descriptions. For CODE-PATH VALIDATION ONLY - never for
hypothesis testing (hypotheses/_LESSONS.md L4).

Shape derived from published facts:
- population: >1,000 participants; ~850+ iPS lines banked (1 clone/donor)
- genotype strata present: sporadic ALS, C9orf72 carriers (>26 repeats), other fALS, controls
- 41 ALS + 4 pre-fALS C9 carriers among 830 WGS'd participants
- spinal neurons from 32-day diMN protocol; multi-omics incl. RNA-seq per donor
Mock defaults are deliberately SMALLER than expected reality so that a real pull which
passes gates cannot be confused with mock output.
"""
import argparse, os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-c9", type=int, default=6, help="mock C9 carriers (default deliberately small)")
    ap.add_argument("--n-nonc9", type=int, default=10)
    ap.add_argument("--n-control", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20260826)
    ap.add_argument("--outdir", default=os.path.join(HERE, "outputs"))
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(a.outdir, exist_ok=True)
    rows = []
    def add(group, n):
        for i in range(n):
            carrier = group == "C9_carrier"
            rows.append({
                "donor_id": f"MOCK-{group[:3].upper()}-{i+1:03d}",
                "sample_id": f"MOCK-RNA-{group[:3].upper()}-{i+1:03d}",
                "genotype_group": group,
                "sex": rng.choice(["M","F"], p=[0.6,0.4]),
                "c9_repeat_size": int(rng.integers(30, 800)) if carrier else np.nan,
                "culture_batch": f"B{rng.integers(1,5)}",
                "rna_qc_pass": bool(rng.random() > 0.05),
            })
    add("C9_carrier", a.n_c9)
    add("nonC9_ALS", a.n_nonc9)
    add("control", a.n_control)
    df = pd.DataFrame(rows)
    out = os.path.join(a.outdir, "mock_sample_metadata.csv")
    df.to_csv(out, index=False)
    print(f"mock metadata written: {out} ({len(df)} rows; prefix 'MOCK-' marks synthetic)")

if __name__ == "__main__":
    main()
