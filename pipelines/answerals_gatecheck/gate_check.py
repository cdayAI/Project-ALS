#!/usr/bin/env python
"""Apply the PRE-REGISTERED H-011/H-012 go/no-go gates to an AnswerALS sample-metadata CSV.

Gates (from hypotheses/H-011.md and H-012.md power_analysis blocks - committed before any
data access):
  G1 (both): post-QC n >= 15 C9 carriers AND >= 15 controls with qc-passing RNA-seq
             (if 10<=n<15: DOWNGRADE flag for H-011 primary endpoint; if n<10: FAIL)
  G2 (H-012 pilot precondition): metadata must identify enough donors to estimate PSI
             variance in stage 1; requires >= 60 total post-QC samples as a floor
  G3 (design sanity): C9 carriers must have repeat-size annotations (>26) for >=90%;
             sex must be annotated for >=95% (sex-concordance QC is L7)

Refuses files whose sample_ids/donor_ids carry the 'MOCK' marker unless --allow-mock.
Exit code: 0 if all applicable gates PASS, 2 otherwise.
"""
import argparse, json, os
import pandas as pd

REQUIRED_COLS = ["donor_id","sample_id","genotype_group","sex","c9_repeat_size",
                 "culture_batch","rna_qc_pass"]

def check(meta_path, allow_mock=False):
    df = pd.read_csv(meta_path)
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        return {"verdict":"FAIL","reason":f"missing required columns: {missing}"}
    if not allow_mock:
        blob = df["donor_id"].astype(str) + df["sample_id"].astype(str)
        if blob.str.contains("MOCK", case=False).any():
            return {"verdict":"FAIL","reason":"MOCK data detected; rerun with --allow-mock for code-path validation"}
    q = df[df.rna_qc_pass.astype(str).str.lower().isin(["true","1","yes"])]
    res = {"total_rows": len(df), "post_qc_rows": len(q)}
    for g in ["C9_carrier","nonC9_ALS","control"]:
        res[f"post_qc_{g}"] = int((q.genotype_group == g).sum())
    n_c9, n_ctrl = res["post_qc_C9_carrier"], res["post_qc_control"]

    g1 = {"name":"G1_min_n","required":"C9>=15 AND control>=15"}
    if n_c9 >= 15 and n_ctrl >= 15:
        g1["status"] = "PASS"
    elif 10 <= min(n_c9, n_ctrl) < 15:
        g1["status"] = "DOWNGRADE"   # H-011 primary endpoint drops to module-level aggregate
    else:
        g1["status"] = "FAIL"
    res["gates"] = [g1]

    g2 = {"name":"G2_pilot_floor_total_samples","required":"post-QC total >= 60"}
    g2["status"] = "PASS" if len(q) >= 60 else "FAIL"
    res["gates"].append(g2)

    carr = q[q.genotype_group == "C9_carrier"]
    rep_annot = carr.c9_repeat_size.notna().mean() if len(carr) else float("nan")
    sex_annot = q.sex.notna().mean()
    g3 = {"name":"G3_design_sanity",
          "repeat_annotation_rate_c9": round(float(rep_annot),3),
          "sex_annotation_rate": round(float(sex_annot),3)}
    g3["status"] = "PASS" if (rep_annot >= 0.9 and sex_annot >= 0.95) else "FAIL"
    res["gates"].append(g3)

    hard_fail = any(g["status"] == "FAIL" for g in res["gates"])
    res["verdict"] = "FAIL" if hard_fail else ("DOWNGRADE" if any(g["status"]=="DOWNGRADE" for g in res["gates"]) else "PASS")
    return res

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--allow-mock", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    res = check(a.metadata, a.allow_mock)
    print(json.dumps(res, indent=1))
    if a.out:
        json.dump(res, open(a.out,"w"), indent=1)
    raise SystemExit(0 if res.get("verdict") == "PASS" else 2)

if __name__ == "__main__":
    main()
