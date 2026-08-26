#!/usr/bin/env python
"""GSE283507 replication arm for exp002/H-007b - count-level reprocessing.

Remediation of reviews/exp002_review.md sec 3:
  - RAW COUNTS (_Read_Count), not the artifact TPM matrix
  - FULL gene universe background (every gene in the deposited quantification matrix)
  - per-family GO enrichment (no union collapse), direction-specific (down-regulated)
  - this file IS the committed analysis code (AGENTS.md rule 6)

Usage: <python> replication_gse283507.py [--with-fcb]
Writes results to results/gse283507_countlevel_enrichment.csv and
results/gse283507_countlevel_de.csv, and (with --with-fcb) FC-B rescue tables.
"""
import argparse, os, sys
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.abspath(os.path.dirname(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import run as exp002  # noqa: E402 (load_go_sets + DPR_CURATED)

DATA = os.path.join(ROOT, "data", "gse283507", "GSE283507_raw_Count_FPKM_TPM.csv.gz")
RESULTS = os.path.join(HERE, "results")
SEED = 20260825

# families reported separately; only PRIMARY gates H-007b FC-A
FAMILIES = {
    "primary_cytoskeleton_organization": ("GO:0007010",),
    "secondary_ECM": ("GO:0031012",),
    "secondary_synapse": ("GO:0045202",),
    "secondary_cell_adhesion": ("GO:0007155",),
    "killed_nucleolar_ribo": ("GO:0005730", "GO:0042254"),
    "killed_speckle_splicing": ("GO:0016607", "GO:0000398"),
}

def median_ratio_size_factors(counts):
    """DESeq2-style size factors on a genes x samples raw-count matrix."""
    valid = counts > 0
    geom = np.exp(np.nanmean(np.log(counts.astype(float)), axis=1))
    ok = geom > 0
    ratios = counts[ok] / geom[ok, None]
    return np.median(ratios, axis=0)

def load_counts():
    df = pd.read_csv(DATA)
    cnt_cols = [c for c in df.columns if c.endswith("_Read_Count")]
    mat = df[["Gene_Symbol"] + cnt_cols].dropna(subset=["Gene_Symbol"])
    mat = mat.groupby("Gene_Symbol")[cnt_cols].sum()          # collapse duplicate symbols
    return mat

def dmso_de(mat):
    cols = [c for c in mat.columns if "DMSO" in c and (c.startswith("RC802") or c.startswith("TDP43"))]
    sub = mat[cols]
    meta = pd.DataFrame({"sample": cols})
    meta["genotype"] = ["ctrl" if s.startswith("RC802") else "mutant" for s in cols]
    meta["time"] = meta["sample"].str.extract(r"(6H|36H|84H)")
    # detection filter for TESTING (documented); universe below uses all genes
    keep = (sub >= 10).sum(axis=1) >= 9                       # >=10 raw counts in half the samples
    tested = sub[keep]
    sf = median_ratio_size_factors(tested.values)
    norm_log = np.log2(tested.div(sf, axis=1) + 1)
    T = pd.get_dummies(meta["time"], drop_first=True).astype(float)
    X = np.column_stack([np.ones(len(meta)), (meta["genotype"] == "mutant").astype(float), T.values])
    n, p_ = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    Y = norm_log.values.T
    resid = Y - X @ (XtXi @ X.T @ Y)
    s2 = (resid ** 2).sum(axis=0) / (n - p_)
    beta = XtXi @ X.T @ Y
    tc = beta[1] / np.sqrt(s2 * XtXi[1, 1])
    pv = 2 * stats.t.sf(np.abs(tc), n - p_)
    res = pd.DataFrame({"symbol": norm_log.index, "log2FC_mut_vs_ctrl": beta[1],
                        "t": tc, "p": pv})
    res["fdr_bh"] = stats.false_discovery_control(res["p"].values) \
        if hasattr(stats, "false_discovery_control") else _bh(res["p"].values)
    return res.sort_values("t", ascending=False).reset_index(drop=True)

def _bh(p):
    p = np.asarray(p); n = len(p)
    o = np.argsort(p); ranked = p[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(ranked, 0, 1)
    return out

def enrich_permutation(universe_syms, module_syms, set_syms, rng, nperm=20000):
    u = np.array(universe_syms)
    S = np.array([s in set(set_syms) for s in u], dtype=bool)
    M = np.array([s in set(module_syms) for s in u], dtype=bool)
    k, m = int(M.sum()), int(S.sum())
    if k == 0 or m == 0:
        return None
    obs = int(S[M].sum())
    exp = m * k / len(u)
    hits = np.flatnonzero(S)
    ge = 0
    folds_null = np.empty(nperm)
    for i in range(nperm):
        draw = rng.choice(len(u), size=k, replace=False)
        c = int(np.isin(draw, hits).sum())
        ge += c >= obs
        folds_null[i] = c / max(exp, 1e-9)
    return {"module_k": k, "set_m_in_universe": m, "observed_overlap": obs,
            "expected_overlap": round(exp, 3),
            "fold_enrichment": round(obs / max(exp, 1e-9), 3),
            "perm_p": round((ge + 1) / (nperm + 1), 5),
            "null_fold_mean": round(float(folds_null.mean()), 3),
            "null_fold_p95": round(float(np.percentile(folds_null, 95)), 3)}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-fcb", action="store_true")
    args = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)

    mat = load_counts()
    universe = list(mat.index)                                 # FULL gene universe of the matrix
    print(f"full gene universe: {len(universe)} genes "
          f"(deposited quantification matrix limit)")

    res = dmso_de(mat)
    res.to_csv(os.path.join(RESULTS, "gse283507_countlevel_de.csv"), index=False)
    down = list(res[(res.fdr_bh < 0.10) & (res.log2FC_mut_vs_ctrl < 0)].symbol)
    up = list(res[(res.fdr_bh < 0.10) & (res.log2FC_mut_vs_ctrl > 0)].symbol)
    # variance sanity check vs artifact TPMs
    sf_all = median_ratio_size_factors(mat.values)
    nl = np.log2(mat.div(sf_all, axis=1) + 1)
    grp_var = np.stack([nl.iloc[:, 3*i:3*i+3].var(axis=1) for i in range(6)]).mean()
    print(f"median within-group log2 variance (counts): {np.nanmedian(grp_var):.3f} "
          f"(artifact TPMs were ~0.019)")
    print(f"DE: {len(res)} tested | FDR<0.10 down={len(down)} up={len(up)}")

    go = exp002.load_go_sets(os.path.join(ROOT, "data", "gse303931", "ref"))
    fam_sets = {}
    for fam, terms in FAMILIES.items():
        s = set().union(*[go[t] for t in terms]) if len(terms) > 1 else set(go[terms[0]])
        fam_sets[fam] = s
    fam_sets["DPR_curated_literature"] = set(exp002.DPR_CURATED)

    rng = np.random.default_rng(SEED)
    rows = []
    modules = {"down_FDR10": down, "up_FDR10": up,
               "top500_absT": list(pd.concat([res.nlargest(500, "t"), res.nsmallest(500, "t")]).symbol)}
    for mname, msyms in modules.items():
        for fam, ss in fam_sets.items():
            r = enrich_permutation(universe, msyms, ss, rng)
            if r:
                rows.append({"module": mname, "family": fam,
                             "gating_family": fam.startswith("primary"), **r})
    enr = pd.DataFrame(rows)
    enr.to_csv(os.path.join(RESULTS, "gse283507_countlevel_enrichment.csv"), index=False)
    gate = enr[(enr.module == "down_FDR10") & enr.gating_family]
    if len(gate):
        g = gate.iloc[0]
        fc_a = bool(g.fold_enrichment >= 2.0 and g.perm_p <= 0.05)
        print(f"[FC-A primary] fold={g.fold_enrichment} P={g.perm_p} -> {'PASS' if fc_a else 'FAIL'}")
    else:
        print("[FC-A primary] no module genes overlap family - FAIL")

    if args.with_fcb:
        fcb_rows = []
        genotypes = {"RC802": "WT_ctrl", "TDP43": "TDP43M337V", "TD-KO": "TDP43M337V_DRD2KO"}
        mod_genes = [g for g in fam_sets["primary_cytoskeleton_organization"]
                     if g in mat.index]
        print(f"FC-B module genes present in matrix: {len(mod_genes)}")
        for gcol, gname in genotypes.items():
            for t in ["6H", "36H", "84H"]:
                d = [f"{gcol}DMSO{t}-{r}_Read_Count" for r in (1, 2, 3)]
                q = [f"{gcol}ROPI{t}-{r}_Read_Count" for r in (1, 2, 3)]
                sub = mat.loc[mod_genes, d + q]
                sf = median_ratio_size_factors(sub.values)
                l = np.log2(sub.div(sf, axis=1) + 1)
                lfc = l[q].mean(axis=1) - l[d].mean(axis=1)
                tt = stats.ttest_1samp(lfc.values, 0.0)
                fcb_rows.append({"genotype": gname, "time": t,
                                 "mean_module_lfc_ropi_dmso": float(lfc.mean()),
                                 "p": float(tt.pvalue)})
        tab = pd.DataFrame(fcb_rows)
        wt = tab[tab.genotype == "WT_ctrl"].mean_module_lfc_ropi_dmso.mean()
        resc = []
        for gn in ["TDP43M337V", "TDP43M337V_DRD2KO"]:
            m = tab[tab.genotype == gn].mean_module_lfc_ropi_dmso.mean()
            resc.append({"genotype": gn, "rescue_index_R": float(m - wt),
                         "rescue_positive": bool((m - wt) > 0)})
        resc = pd.DataFrame(resc)
        tab.to_csv(os.path.join(RESULTS, "fcB_countlevel_per_time.csv"), index=False)
        resc.to_csv(os.path.join(RESULTS, "fcB_countlevel_rescue_index.csv"), index=False)
        n_pos = int(resc.rescue_positive.sum())
        print("[FC-B]", resc.to_string(index=False))
        print(f"[FC-B] rescue-positive: {n_pos}/2 disease genotypes -> "
              f"{'PASS' if n_pos == 2 else 'FAIL'}")

if __name__ == "__main__":
    main()
