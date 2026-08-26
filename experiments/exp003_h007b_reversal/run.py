#!/usr/bin/env python
"""exp003_h007b_reversal - LINCS Level5 reversal scoring of the cytoskeletal/ECM/synaptic
iPSC ALS module (hypothesis H-007b).

Stages (all run by default, in order):
  fc_b   - ropinirole within-genotype module-response test (local GSE283507 TPMs)
  build  - build query vectors q1 (C9orf72, GSE303931) and q2 (TARDBP, GSE283507)
  score  - one streaming pass over data/lincs Level5 MODZ GCTX -> per-signature scores
  agg    - aggregate to perturbagen level, FC-D consistency PC, annotate top-50

Usage:
  <python> experiments/exp003_h007b_reversal/run.py [--only fc_b|build|score|agg]

Requires pandas/numpy/scipy/h5py/requests. Reuses DE code from exp002_c9orf72_module.
"""
import argparse, gzip, json, os, sys, time
from collections import deque
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.abspath(os.path.dirname(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "experiments", "exp002_c9orf72_module"))
sys.path.insert(0, os.path.join(ROOT, "pipelines", "perturbation_signatures"))
import run as exp002                      # noqa: E402  (DE + GO machinery)
from lincs_score import read_gctx_metadata, stream_scores_multi, align_query_to_gctx  # noqa

LINCS_DIR = os.environ.get("LINCS_DIR",
    os.path.expanduser("~/Project-ALS/data/lincs"))
GCTX = os.path.join(LINCS_DIR, "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx")
SEED = 20260825
OUT = os.path.join(HERE, "outputs")
RESULTS = os.path.join(HERE, "results")

MODULE_GO = {"cytoskeleton_organization": "GO:0007010", "extracellular_matrix": "GO:0031012",
             "synapse": "GO:0045202", "cell_adhesion": "GO:0007155"}

# ---------------------------------------------------------------- GO sets
def load_go_descendants(refdir, terms):
    ann, parents, alt = {}, {}, {}
    with gzip.open(os.path.join(refdir, "goa_human.gaf.gz"), "rt") as f:
        for line in f:
            if line.startswith("!"):
                continue
            p_ = line.rstrip("\n").split("\t")
            if len(p_) < 15 or p_[3].startswith("NOT"):
                continue
            ann.setdefault(p_[4], set()).add(p_[2])
    cur = None
    with open(os.path.join(refdir, "go-basic.obo")) as f:
        for line in f:
            line = line.rstrip()
            if line == "[Term]":
                cur = {"is_a": []}
            elif cur is None:
                continue
            elif line == "":
                if cur and "id" in cur:
                    parents[cur["id"]] = cur["is_a"]
                    for a in cur.get("alt_ids", []):
                        alt[a] = cur["id"]
                cur = None
            elif cur is not None:
                if line.startswith("id: "): cur["id"] = line[4:]
                elif line.startswith("alt_id: "): cur.setdefault("alt_ids", []).append(line[8:])
                elif line.startswith("is_a: "): cur["is_a"].append(line.split()[1])
    children = {t: [] for t in parents}
    for t, pas in parents.items():
        for pa in pas:
            if pa in children:
                children[pa].append(t)

    def genes_for(term):
        seen, dq = {term}, deque([term])
        while dq:
            x = dq.popleft()
            for c in children.get(x, []):
                if c not in seen:
                    seen.add(c); dq.append(c)
        gs = set()
        for t in seen:
            gs |= ann.get(t, set()) | ann.get(alt.get(t, t), set())
        return gs

    return {name: genes_for(goid) for name, goid in terms.items()}

# ---------------------------------------------------------------- datasets
def load_exp002_de():
    """Recompute (or reload cached) GSE303931 DE with symbols."""
    cache = os.path.join(OUT, "gse303931_de_with_symbols.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache)
    tpm, meta = exp002.load_tpm()
    res = exp002.differential_expression(tpm, meta)
    res = exp002.map_symbols(res)
    res.to_csv(cache, index=False)
    return res

def load_gse283507():
    p = os.path.join(ROOT, "data", "gse283507", "GSE283507_raw_Count_FPKM_TPM.csv.gz")
    df = pd.read_csv(p)
    return df

def de_tardbp_dmso(df):
    """TDP43M337V/+ DMSO vs WT(RC802) DMSO, y ~ genotype + time factors."""
    cols = [c for c in df.columns if c.endswith("_TPM") and "DMSO" in c
            and (c.startswith("RC802") or c.startswith("TDP43"))]
    sub = df[["Gene_Symbol"] + cols].dropna(subset=["Gene_Symbol"]).set_index("Gene_Symbol")
    meta = pd.DataFrame({"sample": cols})
    meta["genotype"] = ["ctrl" if s.startswith("RC802") else "mutant" for s in cols]
    meta["time"] = meta["sample"].str.extract(r"(6H|36H|84H)")
    keep = (sub >= 1).sum(axis=1) >= 9
    log = np.log2(sub[keep] + 1)
    T = pd.get_dummies(meta["time"], drop_first=True).astype(float)
    X = np.column_stack([np.ones(len(meta)), (meta["genotype"] == "mutant").astype(float), T.values])
    n, pp = X.shape
    XtXi = np.linalg.pinv(X.T @ X)
    Y = log.values.T
    resid = Y - X @ (XtXi @ X.T @ Y)
    s2 = (resid ** 2).sum(axis=0) / (n - pp)
    beta = XtXi @ X.T @ Y
    tc = beta[1] / np.sqrt(s2 * XtXi[1, 1])
    pv = 2 * stats.t.sf(np.abs(tc), n - pp)
    res = pd.DataFrame({"symbol": log.index, "log2FC": beta[1], "t": tc, "p": pv})
    res["fdr"] = stats.false_discovery_control(res["p"].values) \
        if hasattr(stats, "false_discovery_control") else _bh(res["p"].values)
    return res.sort_values("t", ascending=False).reset_index(drop=True)

def _bh(p):
    p = np.asarray(p); n = len(p)
    o = np.argsort(p); ranked = p[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(ranked, 0, 1)
    return out

# ---------------------------------------------------------------- FC-B
def fc_b(module_genes):
    """Ropinirole within-genotype module response on GSE283507."""
    df = load_gse283507()
    genotypes = {"RC802": "WT_ctrl", "TDP43": "TDP43M337V", "TD-KO": "TDP43M337V_DRD2KO"}
    times = ["6H", "36H", "84H"]
    rows = []
    sym_idx = df["Gene_Symbol"].notna()
    dfm = df[sym_idx].set_index("Gene_Symbol")
    mod = [g for g in module_genes if g in dfm.index]
    for gcol, gname in genotypes.items():
        for t in times:
            dmso = [f"{gcol}DMSO{t}-{r}_TPM" for r in (1, 2, 3)]
            ropi = [f"{gcol}ROPI{t}-{r}_TPM" for r in (1, 2, 3)]
            A = np.log2(dfm.loc[mod, dmso].values + 1).mean(axis=1)
            B = np.log2(dfm.loc[mod, ropi].values + 1).mean(axis=1)
            lfc = B - A                                   # per-module-gene ROPI-DMSO
            # paired-ish test across module genes (genes as units)
            tt = stats.ttest_1samp(lfc, 0.0)
            rows.append({"genotype": gname, "time": t, "module_n": len(mod),
                         "mean_lfc_ropi_dmso": float(np.mean(lfc)),
                         "t": float(tt.statistic), "p": float(tt.pvalue)})
    tab = pd.DataFrame(rows)
    # rescue index R = mean_t(M_ROPI-M_DMSO)_genotype - same for WT (WT contrast ~ vehicle dynamics)
    wt = tab[tab.genotype == "WT_ctrl"].mean_lfc_ropi_dmso.mean()
    out = []
    for gname in ["TDP43M337V", "TDP43M337V_DRD2KO"]:
        m = tab[tab.genotype == gname].mean_lfc_ropi_dmso.mean()
        out.append({"genotype": gname, "rescue_index_R": float(m - wt),
                    "rescue_positive": bool((m - wt) > 0)})
    return tab, pd.DataFrame(out)

# ---------------------------------------------------------------- queries
def build_queries():
    res303 = load_exp002_de()
    df = load_gse283507()
    res283 = de_tardbp_dmso(df)
    res283.to_csv(os.path.join(OUT, "gse283507_de_dmso_ranked.csv"), index=False)
    q1 = [(s, t) for s, t in zip(res303.symbol, res303.t_stat) if isinstance(s, str)]
    q2 = [(s, t) for s, t in zip(res283.symbol, res283.t) if isinstance(s, str)]
    return q1, q2

# ---------------------------------------------------------------- stages
def stage_fc_b(module_genes):
    tab, resc = fc_b(module_genes)
    tab.to_csv(os.path.join(RESULTS, "fcB_ropinirole_module_response_per_time.csv"), index=False)
    resc.to_csv(os.path.join(RESULTS, "fcB_rescue_index.csv"), index=False)
    print("[FC-B]", resc.to_string(index=False))
    n_pos = int(resc.rescue_positive.sum())
    verdict = "PASS" if n_pos >= 2 else "FAIL"
    print(f"[FC-B] rescue-positive genotypes: {n_pos}/2 disease genotypes -> {verdict}")
    return verdict

def stage_score(q1, q2):
    gene_info = pd.read_csv(os.path.join(LINCS_DIR, "GSE92742_Broad_LINCS_gene_info.txt.gz"), sep="\t")
    pr_by_sym = dict(zip(gene_info.pr_gene_symbol, gene_info.pr_gene_id.astype(str)))
    shape, gids, sig_ids = read_gctx_metadata(GCTX)
    print(f"GCTX matrix {shape}, genes={len(gids)}, sigs={len(sig_ids)}")
    v1, m1, t1 = align_query_to_gctx(q1, gids, pr_by_sym)
    v2, m2, t2 = align_query_to_gctx(q2, gids, pr_by_sym)
    print(f"query mapping: q1 {m1}/{t1} landmark genes, q2 {m2}/{t2}")
    Q = np.column_stack([v1, v2])
    scores = stream_scores_multi(GCTX, Q)          # (473647, 2), higher=reversal
    sig_meta = pd.read_csv(os.path.join(LINCS_DIR, "GSE92742_Broad_LINCS_sig_info.txt.gz"),
                           sep="\t", low_memory=False)[
        ["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id", "pert_dose", "pert_itime"]]
    assert len(sig_meta) == len(sig_ids)
    sig_meta = sig_meta.set_index("sig_id").loc[sig_ids].reset_index()
    sc = pd.DataFrame({"sig_id": sig_ids,
                       "rev_q1_c9orf72": scores[:, 0], "rev_q2_tardbp": scores[:, 1]})
    merged = sc.merge(sig_meta, on="sig_id")
    pd.to_pickle(merged, os.path.join(OUT, "level5_reversal_per_signature.pkl"))
    print(f"per-signature scores saved: {merged.shape}")
    np.save(os.path.join(OUT, "query_mapping.json"),
            {"q1_mapped": m1, "q1_total": t1, "q2_mapped": m2, "q2_total": t2}, allow_pickle=False)
    return sc

def stage_agg(sc):
    sig = sc.drop(columns=["sig_id"])
    grp = {c: sig.groupby(["pert_iname", "pert_type"])[c] for c in ["rev_q1_c9orf72", "rev_q2_tardbp"]}
    agg = sig.groupby(["pert_iname", "pert_type"]).agg(
        rev_q1_mean=("rev_q1_c9orf72", "mean"), rev_q2_mean=("rev_q2_tardbp", "mean"),
        n_sigs=("rev_q1_c9orf72", "size")).reset_index()
    agg["rank_q1"] = agg.rev_q1_mean.rank(ascending=False)
    agg["rank_q2"] = agg.rev_q2_mean.rank(ascending=False)
    agg.to_csv(os.path.join(RESULTS, "perturbagen_reversal_rankings.csv"), index=False)

    # FC-D positive control: cross-dataset ranking consistency + label-permutation null
    rng = np.random.default_rng(SEED)
    top1000_union = agg.nlargest(1000, "rev_q1_mean") \
        .merge(agg.nlargest(1000, "rev_q2_mean")[["pert_iname"]], on="pert_iname", how="outer")
    sub = agg[agg.pert_iname.isin(top1000_union.pert_iname)]
    rho_obs = stats.spearmanr(sub.rev_q1_mean, sub.rev_q2_mean).statistic
    n_perm, ge = 200, 0
    vals = sub.rev_q1_mean.values.copy()
    for _ in range(n_perm):
        rng.shuffle(vals)
        rho_null = stats.spearmanr(vals, sub.rev_q2_mean.values).statistic
        if rho_null >= rho_obs:
            ge += 1
    fc_d_p = (ge + 1) / (n_perm + 1)
    pass_fc_d = bool(rho_obs >= 0.10 and fc_d_p <= 0.05)
    json.dump({"rho_observed": float(rho_obs), "perm_p": fc_d_p,
               "pass": pass_fc_d, "n_top1000_union": int(len(sub))},
              open(os.path.join(RESULTS, "fcD_cross_dataset_consistency.json"), "w"), indent=1)
    print(f"[FC-D] spearman rho={rho_obs:.3f}, perm P={fc_d_p:.4f} -> "
          f"{'PASS' if pass_fc_d else 'FAIL'}")
    return agg

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", choices=["fc_b", "build", "score", "agg"], default=None)
    args = ap.parse_args()

    go = load_go_descendants(os.path.join(ROOT, "data", "gse303931", "ref"), MODULE_GO)
    module_genes = sorted(set().union(*go.values()))
    pd.Series(sorted(module_genes)).to_csv(os.path.join(RESULTS, "module_genes.txt"),
                                           index=False, header=["gene"])
    print(f"module genes (union): {len(module_genes)}")

    if args.only in (None, "fc_b"):
        stage_fc_b(module_genes)
    if args.only in (None, "build", "score"):
        q1, q2 = build_queries()
        print(f"queries built: q1={len(q1)} genes, q2={len(q2)} genes")
    if args.only in (None, "score"):
        sc = stage_score(q1, q2)
    else:
        pk = os.path.join(OUT, "level5_reversal_per_signature.pkl")
        sc = pd.read_pickle(pk) if os.path.exists(pk) else None
    if args.only in (None, "agg") and sc is not None:
        stage_agg(sc)

if __name__ == "__main__":
    main()
