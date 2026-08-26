#!/usr/bin/env python
"""exp002_c9orf72_module - H-007 DPR transcriptional correlate module.

End-to-end: DE (C9orf72 mutant vs isogenic control, GSE303931) ->
ranked module lists -> GO footprint enrichment (falsification criterion 1) ->
LINCS Level5 metadata coverage cross-ref.

Usage: <project-python> experiments/exp002_c9orf72_module/run.py
Requires: pandas, numpy, scipy, requests. Raw data in data/gse303931/ (see config.yaml).
"""
import glob, gzip, os, sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DATA = os.path.join(ROOT, "data", "gse303931")
OUT = os.path.join(os.path.dirname(__file__), "outputs")
RESULTS = os.path.join(os.path.dirname(__file__), "results")
LINCS_DIR = os.environ.get("LINCS_DIR",
    os.path.expanduser("~/Project-ALS/data/lincs"))  # read-only, owned by another agent
SEED = 20260825

# ---------- 1. load expression ----------
def load_tpm():
    files = sorted(glob.glob(os.path.join(DATA, "*.salmon_quant.sf.gz")))
    assert len(files) == 12, f"expected 12 quant files in {DATA}, found {len(files)}"
    cols = {}
    for f in files:
        s = os.path.basename(f).replace(".salmon_quant.sf.gz", "")
        df = pd.read_csv(f, sep="\t", index_col=0, usecols=["Name", "TPM"])
        cols[s] = df["TPM"]
    tpm = pd.DataFrame(cols)
    meta = pd.DataFrame({"sample": tpm.columns})
    meta["line"] = meta["sample"].str.extract(r"(C9\d\d)")
    meta["condition"] = np.where(meta["sample"].str.contains("Mutant"), "mutant", "control")
    return tpm, meta

# ---------- 2. differential expression ----------
def differential_expression(tpm, meta):
    keep = (tpm >= 1).sum(axis=1) >= 6
    log = np.log2(tpm[keep] + 1)
    X = np.column_stack([
        np.ones(len(meta)),
        (meta["condition"] == "mutant").astype(float),
        (meta["line"] == "C952").astype(float),
    ])
    Y = log.values.T
    n, p = X.shape
    XtX_inv = np.linalg.pinv(X.T @ X)
    resid = Y - X @ (XtX_inv @ X.T @ Y)
    sigma2 = (resid ** 2).sum(axis=0) / (n - p)
    beta = XtX_inv @ X.T @ Y
    se = np.sqrt(sigma2 * XtX_inv[1, 1])
    t_cond = beta[1] / se
    pvals = 2 * stats.t.sf(np.abs(t_cond), n - p)
    res = pd.DataFrame({
        "ENSG": log.index,
        "log2FC_mutant_vs_ctrl": beta[1],
        "t_stat": t_cond,
        "p_value": pvals,
    })
    res["fdr_bh"] = stats.false_discovery_control(res["p_value"].values) \
        if hasattr(stats, "false_discovery_control") else _bh(res["p_value"].values)
    res = res.sort_values("t_stat", ascending=False).reset_index(drop=True)
    return res

def _bh(p):
    p = np.asarray(p); n = len(p)
    o = np.argsort(p); ranked = p[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(ranked, 0, 1)
    return out

# ---------- 3. symbol mapping ----------
def map_symbols(res):
    """Map ENSG -> HGNC symbol via Ensembl BioMart. Cached in outputs/ because the
    BioMart service is flaky (intermittent 'Service unavailable' HTML pages)."""
    import requests
    cache = os.path.join(OUT, "ensembl_symbol_map.csv")
    res = res.copy()
    res["ENSG_base"] = res["ENSG"].str.split(".").str[0]
    if os.path.exists(cache):
        bm = pd.read_csv(cache, dtype=str)
        print(f"    using cached symbol map ({len(bm)} rows)")
    else:
        q = ('<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>'
             '<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" '
             'datasetConfigVersion="0.6"><Dataset name="hsapiens_gene_ensembl" interface="default">'
             '<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/>'
             '<Attribute name="gene_biotype"/></Dataset></Query>')
        last_err = None
        for attempt in range(5):
            try:
                r = requests.get("http://www.ensembl.org/biomart/martservice",
                                 params={"query": q}, timeout=600)
                if "<html" in r.text[:200].lower():
                    raise RuntimeError("BioMart returned an HTML error page "
                                       "(service unavailable)")
                bm = pd.read_csv(pd.io.common.StringIO(r.text), sep="\t", header=None,
                                 names=["ENSG_base", "symbol", "biotype"], dtype=str)
                bm = bm.dropna(subset=["symbol"]).drop_duplicates("ENSG_base")
                if len(bm) < 20000:
                    raise RuntimeError(f"suspiciously small BioMart table ({len(bm)} rows)")
                bm.to_csv(cache, index=False)
                break
            except Exception as e:
                last_err = e
                print(f"    BioMart attempt {attempt+1} failed: {e}")
        else:
            raise RuntimeError(f"symbol mapping unavailable after retries: {last_err}")
    mapped = res["ENSG_base"].isin(set(bm["ENSG_base"])).mean()
    if mapped < 0.7:
        raise RuntimeError(f"symbol map covers only {mapped:.0%} of queried genes - refusing to continue")
    return res.merge(bm, on="ENSG_base", how="left")

# ---------- 4. GO gene sets ----------
def load_go_sets(refdir):
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
    from collections import deque
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

    return {
        "nucleolus": genes_for("GO:0005730"),
        "ribosome_biogenesis": genes_for("GO:0042254"),
        "nuclear_speckle": genes_for("GO:0016607"),
        "mrna_splicing_major": genes_for("GO:0000398"),
    }

DPR_CURATED = {"EIF2S1","EIF2AK2","EIF2AK3","EIF2AK4","ATF4","DDIT3","XBP1","ERN1",
               "HSPA5","DNAJB9","TP53","CDKN1A","MDM2","NPM1","NCL","NOLC1","DDX21",
               "TCOF1","RRN3","POLR1A","POLR1B","RPS19","RPL11","UBTF","SON","SRSF1",
               "SRSF2","HNRNPA1","HNRNPA2B1","TARDBP"}

# ---------- 5. permutation enrichment ----------
def enrich(bg_syms, module_syms, set_syms, rng, nperm=20000):
    S_mask = np.array([s in set(set_syms) for s in bg_syms], dtype=bool)
    M_mask = np.array([s in set(module_syms) for s in bg_syms], dtype=bool)
    k, m = len(module_syms), int(S_mask.sum())
    assert k > 0 and m > 0, "empty module or gene set - cannot run enrichment"
    obs = int(S_mask[M_mask].sum())
    hits_idx = np.flatnonzero(S_mask)
    exp = m * k / len(bg_syms)
    ge = 0
    for _ in range(nperm):
        draw = rng.choice(len(bg_syms), size=k, replace=False)
        if np.isin(draw, hits_idx).sum() >= max(obs, 1):
            ge += 1
    return {"k_module": k, "m_set": m, "observed_overlap": obs, "expected_overlap": round(exp, 3),
            "fold_enrichment": round(obs / exp, 3) if exp else float("inf"),
            "permutation_p": round((ge + 1) / (nperm + 1), 5)}

# ---------- main ----------
def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    tpm, meta = load_tpm()
    print(f"[1/6] loaded TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")

    res = differential_expression(tpm, meta)
    try:
        res = map_symbols(res)
    except Exception as e:  # offline fallback: keep ENSG ids
        print(f"    WARN symbol mapping failed ({e}); proceeding with ENSG ids")
        res["symbol"], res["biotype"], res["ENSG_base"] = np.nan, np.nan, res["ENSG"].str.split(".").str[0]
    res.to_csv(os.path.join(OUT, "de_ranked_full.csv"), index=False)
    print(f"[2/6] DE done: {(res['fdr_bh']<0.05).sum()} sig at FDR<0.05; full table -> outputs/de_ranked_full.csv")

    # ranked module lists (symbol-ranked where available)
    bg = res.dropna(subset=["symbol"]).reset_index(drop=True)
    modules = {
        "module_fdr10_all": bg[bg.fdr_bh < 0.10],
        "module_fdr05_up": bg[(bg.fdr_bh < 0.05) & (bg.log2FC_mutant_vs_ctrl > 0)],
        "module_fdr05_down": bg[(bg.fdr_bh < 0.05) & (bg.log2FC_mutant_vs_ctrl < 0)],
        "module_top500_absT": pd.concat([bg.nlargest(500, "t_stat"), bg.nsmallest(500, "t_stat")]),
    }
    for name, mod in modules.items():
        mod[["symbol", "ENSG", "log2FC_mutant_vs_ctrl", "t_stat", "p_value", "fdr_bh"]] \
            .to_csv(os.path.join(OUT, f"{name}_ranked.csv"), index=False)
    # plain ranked lists (GCT-friendly: symbol \t score)
    with open(os.path.join(OUT, "module_ranked_list.tsv"), "w") as fh:
        fh.write("module\tgene\tlog2FC\tt_stat\tfdr_bh\n")
        for name, mod in modules.items():
            for _, row in mod.iterrows():
                fh.write(f"{name}\t{row.symbol}\t{row.log2FC_mutant_vs_ctrl:.4f}\t"
                         f"{row.t_stat:.4f}\t{row.fdr_bh:.5g}\n")
    print(f"[3/6] wrote 4 ranked module lists -> outputs/")

    # pathway sets
    refdir = os.path.join(DATA, "ref")
    go = load_go_sets(refdir)
    nucleolar_set = go["nucleolus"] | go["ribosome_biogenesis"]
    speckle_set = go["nuclear_speckle"] | go["mrna_splicing_major"]
    rng = np.random.default_rng(SEED)
    rows = []
    for mname, mod in modules.items():
        syms = list(mod.symbol.dropna())
        for sname, ss in [("nucleolar+ribo_biogenesis", nucleolar_set),
                          ("speckle+splicing", speckle_set),
                          ("DPR_curated_literature", DPR_CURATED)]:
            rows.append({"module": mname, "gene_set": sname,
                         **enrich(bg.symbol.values, syms, ss, rng)})
    enr = pd.DataFrame(rows)
    enr.to_csv(os.path.join(RESULTS, "enrichment.csv"), index=False)
    fc1 = enr[(enr.gene_set != "DPR_curated_literature")]
    triggered = ((fc1.fold_enrichment < 2.0) & (fc1.permutation_p > 0.05)).all()
    print(f"[4/6] footprint enrichment -> results/enrichment.csv ; "
          f"H-007 falsification criterion 1 TRIGGERED: {triggered}")

    # LINCS metadata-only cross-ref (read-only; gctx owned by another agent)
    lrows = []
    try:
        gene_info = pd.read_csv(os.path.join(LINCS_DIR, "GSE92742_Broad_LINCS_gene_info.txt.gz"), sep="\t")
        lm = set(gene_info.loc[gene_info.pr_is_lm == 1, "pr_gene_symbol"])
        bing = set(gene_info.pr_gene_symbol)
        for name, mod in modules.items():
            syms = set(mod.symbol.dropna())
            lrows.append({"module": name, "n_module_genes": len(syms),
                          "n_in_lincs_measured": len(syms & bing),
                          "n_in_lincs_landmark978": len(syms & lm)})
    except FileNotFoundError as e:
        print(f"    LINCS metadata not accessible ({e}); dependency noted.")
    lin = pd.DataFrame(lrows)
    lin.to_csv(os.path.join(RESULTS, "lincs_crossref_metadata_only.csv"), index=False)
    print("[5/6] LINCS metadata cross-ref -> results/lincs_crossref_metadata_only.csv "
          "(full reversal scoring deferred to pipelines/perturbation_signatures/)")

    top = res.head(25)[["symbol", "log2FC_mutant_vs_ctrl", "t_stat", "fdr_bh"]]
    top.to_csv(os.path.join(RESULTS, "top_genes_by_t.csv"), index=False)
    print("[6/6] top genes -> results/top_genes_by_t.csv")
    return 0

if __name__ == "__main__":
    sys.exit(main())
