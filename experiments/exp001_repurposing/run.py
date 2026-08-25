#!/usr/bin/env python
"""Experiment 001: transcriptome-based drug-repurposing screen for ALS.

Pipeline:
  1. Parse GSE124439 (post-mortem CNS RNA-seq) sample metadata + counts.
  2. Differential expression ALS vs control (limma-style moderated t, numpy/scipy).
  3. Build ranked disease signature; sanity-check against known ALS biology genes.
  4. Stream-score every LINCS L1000 Phase-1 Level5 signature (GSE92742) against it
     by cosine reversal (-cosine on shared inferred genes).
  5. Positive controls: rank known neuro/immune-active compounds (incl. riluzole).
  6. Annotate top candidates with ChEMBL max_phase + PubChem properties.

Run:  python run.py --config config.yaml   (from this directory)
"""
import argparse, gzip, json, re, sys, time
from urllib.parse import quote
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ----------------------------------------------------------------------------
# 1. Disease data: GSE124439
# ----------------------------------------------------------------------------

def parse_series_matrix(path):
    """Extract per-sample title / group / subregion from GEO series matrix."""
    rows, char_rows = {}, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("!Sample_characteristics_ch1"):
                parts = line.rstrip("\n").split("\t")
                char_rows.append([p.strip().strip('"') for p in parts[1:]])
            elif line.startswith("!Sample_"):
                parts = line.rstrip("\n").split("\t")
                rows[parts[0]] = [p.strip().strip('"') for p in parts[1:]]
    titles = rows["!Sample_title"]
    accs = rows["!Sample_geo_accession"]
    chars = char_rows
    groups, subregions = [], []
    for i in range(len(titles)):
        g = s = None
        for ch in chars:
            if ch.startswith("sample group:"):
                g = ch.split(":", 1)[1].strip()
            elif ch.startswith("cns subregion:"):
                s = ch.split(":", 1)[1].strip()
        groups.append(g)
        subregions.append(s)
    meta = pd.DataFrame({
        "geo_accession": accs, "title": titles,
        "group": groups, "subregion": subregions,
    })
    return meta


def build_count_matrix(raw_dir, meta):
    """Read per-sample gene-symbol count files into one matrix."""
    series = {}
    for _, row in meta.iterrows():
        path = Path(raw_dir) / f"{row.geo_accession}_{row.title}_counts.txt.gz"
        df = pd.read_csv(path, sep="\t", skiprows=1, header=None,
                         names=["gene", row.title], dtype={0: str})
        df["gene"] = df["gene"].str.strip('"')
        series[row.title] = df.set_index("gene")[row.title]
    mat = pd.DataFrame(series)
    return mat


def differential_expression(counts, is_case, min_cpm, min_samples):
    """Log2-CPM OLS + limma-style empirical-Bayes moderated t (Smye 2004)."""
    lib = counts.sum(axis=0).values
    cpm = counts.div(lib, axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    cpm = cpm.loc[keep]
    X = np.log2(cpm + 1)
    n_c, n_t = int((~is_case).sum()), int(is_case.sum())
    y = np.log2(counts.loc[keep].div(lib, axis=1) * 1e6 + 1)

    # design: [1, case]
    D = np.column_stack([np.ones(X.shape[1]), is_case.astype(float)])
    pinv = np.linalg.pinv(D)
    beta = y.values @ pinv.T                      # genes x 2
    fitted = beta @ D.T
    resid = y.values - fitted
    df_res = y.shape[1] - 2
    rss = (resid ** 2).sum(axis=1)
    s2 = rss / df_res

    # empirical Bayes: inverse-gamma prior on s2 via method of moments
    m, v = s2.mean(), s2.var(ddof=1)
    if v <= 0:
        d0, s2_0 = np.inf, m
    else:
        a = 2.0 + m * m / v                       # shape of IG prior (d0/2)
        d0 = 2.0 * a
        s2_0 = m * (a - 1.0) / a
    d0 = max(d0, 1.0)
    s2_post = (d0 * s2_0 + df_res * s2) / (d0 + df_res)

    cov_case = np.linalg.inv(D.T @ D)[1, 1]
    se_post = np.sqrt(s2_post * cov_case)
    tstat = beta[:, 1] / se_post
    pval = 2 * stats.t.sf(np.abs(tstat), df=d0 + df_res)
    adj = bh_adjust(pval)

    res = pd.DataFrame({
        "gene": y.index, "log2FC": beta[:, 1],
        "t": tstat, "p": pval, "padj": adj,
    })
    res = res.sort_values("t", ascending=False).reset_index(drop=True)
    return res


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


# ----------------------------------------------------------------------------
# 3-4. LINCS L1000 scoring
# ----------------------------------------------------------------------------

def read_gctx_metadata(gctx_path):
    import h5py
    with h5py.File(gctx_path, "r") as f:
        shape = _descend_to_matrix(f).shape
        row_ids = _read_meta_str(f, ["metadata", "rows", "id"])
        col_ids = _read_meta_str(f, ["metadata", "cols", "id"])
    return shape, row_ids, col_ids


def _descend_to_matrix(f):
    node = f
    for _ in range(6):
        if isinstance(node, h5py.Dataset):
            return node
        node = node[list(node.keys())[0]]
    raise RuntimeError("could not locate matrix dataset")


def _read_meta_str(f, path):
    node = f
    for p in path:
        node = node[p]
    arr = node[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def stream_scores(gctx_path, query_vec, gene_block=500, progress_every=50):
    """Cosine similarity between query_vec (aligned to gctx rows) and all columns.

    Returns (scores ndarray float32, n_genes_used). Single sequential pass.
    """
    import h5py
    with h5py.File(gctx_path, "r") as f:
        mat = _descend_to_matrix(f)
        g, s = mat.shape
        assert mat.shape[0] == len(query_vec), (mat.shape, len(query_vec))
        q = query_vec.astype(np.float32)
        qn = np.linalg.norm(q)
        dots = np.zeros(s, dtype=np.float64)
        norms = np.zeros(s, dtype=np.float64)
        t0 = time.time()
        for a in range(0, g, gene_block):
            b = min(a + gene_block, g)
            Y = mat[a:b, :].astype(np.float64)      # genes_block x sigs
            dots += q[a:b] @ Y
            norms += (Y * Y).sum(axis=0)
            done = b / g
            el = time.time() - t0
            if (a // gene_block) % progress_every == 0 or b == g:
                eta = el / done - el if done > 0 else float("nan")
                log(f"scoring {b}/{g} gene blocks ({done:.1%}), elapsed {el:.0f}s, ETA {eta:.0f}s")
    scores = -(dots / (qn * np.maximum(norms, 1e-12) ** 0.5))  # reversal strength
    return scores.astype(np.float32), g


# ----------------------------------------------------------------------------
# Positive control + annotation
# ----------------------------------------------------------------------------

ALIASES = {
    "rapamycin": "sirolimus",
    "cyclosporine": "cyclosporin a",
    "aspirin": "aspirin",
}

CNS_KNOWN = {"riluzole", "edaravone", "memantine", "valproic acid", "gabapentin",
             "baclofen", "diazepam", "levetiracetam", "topiramate", "lamotrigine",
             "haloperidol", "clozapine", "fluoxetine", "sertraline", "donepezil",
             "levodopa", "carbidopa", "pramipexole", "amantadine", "zolpidem"}


def positive_control_table(ranked_drugs, compounds):
    """ranked_drugs: DataFrame indexed by lowercase pert_iname with 'rank' column."""
    rows = []
    n = len(ranked_drugs)
    for comp in compounds:
        name = ALIASES.get(comp.lower(), comp.lower())
        if name in ranked_drugs.index:
            r = int(ranked_drugs.loc[name, "rank"])
            pct = 100.0 * r / n
            verdict = "reverses" if pct <= 10 else ("neutral" if pct <= 60 else "mimics")
            rows.append({"query_name": comp, "lincs_pert_iname": name, "rank": r,
                         "percentile": round(pct, 1), "direction": verdict})
        else:
            rows.append({"query_name": comp, "lincs_pert_iname": name, "rank": None,
                         "percentile": None, "direction": "not tested in LINCS phase 1"})
    return pd.DataFrame(rows)


def chembl_lookup(name, session, url, cache):
    key = f"chembl::{name.lower()}"
    if key in cache:
        return cache[key]
    try:
        r = session.get(url, params={"pref_name__iexact": name, "limit": 1}, timeout=30)
        r.raise_for_status()
        mols = r.json().get("molecules", [])
        out = {}
        if mols:
            mol = mols[0]
            out = {"chembl_id": mol.get("molecule_chembl_id"),
                   "max_phase": mol.get("max_phase"),
                   "first_approval": mol.get("first_approval"),
                   "pref_name": mol.get("pref_name")}
    except Exception as e:
        out = {"error": str(e)[:120]}
    cache[key] = out
    time.sleep(0.35)
    return out


PUBCHEM_PROPS = "MolecularFormula,MolecularWeight,CanonicalSMILES,HydrogenBondDonorCount"


def pubchem_lookup(name, session, url, cache):
    key = f"pubchem::{name.lower()}"
    if key in cache:
        return cache[key]
    try:
        u = url.rstrip("/") + f"/name/{quote(name)}/properties/{PUBCHEM_PROPS}/JSON"
        r = session.get(u, timeout=30)
        if r.status_code == 404:
            out = {}
        else:
            r.raise_for_status()
            out = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
    except Exception as e:
        out = {"error": str(e)[:120]}
    cache[key] = out
    time.sleep(0.35)
    return out


# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parents[2]
    cfg = yaml.safe_load((Path(__file__).parent / args.config).read_text())
    outdir = Path(__file__).parent / cfg["outputs"]["dir"]
    outdir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    def rp(p):
        return repo / p

    # ---- 1-2. disease data + DE ----
    log("parsing GSE124439 metadata")
    d = cfg["disease"]
    meta = parse_series_matrix(rp(d["series_matrix"]))
    counts = build_count_matrix(rp(d["raw_dir"]), meta)
    meta.to_csv(outdir / "sample_sheet.csv", index=False)
    log(f"count matrix: {counts.shape[0]} genes x {counts.shape[1]} samples")

    groups = meta.set_index("title")["group"].loc[counts.columns]
    keep_groups = [d["case_label"], d["control_label"]]
    mask_sample = groups.isin(keep_groups).values
    counts = counts.loc[:, mask_sample]
    groups = groups[mask_sample]
    is_case = (groups == d["case_label"]).values
    log(f"DE contrast: {is_case.sum()} ALS vs {(~is_case).sum()} control "
        f"(excluded {mask_sample.size - mask_sample.sum()} other-neuro samples)")

    res = differential_expression(counts, is_case, d["min_count_per_million"], d["min_samples"])
    res.to_csv(outdir / "disease_signature_full.csv", index=False)

    up = res.head(50)
    dn = res.tail(50).iloc[::-1]
    pd.concat([up.assign(direction="up"), dn.assign(direction="down")]) \
        .to_csv(outdir / "als_signature_top.csv", index=False)

    # sanity check against known ALS biology
    sanity_genes = ["NEFH", "NEFM", "MAPT", "TARDBP", "SQSTM1", "OPTN", "TBK1", "C9orf72",
                    "AIF1", "TYROBP", "FCGR3A", "CX3CR1", "C1QA", "C1QB", "CD68", "HLA-DRA",
                    "GFAP", "SOD1", "FUS", "VCP", "MATR3", "TUBA4A"]
    sc = res[res.gene.isin(sanity_genes)].set_index("gene")
    sanity_lines = []
    for g in sanity_genes:
        if g in sc.index:
            r = sc.loc[g]
            sanity_lines.append(f"| {g} | {r.log2FC:+.3f} | {r.t:.2f} | {r.padj:.2e} |")
        else:
            sanity_lines.append(f"| {g} | filtered | - | - |")
    immune_up = res[(res.gene.isin(["AIF1", "TYROBP", "C1QA", "C1QB", "CD68", "FCGR3A"])) & (res.log2FC > 0)]
    neurofil_down = res[res.gene.isin(["NEFH", "NEFM"]) & (res.log2FC < 0)]
    sanity_ok = len(immune_up) >= 3 and True  # microglial activation expected in ALS CNS
    (outdir / "sanity_check.md").write_text(
        "# ALS signature sanity check\n\n"
        "| gene | log2FC | moderated t | padj |\n|---|---|---|---|\n" +
        "\n".join(sanity_lines) +
        f"\n\nMicroglial/immune genes up in ALS: {len(immune_up)}\n"
        f"Neurofilament genes down in ALS: {len(neurofil_down)}\n"
    )
    log(f"sanity check: {len(immune_up)} immune genes up, NEFH/NEFM down={len(neurofil_down)}")

    # ---- 3-4. LINCS scoring ----
    l = cfg["lincs"]
    gctx = rp(l["level5_gctx"])
    size = gctx.stat().st_size
    expected = 21328033748
    if size != expected:
        log(f"WARNING: GCTX file is {size} bytes, expected {expected} - download incomplete?")
    gene_info = pd.read_csv(rp(l["gene_info"]), sep="\t")
    sym2row = pd.Series(gene_info.index.values, index=gene_info["pr_gene_symbol"].str.upper()).to_dict()
    shared = res[res.gene.str.upper().isin(sym2row)]
    q = np.full(12328, np.nan)
    for _, r in shared.iterrows():
        q[sym2row[r.gene.upper()]] = r.t
    valid = ~np.isnan(q)
    q[~valid] = 0.0
    log(f"disease signature mapped to {valid.sum()} L1000 genes")

    log("stream-scoring LINCS Level5 signatures")
    scores, ngene = stream_scores(gctx, q, gene_block=l["gene_block_size"])
    _, row_ids, col_ids = read_gctx_metadata(gctx)
    sig_scores = pd.DataFrame({"sig_id": col_ids, "score_reversal": scores})

    sig_info = pd.read_csv(rp(l["sig_info"]), sep="\t")
    merged = sig_scores.merge(
        sig_info[["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id", "pert_idose"]],
        on="sig_id", how="left")
    cp = merged[merged.pert_type == l["pert_type_filter"]].copy()
    cp["rank_sig"] = cp.score_reversal.rank(method="first").astype(int)
    cp.sort_values("score_reversal").to_csv(outdir / "all_signature_scores_trt_cp.csv", index=False)

    drug = (cp.groupby(cp.pert_iname.str.lower())
              .agg(best_score=("score_reversal", "min"),
                   median_score=("score_reversal", "median"),
                   n_sigs=("sig_id", "count"))
              .sort_values("best_score").reset_index()
              .rename(columns={"pert_iname": "pert_iname_lower"}))
    drug["rank"] = np.arange(1, len(drug) + 1)
    drug.to_csv(outdir / "drug_ranking.csv", index=False)
    log(f"scored {len(merged)} signatures; {len(drug)} distinct small-molecule perturbagens")

    # ---- 5. positive control ----
    pc = positive_control_table(drug.set_index("pert_iname_lower"), cfg["positive_controls"]["compounds"])
    pc.to_csv(outdir / "positive_controls.csv", index=False)
    ril = pc[pc.query_name == "riluzole"].iloc[0]
    log(f"positive control: riluzole rank {ril['rank']} ({ril['percentile']}th pct, {ril['direction']})")

    # ---- 6. annotation ----
    import requests
    session = requests.Session()
    ann_cfg = cfg["annotation"]
    cache_path = Path(__file__).parent / ann_cfg["cache"]
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}

    top = drug.head(ann_cfg["top_n_drugs"]).copy()
    recs = []
    for _, r in top.iterrows():
        disp = r.pert_iname_lower.title() if r.pert_iname_lower.islower() else r.pert_iname_lower
        ch = chembl_lookup(disp, session, ann_cfg["chembl_url"], cache)
        pu = pubchem_lookup(disp, session, ann_cfg["pubchem_url"], cache)
        mw = pu.get("MolecularWeight")
        hbd = pu.get("HydrogenBondDonorCount")
        cns = (r.pert_iname_lower in CNS_KNOWN) or (
            mw is not None and hbd is not None and float(mw) < 450 and int(hbd) <= 3)
        recs.append({
            "pert_iname": r.pert_iname_lower, "rank": int(r["rank"]),
            "best_score": r.best_score, "n_sigs": int(r.n_sigs),
            "chembl_id": ch.get("chembl_id"), "max_phase": ch.get("max_phase"),
            "first_approval": ch.get("first_approval"),
            "MW": mw, "HBD": hbd, "formula": pu.get("MolecularFormula"),
            "cns_flag": bool(cns),
        })
    cache_path.write_text(json.dumps(cache, indent=1))
    annot = pd.DataFrame(recs)
    annot.to_csv(outdir / "top_candidates_annotated.csv", index=False)
    log(f"annotated top {len(annot)} candidates")

    # ---- summary ----
    summary = {
        "n_samples_als": int(is_case.sum()), "n_samples_control": int((~is_case).sum()),
        "n_genes_tested": int(len(res)),
        "n_lincs_signatures_scored": int(len(merged)),
        "n_small_molecule_perturbagens": int(len(drug)),
        "riluzole_rank": int(ril["rank"]) if ril["rank"] is not None else None,
        "riluzole_percentile": ril["percentile"],
        "wall_time_sec": round(time.time() - t_start, 1),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"DONE wall time {summary['wall_time_sec']}s")


if __name__ == "__main__":
    main()
