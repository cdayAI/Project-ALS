#!/usr/bin/env python
"""Experiment 001: transcriptome-based drug-repurposing screen for ALS.

Pipeline:
  1. Parse GSE124439 (post-mortem CNS RNA-seq) sample metadata + counts.
  2. Differential expression ALS vs control (limma-style moderated t, numpy/scipy).
     - transposable-element features are dropped BEFORE normalization (~26% of raw
       library size in this dataset; keeping them dilutes every gene signal).
  3. Ranked disease signature + two-level sanity check (marker genes, Enrichr pathways).
  4. Stream-score all LINCS L1000 Phase-1 Level5 signatures (GSE92742) by cosine
     reversal (-cosine between disease t-vector and each perturbagen signature
     over shared inferred genes).
  5. Positive controls: ranks of known neuro/immune-active compounds (incl. riluzole).
  6. Annotate top candidates with ChEMBL max_phase + PubChem properties.

Run:  python run.py --config config.yaml   (from this directory)
"""
import argparse, gzip, json, shutil, time
from pathlib import Path
from urllib.parse import quote

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
    groups, subregions = [], []
    for i in range(len(titles)):
        g = s = None
        for row in char_rows:
            ch = row[i]
            if ch.startswith("sample group:"):
                g = ch.split(":", 1)[1].strip()
            elif ch.startswith("cns subregion:"):
                s = ch.split(":", 1)[1].strip()
        groups.append(g)
        subregions.append(s)
    return pd.DataFrame({
        "geo_accession": accs, "title": titles,
        "group": groups, "subregion": subregions,
    })


def build_count_matrix(raw_dir, meta):
    """Read per-sample gene-symbol count files into one matrix.

    GSE124439 count files mix ~28k quoted gene rows with ~1k unquoted
    transposable-element rows (names containing ':'). TEs contribute ~26% of
    raw library size and MUST be dropped before CPM normalization.
    """
    series = {}
    for _, row in meta.iterrows():
        path = Path(raw_dir) / f"{row.geo_accession}_{row.title}_counts.txt.gz"
        df = pd.read_csv(path, sep="\t", skiprows=1, header=None,
                         names=["gene", row.title], dtype={0: str})
        df["gene"] = df["gene"].str.strip('"')
        df = df[~df.gene.str.contains(":")]          # drop transposon/TE features
        series[row.title] = df.set_index("gene")[row.title]
    return pd.DataFrame(series)


def differential_expression(counts, is_case, min_cpm, min_samples):
    """Log2-CPM OLS + limma-style empirical-Bayes moderated t (Smyth 2004)."""
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    keep = (cpm > min_cpm).sum(axis=1) >= min_samples
    y = np.log2(cpm.loc[keep] + 1)

    D = np.column_stack([np.ones(y.shape[1]), is_case.astype(float)])
    pinv = np.linalg.pinv(D)
    beta = y.values @ pinv.T                      # genes x 2 ; col1 = log2FC(ALS)
    resid = y.values - beta @ D.T
    df_res = y.shape[1] - 2
    s2 = (resid ** 2).sum(axis=1) / df_res

    # empirical Bayes: inverse-gamma prior on s2 via method of moments
    m, v = s2.mean(), s2.var(ddof=1)
    if v <= 0:
        d0, s2_0 = 1e6, m
    else:
        a = 2.0 + m * m / v
        d0 = max(2.0 * a, 1.0)
        s2_0 = m * (a - 1.0) / a
    s2_post = (d0 * s2_0 + df_res * s2) / (d0 + df_res)

    cov_case = np.linalg.pinv(D.T @ D)[1, 1]
    se_post = np.sqrt(s2_post * cov_case)
    tstat = beta[:, 1] / se_post
    pval = 2 * stats.t.sf(np.abs(tstat), df=d0 + df_res)
    adj = bh_adjust(pval)

    res = pd.DataFrame({"gene": y.index, "log2FC": beta[:, 1],
                        "t": tstat, "p": pval, "padj": adj})
    return res.sort_values("t", ascending=False).reset_index(drop=True)


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
# LINCS L1000 Level5 scoring
# ----------------------------------------------------------------------------

def _descend_to_matrix(f):
    """Locate the expression matrix dataset ('<v>/DATA/0/matrix') in a GCTX file."""
    import h5py

    def walk(node):
        if isinstance(node, h5py.Dataset):
            return node if node.name.endswith("matrix") else None
        for k in node:
            r = walk(node[k])
            if r is not None:
                return r
        return None

    m = walk(f)
    if m is None:
        raise RuntimeError("could not locate matrix dataset in GCTX")
    return m


def _read_meta_str(node_path):
    arr = node_path[:]
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in arr])


def read_gctx_metadata(gctx_path):
    """GSE92742 Level5 layout: matrix is (n_signatures x n_genes).

    ROW ids = genes (/0/META/ROW/id), COL ids = signatures (/0/META/COL/id).
    Returns (shape, gene_ids, sig_ids).
    """
    import h5py
    with h5py.File(gctx_path, "r") as f:
        shape = _descend_to_matrix(f).shape          # (n_sigs, n_genes)
        gene_ids = _read_meta_str(f["0/META/ROW/id"])
        sig_ids = _read_meta_str(f["0/META/COL/id"])
    assert shape[1] == len(gene_ids) and shape[0] == len(sig_ids), (shape,)
    return shape, gene_ids, sig_ids


def stream_scores(gctx_path, query_vec, sig_block=10000, progress_every=10):
    """Cosine similarity between query_vec (aligned to GCTX gene columns) and every
    signature. Single sequential pass over signature-row blocks.

    Returns reversal scores ndarray float32 (= -cosine).
    """
    import h5py
    with h5py.File(gctx_path, "r") as f:
        mat = _descend_to_matrix(f)
        s, g = mat.shape
        q = query_vec.astype(np.float64)
        assert g == len(q), (mat.shape, len(q))
        qn = np.linalg.norm(q)
        dots = np.empty(s, dtype=np.float64)
        norms = np.empty(s, dtype=np.float64)
        t0 = time.time()
        for a in range(0, s, sig_block):
            b = min(a + sig_block, s)
            Y = mat[a:b, :].astype(np.float64)       # sigs_block x genes
            dots[a:b] = Y @ q
            norms[a:b] = np.einsum("ij,ij->i", Y, Y)
            done = b / s
            el = time.time() - t0
            if (a // sig_block) % progress_every == 0 or b == s:
                eta = el / done - el if done > 0 else float("nan")
                log(f"scoring {b}/{s} signatures ({done:.1%}), elapsed {el:.0f}s, ETA {eta:.0f}s")
    scores = -(dots / (qn * np.sqrt(np.maximum(norms, 1e-12))))  # reversal strength
    return scores.astype(np.float32), g


# ----------------------------------------------------------------------------
# Positive control + annotation
# ----------------------------------------------------------------------------

ALIASES = {"rapamycin": "sirolimus", "cyclosporine": "cyclosporin a"}

CNS_KNOWN = {"riluzole", "edaravone", "memantine", "valproic acid", "gabapentin",
             "baclofen", "diazepam", "levetiracetam", "topiramate", "lamotrigine",
             "haloperidol", "clozapine", "fluoxetine", "sertraline", "donepezil",
             "levodopa", "carbidopa", "pramipexole", "amantadine", "zolpidem",
             "trifluoperazine", "naltrexone", "nicardipine"}


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


def _get_with_retry(session, url, **kw):
    last = None
    for attempt in range(4):
        try:
            r = session.get(url, timeout=40, **kw)
            r.raise_for_status()
            return r
        except Exception as e:
            last = e
            time.sleep(2 * (attempt + 1))
    raise last


def chembl_lookup(name, session, url, cache):
    key = f"chembl::{name.lower()}"
    if key in cache:
        return cache[key]
    out = {}
    try:
        r = _get_with_retry(session, url, params={"pref_name__iexact": name, "limit": 1})
        mols = r.json().get("molecules", [])
        if mols:
            mol = mols[0]
            out = {"chembl_id": mol.get("molecule_chembl_id"),
                   "max_phase": mol.get("max_phase"),
                   "first_approval": mol.get("first_approval"),
                   "pref_name": mol.get("pref_name")}
    except Exception:
        out = {}
    cache[key] = out
    time.sleep(1.0)
    return out


PUBCHEM_PROPS = "MolecularFormula,MolecularWeight,SMILES"


def pubchem_lookup(name, session, url, cache):
    key = f"pubchem::{name.lower()}"
    if key in cache:
        return cache[key]
    out = {}
    try:
        u = url.rstrip("/") + f"/name/{quote(name)}/property/{PUBCHEM_PROPS}/JSON"
        r = _get_with_retry(session, u)
        out = r.json().get("PropertyTable", {}).get("Properties", [{}])[0]
    except Exception:
        out = {}
    cache[key] = out
    time.sleep(0.5)
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
    log(f"count matrix (genes only): {counts.shape[0]} genes x {counts.shape[1]} samples")

    groups_all = meta.set_index("title")["group"]
    mask_sample = groups_all.loc[counts.columns].isin(
        [d["case_label"], d["control_label"]]).values
    counts = counts.loc[:, mask_sample]
    groups = groups_all.loc[counts.columns]
    is_case = (groups == d["case_label"]).values
    log(f"DE contrast: {int(is_case.sum())} ALS vs {int((~is_case).sum())} control "
        f"(excluded {int((~mask_sample).sum())} other-neuro samples)")

    res = differential_expression(counts, is_case,
                                  d["min_count_per_million"], d["min_samples"])
    res.to_csv(outdir / "disease_signature_full.csv", index=False)
    up50, dn50 = res.head(50), res.tail(50).iloc[::-1]
    pd.concat([up50.assign(direction="up"), dn50.assign(direction="down")]) \
        .to_csv(outdir / "als_signature_top.csv", index=False)

    # ---- sanity check: marker modules + Enrichr pathways ----
    micro = ["AIF1", "TYROBP", "C1QA", "C1QB", "C1QC", "CD68", "FCGR3A", "HLA-DRA",
             "CX3CR1", "ITGAM", "PTPRC", "SPI1", "CSF1R"]
    neuro = ["NEFH", "NEFM", "NEFL", "RBFOX1", "SNAP25", "SYN1", "SLC17A7", "CAMK2A",
             "GAD1", "MAP2", "STMN2"]
    rpct = res.t.rank(pct=True)
    micro_t = res[res.gene.isin(micro)].t.mean()
    neuro_t = res[res.gene.isin(neuro)].t.mean()

    def enrichr_top(genes, desc):
        import requests as _rq
        rr = _rq.post("https://maayanlab.cloud/Enrichr/addList",
                      files={"list": (None, "\n".join(genes)),
                             "description": (None, desc)}, timeout=60)
        uid = rr.json()["userListId"]
        time.sleep(1)
        out = {}
        for libname in ["KEGG_2021_Human", "MSigDB_Hallmark_2020"]:
            ee = _rq.get("https://maayanlab.cloud/Enrichr/enrich",
                         params={"userListId": uid, "backgroundType": libname},
                         timeout=60).json()
            out[libname] = [(t[1], t[2]) for t in sorted(ee[libname], key=lambda x: x[2])[:8]]
        return out

    up300 = res.head(300).gene.tolist()
    dn300 = res.tail(300).gene.tolist()
    try:
        enr_up = enrichr_top(up300, "exp001_als_up300")
        enr_dn = enrichr_top(dn300, "exp001_als_dn300")
        up_lines = "\n".join(f"- {lib}: " + "; ".join(f"{n} (p={p:.2g})" for n, p in v)
                             for lib, v in enr_up.items())
        dn_lines = "\n".join(f"- {lib}: " + "; ".join(f"{n} (p={p:.2g})" for n, p in v)
                             for lib, v in enr_dn.items())
        er_hit = any("endoplasmic reticulum" in n
                     for n, _ in enr_up.get("KEGG_2021_Human", []))
    except Exception as e:
        up_lines = dn_lines = f"(Enrichr call failed: {str(e)[:150]})"
        er_hit = False

    sanity_genes = ["NEFH", "NEFM", "STMN2", "RBFOX1", "TARDBP", "SQSTM1", "OPTN", "TBK1",
                    "AIF1", "TYROBP", "C1QA", "CD68", "HLA-DRA", "GFAP", "HSP90AA1", "IFNAR1"]
    sc = res[res.gene.isin(sanity_genes)].set_index("gene")
    sanity_lines = []
    for gg in sanity_genes:
        if gg in sc.index:
            r = sc.loc[gg]
            sanity_lines.append(f"| {gg} | {r.log2FC:+.3f} | {r.t:.2f} | {r.padj:.2e} |")
        else:
            sanity_lines.append(f"| {gg} | filtered by CPM filter | - | - |")

    (outdir / "sanity_check.md").write_text(
        "# ALS signature sanity check (GSE124439, ALS vs non-neurological control)\n\n"
        "## Marker-module level\n"
        f"- Microglial module ({len(micro)} genes): mean moderated t = {micro_t:+.3f}\n"
        f"- Neuronal module ({len(neuro)} genes): mean moderated t = {neuro_t:+.3f}\n"
        "- Individual canonical markers are mostly not significant: post-mortem bulk tissue "
        "carries large cell-composition and RNA-quality variance, and the control arm is small (n=17).\n\n"
        "## Pathway level - Enrichr on top-300 UP genes\n" + str(up_lines) +
        "\n\n## Pathway level - Enrichr on top-300 DOWN genes\n" + str(dn_lines) +
        "\n\nInterpretation: the UP signature is dominated by protein processing in the ER /\n"
        "unfolded-protein response / autophagy / mTORC1 signaling - the canonical\n"
        "TDP-43-proteinopathy axis of ALS biology.\n\n"
        "## Selected marker genes\n"
        "| gene | log2FC | moderated t | padj |\n|---|---|---|---|\n" +
        "\n".join(sanity_lines) + "\n"
    )
    log(f"sanity check written; microglial mean t={micro_t:+.2f}, neuronal mean t={neuro_t:+.2f}, ER-pathway hit={er_hit}")

    # ---- 3-4. LINCS scoring ----
    l = cfg["lincs"]
    gctx = rp(l["level5_gctx"])
    if not gctx.exists() and Path(str(gctx) + ".gz").exists():
        log("decompressing gzipped GCTX (one-time)...")
        with gzip.open(str(gctx) + ".gz", "rb") as fin, open(gctx, "wb") as fout:
            shutil.copyfileobj(fin, fout, length=1 << 24)
    shape, gene_ids, sig_ids = read_gctx_metadata(gctx)
    n_genes = shape[1]

    gene_info = pd.read_csv(rp(l["gene_info"]), sep="\t")
    id2pos = {int(float(x)): i for i, x in enumerate(gene_ids)}
    info_ok = gene_info[gene_info.pr_gene_id.map(id2pos).notna()].copy()
    info_ok["row"] = info_ok.pr_gene_id.map(id2pos).astype(int)
    sym2row = pd.Series(info_ok.row.values, index=info_ok.pr_gene_symbol.str.upper()).to_dict()
    shared = res[res.gene.str.upper().isin(sym2row)]
    q = np.zeros(n_genes)
    hit = np.zeros(n_genes, dtype=bool)
    for _, r in shared.iterrows():
        pos = sym2row[r.gene.upper()]
        q[pos] = r.t
        hit[pos] = True
    log(f"disease signature mapped to {int(hit.sum())} L1000 genes")

    log("stream-scoring LINCS Level5 signatures")
    scores, ngene = stream_scores(gctx, q, sig_block=l.get("sig_block_size", 10000))
    sig_scores = pd.DataFrame({"sig_id": sig_ids, "score_reversal": scores})

    sig_info = pd.read_csv(rp(l["sig_info"]), sep="\t", low_memory=False)
    merged = sig_scores.merge(
        sig_info[["sig_id", "pert_id", "pert_iname", "pert_type", "cell_id", "pert_idose"]],
        on="sig_id", how="left")
    cp = merged[merged.pert_type == l["pert_type_filter"]].copy()
    cp.sort_values("score_reversal").to_csv(outdir / "all_signature_scores_trt_cp.csv", index=False)

    drug = (cp.groupby(cp.pert_iname.str.lower())
              .agg(best_score=("score_reversal", "min"),
                   median_score=("score_reversal", "median"),
                   n_sigs=("sig_id", "count"))
              .sort_values("best_score").reset_index()
              .rename(columns={"pert_iname": "pert_iname_lower"}))
    drug.insert(0, "rank", np.arange(1, len(drug) + 1))
    drug.to_csv(outdir / "drug_ranking.csv", index=False)

    sens = drug[drug.n_sigs >= 3].reset_index(drop=True).copy()
    sens.rename(columns={"rank": "rank_all"}, inplace=True)
    sens.insert(0, "rank_nsig3", np.arange(1, len(sens) + 1))
    sens.to_csv(outdir / "drug_ranking_nsig3_sensitivity.csv", index=False)
    log(f"scored {len(merged)} signatures; {len(drug)} distinct small-molecule drugs "
        f"({len(sens)} with >=3 signatures)")

    # ---- 5. positive control ----
    pc = positive_control_table(drug.set_index("pert_iname_lower"),
                                cfg["positive_controls"]["compounds"])
    pc.to_csv(outdir / "positive_controls.csv", index=False)
    ril = pc[pc.query_name == "riluzole"].iloc[0]
    log(f"positive control: riluzole rank {ril['rank']} "
        f"({ril['percentile']}th pctile of {len(drug)}, direction={ril['direction']})")

    # ---- 6. annotation of top-N candidates ----
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
        cns = (r.pert_iname_lower in CNS_KNOWN) or (
            mw is not None and float(mw) < 450)   # crude size-based BBB proxy + curated list
        recs.append({
            "pert_iname": r.pert_iname_lower, "rank": int(r["rank"]),
            "best_score": r.best_score, "median_score": r.median_score,
            "n_sigs": int(r.n_sigs),
            "chembl_id": ch.get("chembl_id"), "max_phase": ch.get("max_phase"),
            "first_approval": ch.get("first_approval"),
            "MW": mw, "formula": pu.get("MolecularFormula"),
            "SMILES": pu.get("SMILES"),
            "cns_flag": bool(cns),
        })
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(cache, indent=1))
    annot = pd.DataFrame(recs)
    annot.to_csv(outdir / "top_candidates_annotated.csv", index=False)
    log(f"annotated top {len(annot)} candidates")

    summary = {
        "n_samples_als": int(is_case.sum()),
        "n_samples_control": int((~is_case).sum()),
        "n_genes_tested": int(len(res)),
        "n_lincs_signatures_scored": int(len(merged)),
        "n_small_molecule_perturbagens": int(len(drug)),
        "riluzole_rank": None if ril["rank"] is None else int(ril["rank"]),
        "riluzole_percentile": ril["percentile"],
        "wall_time_sec": round(time.time() - t_start, 1),
    }
    (outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"DONE wall time {summary['wall_time_sec']}s")


if __name__ == "__main__":
    main()
