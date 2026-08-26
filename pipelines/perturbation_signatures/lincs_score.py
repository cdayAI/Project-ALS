#!/usr/bin/env python
"""LINCS L1000 Level5 streaming reversal scorer for pipelines/perturbation_signatures.

Core scoring functions vendored from sprint1-repurposing's salvaged exp001 code
(branch exp001-sprint1-handoff, commit fcfff64, file experiments/exp001_repurposing/run.py,
MIT-style repo license), with attribution per board coordination. Adaptations: multi-query
batch scoring (single pass for Q columns), gene-id alignment helper.
"""
import time
import numpy as np


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


def stream_scores_multi(gctx_path, query_mat, log=print, sig_block=10000, progress_every=10):
    """Cosine similarity between each column of query_mat (n_genes x n_queries, aligned
    to GCTX gene ROW order) and every signature. ONE sequential pass over signature-row
    blocks regardless of the number of queries.

    Returns scores ndarray float32 of shape (n_signatures, n_queries), value = -cosine
    (higher = stronger REVERSAL of the query signature).
    """
    import h5py
    with h5py.File(gctx_path, "r") as f:
        mat = _descend_to_matrix(f)
        s, g = mat.shape
        Q = query_mat.astype(np.float64)
        assert g == Q.shape[0], (mat.shape, Q.shape)
        qn = np.linalg.norm(Q, axis=0)                      # per-query norms
        dots = np.empty((s, Q.shape[1]), dtype=np.float64)
        norms = np.empty(s, dtype=np.float64)
        t0 = time.time()
        for a in range(0, s, sig_block):
            b = min(a + sig_block, s)
            Y = mat[a:b, :].astype(np.float64)              # sig_block x genes
            dots[a:b] = Y @ Q
            norms[a:b] = np.einsum("ij,ij->i", Y, Y)
            done = b / s
            el = time.time() - t0
            if (a // sig_block) % progress_every == 0 or b == s:
                eta = el / done - el if done > 0 else float("nan")
                log(f"scoring {b}/{s} signatures ({done:.1%}), elapsed {el:.0f}s, ETA {eta:.0f}s")
    denom = np.sqrt(np.maximum(norms, 1e-12))[:, None] * qn[None, :]
    scores = -(dots / np.maximum(denom, 1e-300))
    return scores.astype(np.float32)


def align_query_to_gctx(symbol_tstat_pairs, gene_ids, pr_id_by_symbol):
    """Build a dense query vector in GCTX ROW order.

    symbol_tstat_pairs: iterable of (hgnc_symbol, signed_statistic)
    gene_ids: GCTX /0/META/ROW/id values (pr_gene_id strings for GSE92742)
    pr_id_by_symbol: dict hgnc_symbol -> pr_gene_id
    Returns (vector float64, n_mapped, n_total)
    """
    idx_of_gid = {gid: i for i, gid in enumerate(gene_ids)}
    vec = np.zeros(len(gene_ids), dtype=np.float64)
    mapped = 0
    for sym, stat in symbol_tstat_pairs:
        gid = pr_id_by_symbol.get(sym)
        if gid is not None and gid in idx_of_gid:
            vec[idx_of_gid[gid]] = stat
            mapped += 1
    return vec, mapped, len(list(symbol_tstat_pairs))
