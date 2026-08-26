#!/usr/bin/env python
"""te_contamination_check.py - pre-flight QC for raw count matrices.

MOTIVATION
GEO RNA-seq supplementary count files sometimes mix gene rows with non-gene
feature classes - transposable elements / repeat-consensus rows are the worst
offender seen so far: in GSE124439 (ALS post-mortem CNS), ~1k TE rows carry
~26% of all raw counts. If those rows are present during CPM/TPM/library-size
normalization, EVERY real gene signal is diluted by a global factor and the
downstream analysis silently loses power (this killed the first exp001
signature attempt before diagnosis).

WHAT THIS TOOL DOES
1. Loads a count matrix either from a directory of GEO-style per-sample files
   (`*_counts.txt[.gz]`, tab-separated, optional leading metadata line) or from
   a single matrix file (tsv/csv, feature IDs in the first column).
2. Classifies every feature ID against an optional gene-reference whitelist
   (one symbol per line; e.g. NCBI `Homo_sapiens.gene_info` Symbol column,
   HGNC complete set). Non-whitelisted IDs are subclassified by pattern rules:
     - IDs containing ':'           -> transposable_element/repeat consensus
     - anything else not whitelisted -> other_non_gene
   Without a whitelist only the pattern heuristics run (mode is flagged).
3. Quantifies each feature class's share of total raw counts (per-sample and
   aggregate). This is the number that matters: a few thousand near-zero TE
   rows are harmless, but TE rows carrying >1% of library reads corrupt
   size-factor normalization.
4. Emits a GO / CAUTION / NO-GO normalization recommendation.

VERDICT RULES (non-gene share of total raw reads)
    >= 5%   NO-GO   - strip non-gene rows BEFORE any normalization; do not
                      trust analyses normalized on the raw matrix
    >= 1%   CAUTION - strip rows and re-run; flag affected samples
    <  1%   GO      - proceed (still record the finding)

USAGE
    python te_contamination_check.py \
        --input-dir data/gse124439 --pattern '*_counts.txt.gz' \
        --gene-whitelist data/ref/Homo_sapiens.gene_info_symbols.txt \
        --outdir pipelines/data_qc/demo_output

    python te_contamination_check.py --matrix counts.tsv --gene-whitelist symbols.txt

Exit codes: 0 = GO/CAUTION finished, 2 = NO-GO (only with --fail-on-no-go),
1 = operational error.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

TE_PATTERN = ":"            # repeat consensus names, e.g. "X8_LINE:CR1:LINE", "Zaphod:hAT-Tip100:DNA"
GENE_LIKE_RE = re.compile(r"^[A-Za-z][A-Za-z0-9._@|-]*$")

CLASS_TE = "transposable_element_repeat"
CLASS_OTHER = "other_non_gene"
CLASS_GENE = "gene"
CLASS_GENE_UNVERIFIED = "gene_like_unverified"


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------- #
# input loaders
# --------------------------------------------------------------------------- #

def load_counts_from_dir(indir: Path, pattern: str, max_samples: int | None):
    """GEO-style directory: one `<acc>_<title>_counts.txt[.gz]` per sample."""
    paths = sorted(Path(indir).glob(pattern))
    if not paths:
        raise FileNotFoundError(f"no files match {pattern!r} under {indir}")
    if max_samples:
        paths = paths[:max_samples]
    series = {}
    for p in paths:
        opener = gzip.open if p.suffix == ".gz" else open
        with opener(p, "rt") as fh:
            first = fh.readline()
            # GEO files start with a metadata line like "gene/TE\t<bam path>"
            skip = 0 if "\t" not in first or not first.split("\t")[1][:2].startswith(("..", "/", "~")) else 1
        df = pd.read_csv(p, sep="\t", skiprows=skip, header=None,
                         names=["fid", p.stem], dtype={0: str})
        df["fid"] = df["fid"].str.strip('"')
        series[p.stem] = df.set_index("fid")[p.stem]
    mat = pd.DataFrame(series).apply(pd.to_numeric, errors="coerce").fillna(0)
    if mat.values.sum() <= 0:
        raise ValueError(f"{indir}: parsed count files contain zero raw counts - "
                         f"check parsing assumptions (skiprows/delimiter)")
    return mat


def load_matrix_file(path: Path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        first = fh.readline()
    sep = "\t" if "\t" in first else ","
    mat = pd.read_csv(path, sep=sep, index_col=0)
    if mat.shape[1] == 0 or first.count(sep) < 1:
        raise ValueError(f"{path}: could not parse columns with detected delimiter "
                         f"{sep!r} - check file format")
    mat = mat.apply(pd.to_numeric, errors="coerce").fillna(0)
    if mat.values.sum() <= 0:
        raise ValueError(f"{path}: matrix contains zero raw counts - wrong "
                         f"delimiter/format would silently produce GO verdicts")
    return mat


def load_whitelist(path: Path) -> set:
    """Accepts a plain one-symbol-per-line file OR an NCBI `gene_info`-style TSV.

    For gene_info-style input, Symbol, Symbol_from_nomenclature_authority and
    Synonyms are all added - older count matrices often use outdated symbols
    (e.g. SEPT7 instead of SEPTIN7) and they must not be misclassified as
    non-gene features.
    """
    syms: set = set()
    with open(path) as fh:
        cols = fh.readline().rstrip("\n").split("\t")
        sym_idx = cols.index("Symbol") + 1 if "Symbol" in cols else None
        extra_idx = [cols.index(c) + 1 for c in ("Synonyms", "Aliases",
                                                 "Symbol_from_nomenclature_authority")
                     if c in cols]
        fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")

            def add(s):
                s = s.strip().strip('"').upper()
                if s:
                    syms.add(s)

            if sym_idx:
                add(parts[sym_idx - 1])
                for ix in extra_idx:
                    if len(parts) >= ix:
                        for tok in parts[ix - 1].split("|"):
                            if tok and tok != "-":
                                add(tok)
            else:
                add(parts[0])
    return syms


# --------------------------------------------------------------------------- #
# classification + metrics
# --------------------------------------------------------------------------- #

def classify_features(fids: pd.Index, whitelist: set | None):
    """Return dict class -> set of feature ids."""
    out: dict[str, set] = {}
    upper = {f: str(f).upper() for f in fids}

    def bucket(fid: str, cls: str):
        out.setdefault(cls, set()).add(fid)

    for f in fids:
        s = str(f)
        if whitelist is not None:
            if upper[f] in whitelist:
                bucket(f, CLASS_GENE)
            elif TE_PATTERN in s:
                bucket(f, CLASS_TE)
            else:
                bucket(f, CLASS_OTHER)
        else:  # heuristic-only mode
            if TE_PATTERN in s:
                bucket(f, CLASS_TE)
            elif GENE_LIKE_RE.match(s):
                bucket(f, CLASS_GENE_UNVERIFIED)
            else:
                bucket(f, CLASS_OTHER)
    return out


def class_metrics(mat: pd.DataFrame, classes: dict[str, set]):
    total = float(mat.values.sum())
    rows = []
    for cls, ids in sorted(classes.items()):
        sub = mat.loc[list(ids)]
        per_sample = sub.sum(axis=0) / max(total, 1e-12)
        rows.append({
            "class": cls,
            "n_features": int(len(ids)),
            "frac_features": round(len(ids) / len(mat), 5),
            "frac_raw_reads": round(float(sub.values.sum()) / max(total, 1e-12), 5),
            "per_sample_read_frac_min": round(float(per_sample.min()), 5),
            "per_sample_read_frac_median": round(float(per_sample.median()), 5),
            "per_sample_read_frac_max": round(float(per_sample.max()), 5),
            "top_features_by_reads": [
                {"feature": str(k), "reads": int(v)}
                for k, v in sub.sum(axis=1).sort_values(ascending=False)
                .head(10).items()
            ],
        })
    return rows, total


def make_verdict(metrics_rows, used_whitelist: bool):
    non_gene_frac = sum(
        r["frac_raw_reads"] for r in metrics_rows if r["class"] != CLASS_GENE)
    if non_gene_frac >= 0.05:
        verdict, code = "NO-GO", 2
        rec = ("Non-gene features carry {:.1%} of raw reads. Strip them BEFORE "
               "CPM/TPM/size-factor normalization; any analysis normalized on "
               "the raw matrix is globally diluted and should be re-run."
               ).format(non_gene_frac)
    elif non_gene_frac >= 0.01:
        verdict, code = "CAUTION", 0
        rec = ("Non-gene features carry {:.1%} of raw reads (>1%). Remove them "
               "before normalization and note affected samples."
               ).format(non_gene_frac)
    else:
        verdict, code = "GO", 0
        rec = "Non-gene features carry only {:.2%} of raw reads; safe to normalize.".format(non_gene_frac)
    return {
        "verdict": verdict,
        "exit_code": code,
        "non_gene_read_fraction": round(non_gene_frac, 5),
        "whitelist_used": used_whitelist,
        "recommendation": rec,
    }


# --------------------------------------------------------------------------- #
# reporting
# --------------------------------------------------------------------------- #

def write_report_md(res: dict, path: Path):
    L = ["# Data QC: non-gene feature contamination check",
         "",
         f"- Generated: {res['meta']['generated']}",
         f"- Input: `{res['meta']['input']}`",
         f"- Features analysed: {res['meta']['n_features']:,} x {res['meta']['n_samples']} samples",
         f"- Gene reference whitelist: **{'YES' if res['verdict']['whitelist_used'] else 'NO (heuristic mode - results less specific)'}**",
         "",
         f"## Verdict: **{res['verdict']['verdict']}**",
         "",
         res["verdict"]["recommendation"],
         "",
         "| class | n features | % features | % raw reads | per-sample read % (min/med/max) |",
         "|---|---|---|---|---|"]
    for r in res["classes"]:
        L.append("| {} | {:,} | {:.2%} | {:.2%} | {:.2%} / {:.2%} / {:.2%} |".format(
            r["class"], r["n_features"], r["frac_features"], r["frac_raw_reads"],
            r["per_sample_read_frac_min"], r["per_sample_read_frac_median"],
            r["per_sample_read_frac_max"]))
    L += ["", "## Top non-gene features by raw reads", "",
          "| feature | class | reads |", "|---|---|---|"]
    top = sorted(((f, r["class"], d["reads"])
                  for r in res["classes"] if r["class"] != CLASS_GENE
                  for f, d in ((x["feature"], x) for x in r["top_features_by_reads"])),
                 key=lambda t: -t[2])[:20]
    for f, cls, n in top:
        L.append(f"| {f} | {cls} | {n:,} |")
    path.write_text("\n".join(L) + "\n")


# --------------------------------------------------------------------------- #

def run_check(matrix: pd.DataFrame, whitelist: set | None, meta: dict):
    classes = classify_features(matrix.index, whitelist)
    rows, _total = class_metrics(matrix, classes)
    verdict = make_verdict(rows, whitelist is not None)
    return {
        "meta": {**meta,
                 "n_features": int(matrix.shape[0]),
                 "n_samples": int(matrix.shape[1]),
                 "generated": time.strftime("%Y-%m-%d %H:%M:%S")},
        "classes": rows,
        "verdict": verdict,
    }


def self_test():
    """Synthetic sanity checks of classifier + verdict logic."""
    mat = pd.DataFrame(
        [[100, 90], [50, 60], [2000, 2200], [5, 4]],
        index=["A1BG", "X8_LINE:CR1:LINE", "Zaphod:hAT-Tip100:DNA", "NOT_A_SYMBOL!"],
        columns=["s1", "s2"])
    wl = {"A1BG"}
    res = run_check(mat, wl, {"input": "synthetic"})
    assert res["verdict"]["verdict"] == "NO-GO", res["verdict"]
    cls = {r["class"]: r for r in res["classes"]}
    assert cls[CLASS_GENE]["n_features"] == 1
    assert cls[CLASS_TE]["n_features"] == 2
    assert cls[CLASS_OTHER]["n_features"] == 1
    assert abs(sum(r["frac_raw_reads"] for r in res["classes"]) - 1.0) < 1e-3
    # clean matrix -> GO
    res2 = run_check(mat.loc[["A1BG"]], wl, {"input": "synthetic"})
    assert res2["verdict"]["verdict"] == "GO"
    print("self-test PASSED")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("MOTIVATION")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    effective_argv = argv if argv is not None else sys.argv[1:]
    self_test_requested = "--self-test" in effective_argv
    src = ap.add_mutually_exclusive_group(required=not self_test_requested)
    src.add_argument("--input-dir", help="directory of per-sample count files")
    src.add_argument("--matrix", help="single count-matrix tsv/csv (feature ids in col 0)")
    ap.add_argument("--pattern", default="*_counts.txt.gz",
                    help="glob for --input-dir (default %(default)s)")
    ap.add_argument("--max-samples", type=int, default=None,
                    help="analyse at most N samples (speed)")
    ap.add_argument("--gene-whitelist", default=None,
                    help="one symbol per line, or NCBI gene_info TSV with Symbol column")
    ap.add_argument("--outdir", default="qc_output", help="where to write report.json/md")
    ap.add_argument("--fail-on-no-go", action="store_true",
                    help="exit code 2 when verdict is NO-GO (CI-friendly)")
    ap.add_argument("--self-test", action="store_true", help="run synthetic checks and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        self_test()
        return 0

    wl = load_whitelist(Path(args.gene_whitelist)) if args.gene_whitelist else None
    log(f"whitelist: {len(wl):,} symbols" if wl else "whitelist: NONE - heuristic mode")

    if args.input_dir:
        log(f"loading per-sample files from {args.input_dir}")
        mat = load_counts_from_dir(Path(args.input_dir), args.pattern, args.max_samples)
        meta = {"input": str(args.input_dir)}
    else:
        log(f"loading matrix {args.matrix}")
        mat = load_matrix_file(Path(args.matrix))
        meta = {"input": str(args.matrix)}

    log(f"{mat.shape[0]:,} features x {mat.shape[1]} samples")
    res = run_check(mat, wl, meta)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "qc_report.json").write_text(json.dumps(res, indent=1))
    write_report_md(res, outdir / "qc_report.md")
    log(f"verdict: {res['verdict']['verdict']} (non-gene read fraction "
        f"{res['verdict']['non_gene_read_fraction']:.2%}) -> {outdir/'qc_report.md'}")
    print(res["verdict"]["recommendation"])

    if args.fail_on_no_go and res["verdict"]["exit_code"] == 2:
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
