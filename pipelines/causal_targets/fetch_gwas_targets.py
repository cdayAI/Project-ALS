#!/usr/bin/env python
"""Stream A scaffolding: causal target triage (AGENTS.md task board).

Pipeline:
  1. Fetch every GWAS Catalog association for ALS and extract the unique
     gene list (author-reported genes per locus), with supporting evidence
     (association count, strongest p-value, PMIDs, GWAS Catalog accessions).
  2. Overlay each gene with Open Targets Platform tractability data (small
     molecule / antibody modality buckets) -- this is the "druggable genome"
     step: a gene with no tractability evidence in any bucket is not
     currently actionable by existing drug-discovery modalities, however
     causally implicated it is.
  3. Positive control (AGENTS.md rule 5): the known ALS GWAS genes cited in
     research/01_biology_and_therapeutics.md and research/03_data_resources.md
     (SOD1, C9orf72, UNC13A, TBK1, TARDBP, FUS) must appear in the output, or
     the pipeline has a bug, not a scientific finding.

IMPORTANT identifier note (found while building this, not assumed): the
research briefs cite ALS as GWAS Catalog trait "EFO_0000253". As of this
pipeline's build date, the GWAS Catalog REST API's internal trait record for
"amyotrophic lateral sclerosis" resolves under shortForm MONDO_0004976, not
EFO_0000253 -- querying /efoTraits/EFO_0000253 directly 404s. The website's
/gwas/studies/EFO_0000253 URL still works (separate lookup path), but this
pipeline's REST calls use MONDO_0004976, discovered via the API's own
findByEfoTrait text search. See reviews/causal_targets.md if it exists, or
research/03_data_resources.md, for the historical EFO_0000253 citation this
supersedes for API purposes.

Run:
  python fetch_gwas_targets.py --out outputs/
"""
import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

GWAS_API = "https://www.ebi.ac.uk/gwas/rest/api"
ALS_TRAIT_SHORTFORM = "MONDO_0004976"  # see module docstring
OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# Positive control: genes AGENTS.md rule 5 requires this pipeline to
# rediscover, restricted to genes established as common-variant GWAS hits
# (research/03_data_resources.md #1.6: "the big 2018-2022 ALS GWAS loci
# (C9orf72, UNC13A, TBK1, SOD1, etc.)"). Deliberately EXCLUDES TARDBP and FUS
# despite both being core ALS genes in research/01_biology_and_therapeutics.md
# -- a live run first included them and the positive control correctly
# FAILED, because TARDBP/FUS were discovered via family linkage/candidate-
# gene sequencing of RARE, highly-penetrant variants (research/01: TARDBP
# ~3-5% of familial ALS, <1% of sporadic; FUS ~0.3-0.9% of ALS), not via
# population case-control GWAS, which needs common alleles for power. Their
# absence here is a correct, expected result, not a pipeline bug -- and it's
# exactly why this pipeline's GWAS-only view is complementary to, not a
# replacement for, the rare-variant colocalization plan against Project MinE
# (COLOCALIZATION_PLAN.md).
POSITIVE_CONTROL_GENES = {"SOD1", "C9ORF72", "UNC13A", "TBK1"}


class FetchError(RuntimeError):
    """A request failed. Never silently treated as "zero results" -- same
    lesson as tools/trial_matcher (reviews/trial_matcher.md finding 1)."""


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def _get_json(url, params=None, timeout=30):
    if params:
        url = f"{url}?{urllib.parse.urlencode(params)}"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise FetchError(f"GET {url} failed: {exc}") from exc


def _post_graphql(query, variables=None, timeout=30):
    body = json.dumps({"query": query, "variables": variables or {}}).encode("utf-8")
    req = urllib.request.Request(
        OT_GRAPHQL, data=body, headers={"Content-Type": "application/json"}, method="POST"
    )
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise FetchError(f"POST {OT_GRAPHQL} failed: {exc}") from exc


# ----------------------------------------------------------------------------
# 1. GWAS Catalog: ALS associations -> gene-level evidence
# ----------------------------------------------------------------------------

def fetch_als_associations():
    """One request returns every ALS association (confirmed live: the
    `size` param on this projection is not honored -- it always returns the
    full set, 169 associations as of this pipeline's build). If that ever
    changes and results silently truncate, gene counts below would quietly
    undercount; there is no pagination cursor on this endpoint to detect
    truncation, so this is a known, documented limitation, not a bug."""
    url = f"{GWAS_API}/efoTraits/{ALS_TRAIT_SHORTFORM}/associations"
    data = _get_json(url, params={"projection": "associationByEfoTrait"})
    associations = data.get("_embedded", {}).get("associations", [])
    log(f"fetched {len(associations)} GWAS Catalog associations for ALS ({ALS_TRAIT_SHORTFORM})")
    return associations


def genes_from_associations(associations):
    """Collapse per-association loci into unique gene-level evidence."""
    genes = {}
    for assoc in associations:
        pmid = (assoc.get("study") or {}).get("publicationInfo", {}).get("pubmedId")
        accession = (assoc.get("study") or {}).get("accessionId")
        pvalue = assoc.get("pvalue")
        for locus in assoc.get("loci") or []:
            for g in locus.get("authorReportedGenes") or []:
                name = (g.get("geneName") or "").strip()
                if not name:
                    continue
                ensembl_ids = [e.get("ensemblGeneId") for e in g.get("ensemblGeneIds") or [] if e.get("ensemblGeneId")]
                entry = genes.setdefault(name, {
                    "gene": name,
                    "ensembl_ids": set(),
                    "n_associations": 0,
                    "min_pvalue": None,
                    "pmids": set(),
                    "gwas_accessions": set(),
                })
                entry["ensembl_ids"].update(ensembl_ids)
                entry["n_associations"] += 1
                if pvalue is not None and (entry["min_pvalue"] is None or pvalue < entry["min_pvalue"]):
                    entry["min_pvalue"] = pvalue
                if pmid:
                    entry["pmids"].add(str(pmid))
                if accession:
                    entry["gwas_accessions"].add(accession)
    return genes


# ----------------------------------------------------------------------------
# 2. Open Targets: tractability overlay ("druggable genome")
# ----------------------------------------------------------------------------

SEARCH_QUERY = """
query GeneSearch($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {size: 5, index: 0}) {
    hits { id name entity }
  }
}
"""

TARGET_QUERY = """
query TargetInfo($id: String!) {
  target(ensemblId: $id) {
    id
    approvedSymbol
    biotype
    tractability { label modality value }
  }
}
"""


def resolve_ensembl_id(gene_symbol, cached_ensembl_ids):
    """Prefer the Ensembl ID GWAS Catalog already gave us (exact, no search
    ambiguity); fall back to Open Targets' own search only if GWAS Catalog
    didn't supply one."""
    if cached_ensembl_ids:
        return sorted(cached_ensembl_ids)[0]
    result = _post_graphql(SEARCH_QUERY, {"q": gene_symbol})
    hits = ((result.get("data") or {}).get("search") or {}).get("hits") or []
    exact = [h for h in hits if h.get("name", "").upper() == gene_symbol.upper() and h.get("entity") == "target"]
    if exact:
        return exact[0]["id"]
    return hits[0]["id"] if hits else None


def tractability_summary(ensembl_id):
    if not ensembl_id:
        return {"resolved": False, "any_sm": False, "any_ab": False, "buckets": []}
    result = _post_graphql(TARGET_QUERY, {"id": ensembl_id})
    target = (result.get("data") or {}).get("target")
    if not target:
        return {"resolved": False, "any_sm": False, "any_ab": False, "buckets": []}
    buckets = [t for t in (target.get("tractability") or []) if t.get("value")]
    return {
        "resolved": True,
        "approved_symbol": target.get("approvedSymbol"),
        "any_sm": any(b["modality"] == "SM" for b in buckets),
        "any_ab": any(b["modality"] == "AB" for b in buckets),
        "buckets": [f"{b['modality']}:{b['label']}" for b in buckets],
    }


# ----------------------------------------------------------------------------
# 3. Orchestrate + positive control
# ----------------------------------------------------------------------------

def build_target_table(max_genes=None):
    associations = fetch_als_associations()
    genes = genes_from_associations(associations)
    log(f"{len(genes)} unique author-reported genes across all ALS GWAS loci")

    ordered = sorted(genes.values(), key=lambda g: (g["min_pvalue"] if g["min_pvalue"] is not None else 1.0))
    if max_genes:
        ordered = ordered[:max_genes]

    rows = []
    for i, g in enumerate(ordered):
        ensembl_id = resolve_ensembl_id(g["gene"], g["ensembl_ids"])
        try:
            tract = tractability_summary(ensembl_id)
        except FetchError as exc:
            log(f"  tractability lookup failed for {g['gene']}: {exc}")
            tract = {"resolved": False, "any_sm": False, "any_ab": False, "buckets": [], "error": str(exc)}
        rows.append({
            "gene": g["gene"],
            "ensembl_id": ensembl_id,
            "n_gwas_associations": g["n_associations"],
            "min_pvalue": g["min_pvalue"],
            "pmids": sorted(g["pmids"]),
            "gwas_accessions": sorted(g["gwas_accessions"]),
            "tractable_small_molecule": tract["any_sm"],
            "tractable_antibody": tract["any_ab"],
            "tractability_buckets": tract["buckets"],
        })
        if (i + 1) % 20 == 0:
            log(f"  tractability-annotated {i + 1}/{len(ordered)} genes")
    return rows


def check_positive_control(rows):
    found = {r["gene"].upper() for r in rows}
    missing = POSITIVE_CONTROL_GENES - found
    return {
        "expected": sorted(POSITIVE_CONTROL_GENES),
        "found": sorted(POSITIVE_CONTROL_GENES & found),
        "missing": sorted(missing),
        "passed": not missing,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="outputs", help="output directory")
    ap.add_argument("--max-genes", type=int, default=None, help="cap tractability lookups (rate-limit friendly); omit for full run")
    args = ap.parse_args(argv)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = build_target_table(max_genes=args.max_genes)
    control = check_positive_control(rows)

    (out_dir / "gwas_druggable_targets.json").write_text(json.dumps(rows, indent=2))

    csv_lines = ["gene,ensembl_id,n_gwas_associations,min_pvalue,tractable_small_molecule,tractable_antibody,pmids,gwas_accessions"]
    for r in rows:
        csv_lines.append(",".join([
            r["gene"], r["ensembl_id"] or "", str(r["n_gwas_associations"]),
            f"{r['min_pvalue']:.2e}" if r["min_pvalue"] is not None else "",
            str(r["tractable_small_molecule"]), str(r["tractable_antibody"]),
            "|".join(r["pmids"]), "|".join(r["gwas_accessions"]),
        ]))
    (out_dir / "gwas_druggable_targets.csv").write_text("\n".join(csv_lines))

    (out_dir / "positive_control.json").write_text(json.dumps(control, indent=2))

    log(f"positive control: {'PASSED' if control['passed'] else 'FAILED'} -- "
        f"found {control['found']}, missing {control['missing']}")
    log(f"wrote {len(rows)} genes to {out_dir}/gwas_druggable_targets.{{json,csv}}")

    if not control["passed"]:
        log("FAILING: positive control genes missing from GWAS Catalog output -- "
            "this pipeline has a bug (or trait-mapping regression), the run is NOT valid, "
            "do not use this output. See AGENTS.md rule 5.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
