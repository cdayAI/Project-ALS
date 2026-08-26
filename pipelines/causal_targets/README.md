# Causal target triage (Stream A scaffolding)

Fetches ALS GWAS Catalog loci, extracts the gene list, and overlays
druggability (tractability) data. Built for `AGENTS.md`'s Stream A. Unblocks
hypotheses H-003/H-004 in the hypothesis ledger.

## Usage

```
cd pipelines/causal_targets
python fetch_gwas_targets.py --out outputs
```

Takes ~25s live (169 GWAS associations, 79 unique genes after filtering
non-gene placeholder tokens, one Open Targets tractability lookup per
gene). No API key. Exit code 1 if the positive control fails, or if the
GWAS Catalog / Open Targets fetch itself fails outright (a network/API
failure is never reported as "zero genes found" -- see Known limitations
in `reviews/causal_targets.md`) -- either way, treat exit 1 as "do not use
this run's output," not a warning to ignore.

## Data sources / reproducibility contract (AGENTS.md rule 6)

- **GWAS Catalog REST API**: `GET https://www.ebi.ac.uk/gwas/rest/api/efoTraits/MONDO_0004976/associations?projection=associationByEfoTrait`.
  **Identifier note**: the research briefs (`research/01`, `research/03`)
  cite ALS as `EFO_0000253`. As of this pipeline's build, the GWAS Catalog
  REST API's internal record for ALS resolves under shortForm
  `MONDO_0004976` -- `GET .../efoTraits/EFO_0000253` 404s. Found live via
  the API's own `/efoTraits/search/findByEfoTrait?trait=amyotrophic+lateral+sclerosis`
  text search, not assumed. The public website URL
  `https://www.ebi.ac.uk/gwas/studies/EFO_0000253` still resolves (a
  separate lookup path from the REST API), so this is not a claim that
  EFO_0000253 is wrong -- only that the REST API needs the MONDO id.
  Confirmed live: 169 associations, `size` query param not honored by this
  endpoint (always returns the full set in one response; verified by
  comparing `size=1` and `size=2000` -- both returned 169).
- **Open Targets Platform GraphQL API**: `POST https://api.platform.opentargets.org/api/v4/graphql`,
  `target(ensemblId: ...) { tractability { label modality value } }`.
  Ensembl gene IDs are taken directly from GWAS Catalog's own
  `authorReportedGenes.ensemblGeneIds` where available (exact match, no
  search ambiguity); only genes GWAS Catalog didn't supply an Ensembl ID
  for fall back to an Open Targets text search.
- **Why not ChEMBL directly**: ChEMBL's public REST API
  (`www.ebi.ac.uk/chembl/api/data/...`) returned HTTP 500 on every endpoint
  tested during this build, including its own status endpoint -- a live,
  server-side EBI outage, not a URL/query error (confirmed:
  `research/04_pilot_access_report.md` shows the same ChEMBL endpoint
  returning 200 on an earlier date). Open Targets' tractability data is
  itself partly built from ChEMBL plus SureChEMBL/PDB/other sources, so it
  is not a downgrade -- it is arguably the more direct source for
  "druggable genome" specifically (tractability buckets), vs. ChEMBL target
  search which only confirms a target record exists.

## Positive control (AGENTS.md rule 5)

`SOD1, C9ORF72, UNC13A, TBK1` -- established common-variant ALS GWAS loci
per `research/03_data_resources.md` §1.6. **Deliberately excludes TARDBP
and FUS**: a first version of this pipeline included them, the positive
control correctly FAILED, and investigation confirmed their absence from
GWAS Catalog is real (they're rare-variant/familial-linkage genes, not
GWAS hits) -- not a pipeline bug. See `COLOCALIZATION_PLAN.md` for why that
finding is itself the reason a separate rare-variant analysis axis exists.

## Known data-quality notes (from the live GWAS Catalog data itself, not this code)

- Some loci have `authorReportedGenes` of literal `"NR"` (not reported by
  the study authors), `"intergenic"`, or similar placeholder tokens (see
  `NON_GENE_TOKENS` in `fetch_gwas_targets.py`) -- **these are filtered out
  before gene extraction**, not passed through. An earlier version let them
  fall through to Open Targets' fuzzy search, which resolved "NR" to the
  real oncogene NRAS and attached NRAS's real tractability data to a
  16-association, p=2e-14 row that was actually about no specific gene at
  all -- caught by independent review, see `reviews/causal_targets.md`. The
  pipeline logs how many locus mentions were skipped this way.
- `n_gwas_associations` counts independent GWAS Catalog association
  records mentioning a gene, not independent discovery cohorts -- a gene
  hit by the same locus in multiple studies of overlapping samples will
  still count once per study/record.
- **`tractable_*_any_evidence=True` does NOT mean a drug exists.** Open
  Targets reports 8 small-molecule and 10 antibody evidence tiers from
  "Approved Drug" down to weak structural-homology signal ("Structure with
  Ligand", "Druggable Family"). `*_any_evidence` is true if ANY tier is
  true -- verified live that TBK1/SOD1/ACSL5 all show
  `tractable_small_molecule_any_evidence=True` from the weakest tiers only,
  with no approved or clinical-stage compound. Use
  `has_approved_drug_small_molecule`/`has_approved_drug_antibody` for that
  specific, much stronger claim, or read `tractability_buckets` (now in
  both the JSON and the CSV, not JSON-only) for the real tier breakdown.

## Outputs

- `outputs/gwas_druggable_targets.{json,csv}` -- one row per gene: GWAS
  evidence (association count, min p-value, PMIDs, GWAS Catalog
  accessions) + tractability (`tractable_*_any_evidence`,
  `has_approved_drug_*`, and the full `tractability_buckets` list -- see
  the tier-honesty note above before trusting the any_evidence booleans
  alone).
- `outputs/positive_control.json` -- pass/fail + which control genes were
  found/missing on this run.
