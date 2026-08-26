# Causal target triage (Stream A scaffolding)

Fetches ALS GWAS Catalog loci, extracts the gene list, and overlays
druggability (tractability) data. Built for `AGENTS.md`'s Stream A. Unblocks
hypotheses H-003/H-004 in the hypothesis ledger.

## Usage

```
cd pipelines/causal_targets
python fetch_gwas_targets.py --out outputs
```

Takes ~25s live (169 GWAS associations, 81 unique genes, one Open Targets
tractability lookup per gene). No API key. Exit code 1 if the positive
control fails (see below) -- treat that as "do not trust this run's output,"
not a warning to ignore.

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
  the study authors) or `"intergenic"` -- these pass through to the output
  as genes named "NR"/"intergenic" rather than being silently dropped, so
  the row counts are honest, but don't treat those two labels as real gene
  symbols when reading `gwas_druggable_targets.csv`.
- `n_gwas_associations` counts independent GWAS Catalog association
  records mentioning a gene, not independent discovery cohorts -- a gene
  hit by the same locus in multiple studies of overlapping samples will
  still count once per study/record.

## Outputs

- `outputs/gwas_druggable_targets.{json,csv}` -- one row per gene: GWAS
  evidence (association count, min p-value, PMIDs, GWAS Catalog
  accessions) + tractability (small-molecule / antibody bucket booleans
  and the specific buckets that were true).
- `outputs/positive_control.json` -- pass/fail + which control genes were
  found/missing on this run.
