# Review: pipelines/causal_targets (Stream A scaffolding)

Two-stage adversarial process per `AGENTS.md` rule 4 ("no self-certification"):
self-testing during development, then an independent review agent.

## Stage 1 — self-testing during development

A first version's positive control included TARDBP and FUS alongside
SOD1/C9orf72/UNC13A/TBK1, and correctly **FAILED** -- neither gene appeared
in the GWAS Catalog output at all. Investigated rather than silenced:
checked the actual per-gene GWAS output and cross-referenced against
`research/01_biology_and_therapeutics.md` (TARDBP mutations ~3-5% of
familial ALS, <1% sporadic; FUS ~0.3-0.9% of all ALS -- both discovered via
family-linkage/candidate-gene sequencing of rare, highly-penetrant
variants, not population case-control GWAS, which needs common allele
frequencies for statistical power). Concluded their absence is a correct,
expected result, not a bug, and removed them from the positive control with
that reasoning documented in code and README. This judgment call itself
went to Stage 2 for independent scrutiny rather than being trusted on the
builder's own say-so.

## Stage 2 — independent adversarial review

A separate agent, told to actually run things and read real data rather
than just review code, reported:

### TARDBP/FUS exclusion — verified independently, holds up
The reviewer queried GWAS Catalog directly for TARDBP (106 SNPs, 45
distinct traits) and FUS (186 SNPs, traits sampled) **without going through
this pipeline's code**. Both genes are well-indexed in GWAS Catalog
generally -- ruling out a broad indexing gap -- but zero of their
associated traits, across every SNP checked, is ALS or a synonym. This
confirms Stage 1's conclusion from an independent angle: not just "the ALS
query returns nothing" but "these genes' real GWAS hits are for other
diseases entirely." The positive-control judgment call stands.

### Real bug: GWAS Catalog placeholder tokens resolved to fabricated genes -- fixed
GWAS Catalog uses literal `"NR"` (no gene reported) and `"intergenic"` as
`authorReportedGenes` values for loci without a clean single-gene
assignment. The pipeline had no filter for these, so they fell through to
the Ensembl-ID-resolution fallback (a fuzzy Open Targets text search),
which resolved `"NR"` to the real oncogene **NRAS** and `"intergenic"` to a
random lncRNA (**GAPLINC**) -- then attached those genes' real
tractability data to rows that were not actually about NRAS or GAPLINC at
all. Concretely: the `"NR"` row had 16 associations at p=2e-14 (ranked #3
by significance in the output) and showed `tractable_small_molecule=True`
-- NRAS's real tractability, presented as if it were evidence about an
ALS-relevant target. The README already warned "NR"/"intergenic" aren't
real gene symbols, but never said their ensembl_id/tractability columns
were *also* fabricated -- a reader following the README's own advice would
still have trusted those columns.
**Fix:** added a `NON_GENE_TOKENS` filter in `genes_from_associations`,
applied before any Ensembl resolution is attempted, so no downstream code
path can look up tractability for a placeholder. Verified: re-running
after the fix drops from 81 to 79 genes (40 locus mentions skipped, logged
explicitly), and neither NRAS nor GAPLINC appear in the output. The
Ensembl-ID-resolution fallback itself was independently spot-checked
against 10 other genes with legitimate old HGNC aliases (KIAA1727,
MOBKL2B, JMJD2A, KIAA0196, C21orf29, BOD1L) and resolved all of them
correctly -- the bug was specific to the two non-gene placeholder tokens,
not the resolution logic generally.

### Real gap: tractability tier-collapsing overstated druggability -- fixed
The reviewer pulled live Open Targets data for TBK1, SOD1, and ACSL5 and
found `tractable_small_molecule=True` in all three cases came only from
the weakest evidence tiers ("Structure with Ligand", "High-Quality Ligand",
"Druggable Family" -- homology/structural-biology signal), with
"Approved Drug"/"Advanced Clinical"/"Phase 1 Clinical" all false. The
boolean was an `any()` across all 8 SM sub-buckets, so a reader could not
tell "an approved drug exists" from "this protein family is structurally
druggable in principle, no chemical matter yet" -- and the CSV (the more
commonly consumed artifact) only had the collapsed boolean, not the detail.
**Fix:** added `has_approved_drug_small_molecule`/`has_approved_drug_antibody`
as separate, stronger fields; renamed the original booleans to
`tractable_*_any_evidence` to make the weaker claim explicit in the field
name itself; added the full `tractability_buckets` list to the CSV (was
JSON-only before). Verified: TBK1 now shows
`tractable_small_molecule_any_evidence=true,
has_approved_drug_small_molecule=false`, with the buckets visible directly
in the CSV row.

### Minor: uncaught top-level fetch failure -- fixed
If the initial GWAS Catalog fetch itself failed (network down), `FetchError`
propagated uncaught out of `main()` as a raw traceback rather than a clean
message. Reviewer noted this does NOT repeat `tools/trial_matcher`'s
original bug class (reviews/trial_matcher.md finding 1: silently reporting
a failure as "zero results" with exit 0) -- exit code was still correctly
1 and nothing was silently treated as valid, so this was cosmetic, not a
correctness bug. **Fix:** wrapped the top-level `build_target_table()` call
in `main()` with the same try/except-and-clear-message pattern already
used for per-gene tractability failures.

### Confirmed clean (reviewer checked, no finding)
- Reproducibility: two independent live runs produced byte-identical
  output; no nondeterminism.
- Core extraction logic (`genes_from_associations`, min-p-value selection):
  independently reimplemented from the raw API response in a separate
  script by the reviewer; zero mismatches across all 81 genes from the
  pre-fix run.
- AGENTS.md rule 1 (no raw data committed): output files are 36KB/8KB/180B,
  nowhere near the 50MB threshold.
- AGENTS.md rule 6 (cite or it didn't happen): all data-source claims in
  README.md verified live by the reviewer (EFO_0000253 404s, MONDO_0004976
  200s, `size` param genuinely ignored, ChEMBL genuinely down).

## Verdict

**SURVIVES**, after fixing one real correctness bug (fabricated gene
identities from placeholder tokens) and one real honesty gap (tier-collapsed
druggability claims) that self-testing missed. The TARDBP/FUS
positive-control judgment call -- the one place Stage 1 made a genuinely
debatable call rather than following a clear rule -- was independently
re-derived from primary data by Stage 2, not just accepted on the builder's
word. Per README operating rule 3 ("a method that can't rediscover known
truth isn't trusted to find new truth"): this pipeline's own positive
control failing on a first attempt, being investigated rather than
patched-until-green, and then holding up under independent re-derivation,
is the process working as intended, not a story to sand down before
merging.
