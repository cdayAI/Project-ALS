# Trial matcher (Stream C)

Ranks currently-recruiting ALS trials from ClinicalTrials.gov by country/region,
known causal mutation, and rough disease stage. Built for `AGENTS.md`'s Stream C
("fastest human impact" per `research/00_synthesis.md`).

## Usage

```
python trial_matcher.py --country "United States" --mutation c9orf72
python trial_matcher.py --country France --stage late
python trial_matcher.py --mutation sod1 --json > results.json
```

Flags: `--condition` (default "Amyotrophic Lateral Sclerosis"), `--country`,
`--mutation` (common: `sod1`, `c9orf72`, `tardbp`, `fus`, `atxn2`, `tbk1`,
`optn`, `vcp`, `ubqln2`, `ang`, `pfn1`, `nek1`, `matr3`, `sporadic` -- but
any gene name works, matched verbatim if not in that list), `--stage`
(common: `early`, `late` -- any free-text cue also works), `--status`
(default `RECRUITING`), `--limit`, `--json`. No API key, no dependencies
beyond the Python 3 standard library. Exit code is 1 if ClinicalTrials.gov
could not be reached at all (a real network/API failure); 0 with a
possibly-empty result list otherwise.

## Data source / reproducibility contract (AGENTS.md rule 6)

- Endpoint: `GET https://clinicaltrials.gov/api/v2/studies` (NIH/NLM
  ClinicalTrials.gov API v2, public, no key required).
- Params used: `query.cond`, `query.locn`, `filter.overallStatus`,
  `pageSize`, `pageToken`, `fields`.
- Fields fetched: see `FIELDS` in `trial_matcher.py` (NCTId, BriefTitle,
  OverallStatus, Phase, StudyType, LeadSponsorName, EnrollmentCount,
  StartDate, PrimaryCompletionDate, EligibilityCriteria, MinimumAge,
  MaximumAge, LocationCountry/City/Facility).
- No local data is stored; every run hits the live API.

## What this tool does NOT do

- **It does not determine eligibility.** It ranks trials using public
  metadata and keyword proximity in free-text `EligibilityCriteria`, which
  cannot reliably parse negation or clause scope. See "Known limitations."
- `--country` filtering trusts ClinicalTrials.gov's `query.locn` text
  match, which can match a US state name against a same-named country (see
  finding 3 below); the tool cross-checks each result's actual location
  list before claiming a country match, and widens to a global search with
  an explicit note if nothing genuinely matches.
- It does not check drug interactions, insurance, travel feasibility, or
  anything beyond what ClinicalTrials.gov publishes.
- It is not a substitute for talking to a neurologist/ALS clinic or the
  trial site's own contact (listed at each result's URL).

## Known limitations

This tool went through two rounds of adversarial testing: self-testing
during development, then an independent review (per `AGENTS.md` rule 4 --
"no self-certification"). The independent review found the self-tested
version's own positive-control trial (see below) was itself badly
mis-scored, which is why mutation/stage matching no longer emits a
confident "this trial EXCLUDES you" verdict. Full findings and disposition
are in `reviews/trial_matcher.md`. Summary of what's fixed vs. still true:

**Fixed:**
1. A trial listed "Confirmed mutation in the SOD1, FUS or C9orf72 gene"
   under *Exclusion* Criteria; naive whole-text keyword matching scored it
   as a good SOD1 match.
2. A trial titled "Tofersen in Non-SOD1 ALS" required inclusion criteria
   "testing negative for SOD1" -- negation *inside* the Inclusion section.
3. A country filter for "Georgia" (the country) returned only trials
   located in the US state of Georgia, with no indication of the mismatch.
4. `--mutation`/`--stage` were restricted to a fixed choice list, hard-
   blocking real ALS genes (TBK1, OPTN, VCP, ...) not on it.
5. A ClinicalTrials.gov network/API failure was indistinguishable from "0
   trials found" -- meaning an outage could be reported as "no ALS trials
   exist anywhere". Now raises a distinct error with a non-zero exit code.

**Still true, not fixable by keyword heuristics (fundamental, not a bug):**
Keyword + negation-window matching cannot parse *why* a mutation is
mentioned near "exclusion". Two independently confirmed shapes:
- A trial explicitly *for* SOD1-ALS patients (a tofersen safety study)
  mentions "SOD1-ALS" in its exclusion criteria only to exclude
  *concurrent enrollment in another trial* -- not the mutation itself.
- A trial excluding SOD1/FUS carriers explains why with "known absence of
  TDP-43 pathology" -- for a `--mutation tardbp` search this reads as a
  TARDBP exclusion signal when it's actually irrelevant to TARDBP
  patients.
- A trial with an explicit carve-out ("...except for TARDBP gene
  variants") inside otherwise-exclusionary text needs a human to read the
  actual clause to see the carve-out.

**Because of this, mutation and stage notes now quote the actual source
sentence instead of asserting a verdict, and scores use small +-1 weights
for text-derived signals (vs. larger weights for structured data like
country/phase) so a wrong call nudges rank instead of burying a genuinely
strong match.** There is deliberately no hide-non-matches / `--strict`
mode: an earlier version had one, and the independent review proved it
could silently hide the single best-matching trial for a real scenario.

## Positive control

Tofersen is the FDA-approved SOD1-targeted ASO with the clearest public
trial footprint. `--mutation sod1` reliably surfaces a tofersen/SOD1 study
(NCT07259980) in its results -- but note its *score* is intentionally
unremarkable (not high) precisely because of the fundamental limitation
above: this is the trial that exposed the problem, not one the scoring
gets confidently right. The real check is "does it appear, with its
mutation note quoting the actual sentence for a human to read" -- not
"does it score highest". If a future change makes it stop appearing at
all, or makes it score a confident negative, treat that as a regression.
