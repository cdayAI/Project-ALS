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
`--mutation` (`sod1`, `c9orf72`, `tardbp`, `fus`, `atxn2`, `sporadic`),
`--stage` (`early`, `late`), `--status` (default `RECRUITING`), `--limit`,
`--strict`, `--json`. No API key, no dependencies beyond the Python 3
standard library.

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

- **It does not determine eligibility.** It ranks trials by keyword
  proximity in free-text `EligibilityCriteria`, which cannot reliably parse
  negation or scope. See "Known limitations" below for two real, confirmed
  failure cases found during validation.
- It does not check drug interactions, insurance, travel feasibility, or
  anything beyond what ClinicalTrials.gov publishes.
- It is not a substitute for talking to a neurologist/ALS clinic or the
  trial site's own contact (listed at each result's URL).

## Known limitations (found during validation, not theoretical)

Mutation/stage matching is plain keyword + negation-window heuristics over
free text (see `_find_with_negation`, `split_eligibility` in
`trial_matcher.py`), not real NLP. Three real trials exposed exactly why
this is a floor, not a solved problem -- full detail in
`reviews/trial_matcher.md`:

1. **False positive** (fixed): a trial listed "Confirmed mutation in the
   SOD1, FUS or C9orf72 gene" under *Exclusion* Criteria; naive keyword
   matching scored it as a good SOD1 match. Fixed by splitting eligibility
   text on the "Exclusion Criteria" heading and scoring exclusion-section
   hits as strongly negative.
2. **False positive, different shape** (fixed): a trial titled "Tofersen in
   **Non-SOD1** ALS" required inclusion criteria "testing **negative for**
   SOD1" -- negation *inside* the Inclusion section, invisible to a
   section-based split alone. Fixed with a negation-cue window check
   around each keyword match.
3. **False negative** (documented, not fixed -- fundamentally hard): a
   trial explicitly *for* SOD1-ALS patients scored strongly *negative*
   because its exclusion criteria mention "SOD1-ALS" while excluding
   *concurrent enrollment in another trial*, not excluding the mutation.
   No keyword-window heuristic distinguishes "excluding people who have
   X and are also doing Y" from "excluding people who have X". This is a
   real semantic-parsing problem; `--strict` mode can silently hide a
   genuinely strong match because of it. **Do not trust `--strict` mode to
   be complete. Always read the full ranked list, non-strict, or the
   criteria text yourself.**

Given (1)-(3), the scores and notes are a reading-priority aid, not a
filter. This is why the CLI prints a disclaimer to stderr (and embeds one
in `--json` output) on every run, and why `--strict` carries an explicit
danger warning in `--help`.

## Positive control

Riluzole and tofersen are the two FDA-approved ALS drugs with the clearest
public trial footprint (tofersen post-marketing/safety studies are
explicitly SOD1-targeted). A working query for `--mutation sod1` should
surface a tofersen-related SOD1 study somewhere in its (non-strict)
results -- confirmed manually during validation (NCT07259980). If a future
change to this tool makes that stop being true, treat it as a regression.
