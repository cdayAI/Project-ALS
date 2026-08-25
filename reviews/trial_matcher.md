# Review: tools/trial_matcher (Stream C)

Two-stage adversarial process per `AGENTS.md` rule 4 ("no self-certification"):
self-testing during development, then an independent review agent with no
visibility into the builder's reasoning, tasked with finding problems beyond
what was already documented.

## Stage 1 — self-testing during development

Found and fixed two real false-positive matching bugs by running the tool
against the live API and manually reading the eligibility text of top-ranked
results:

1. NCT07322003 ("Pridopidine Phase 3... ALS") lists "Confirmed mutation in
   the SOD1, FUS or C9orf72 gene" under *Exclusion* Criteria. Whole-text
   keyword matching scored it as a good SOD1 match. Fixed by splitting
   eligibility text on the "Exclusion Criteria" heading and treating
   exclusion-section hits as negative.
2. NCT07294144 ("Tofersen in Non-SOD1 ALS") requires inclusion criteria
   "Prior confirmed genetic testing negative for SOD1 and FUS mutations" --
   negation *inside* the Inclusion section, invisible to a section split
   alone. Fixed with a negation-cue window check around each keyword match.

Also manually confirmed the fix by checking NCT07259980 (a tofersen SOD1
safety study) appeared in `--mutation sod1` results at all, and declared this
the tool's positive control. **This check was too shallow** -- it verified
presence, not the score/notes a user would actually see. Stage 2 found the
gap.

## Stage 2 — independent adversarial review

A separate agent, given only the code, the README, and the framing that a
real family may use this tool today, was told to find problems beyond what
Stage 1 already documented, by actually running the tool against the live
API (not just reading code). Findings, most severe first:

### 1. Network/API failure silently reported as "zero trials exist" -- BLOCKED MERGE
Pointing the API base at an unreachable host and running the tool produced:
*"No recruiting trial currently lists a site in France. Showing all 0
recruiting trials worldwide instead."* with exit code 0 -- a network blip or
ClinicalTrials.gov outage would tell a family, confidently, that zero ALS
trials exist anywhere.
**Fix:** `fetch_trials` now raises `FetchError` on any request failure
instead of logging and returning `[]`; `main()` catches it, prints an
explicit "this is a network failure, not a result" message, and returns
exit code 1. Verified: pointing `API_BASE` at an invalid host now raises
`FetchError` instead of returning an empty list.

### 2. The tool's own cited positive-control trial got its most severe false-exclusion warning -- BLOCKED MERGE
`--mutation sod1 --json` ranked NCT07259980 (the tofersen/SOD1 safety study
Stage 1 cited as proof the query works) **dead last of 182 results**, score
-5, with "WARNING: ... this trial likely EXCLUDES this genotype." The actual
exclusion text only bars concurrent enrollment in another interventional
trial -- unrelated to the mutation. Stage 1's validation checked only
presence in the result set, not this. This is Stage 1's finding #3's failure
shape (exclusion mentions the gene for an unrelated reason) hitting the
flagship validation case, which Stage 1 had marked "documented, not fixed."
**Fix:** rearchitected mutation/stage scoring. Confident +4/-5 verdicts
("likely EXCLUDES") were replaced with small +-1 nudges plus a quoted
excerpt of the actual source sentence in the note, so a human reads the real
clause instead of trusting an assertion the tool cannot reliably make.
Verified: NCT07259980 now scores +2 (from country match; mutation
contributes 0 net) and its note quotes the real sentence verbatim ("data
collected while a person with SOD1-ALS is participating in an interventional
clinical trial ... will be excluded") instead of asserting a verdict.

### 3. Country/US-state name collision silently defeats the "no results -> widen" safety net -- fixed
`--country "Georgia"` (the country) returned 8 studies, all located in
Atlanta, Georgia, USA -- ClinicalTrials.gov's `query.locn` does fuzzy
text matching, not strict country filtering. Because `studies` was
non-empty, the zero-results widen-and-explain fallback never fired, so the
user would see only these 8 US-heavy trials with no explanation of the
state/country mix-up.
**Fix:** the fallback now triggers whenever zero *fetched* results actually
have the requested country in their real location list (`trial_matches_country`),
not just when the raw fetch is empty. Verified live: `--country "Georgia"`
now triggers "No recruiting trial has a confirmed site in Georgia. Showing
all 182 recruiting trials worldwide instead."

### 4 & 5. Two more false-negative shapes beyond the one already documented -- addressed by the same rescoring fix as #2, not independently fixable
- NCT07401121 excludes genetic ALS "except for TARDBP gene variants" -- an
  explicit carve-out *for* TARDBP patients, previously scored ambiguous/low.
- NCT06891716 excludes SOD1/FUS "due to known absence of TDP-43 pathology"
  -- for a `--mutation tardbp` search this read backwards as a TARDBP
  exclusion signal.
No keyword-window heuristic reliably distinguishes these from genuine
exclusions; the same rescoring fix as #2 applies here: small score deltas,
quoted source sentence, no confident verdict. Verified live: both trials now
show mild scores (1 and 0, not -4/-5) with the actual sentence quoted in the
note.

### 6. `--mutation`/`--stage` hard-restricted to a fixed choice list, blocking real ALS genes -- fixed
`argparse choices=` rejected `--mutation tbk1` outright, even though TBK1 is
a real ALS gene cited in this repo's own `research/03_data_resources.md`
(GWAS locus). OPTN, VCP, UBQLN2, ANG, PFN1, NEK1, MATR3 had the same
problem.
**Fix:** dropped the `choices=` restriction; any string is now accepted and
matched verbatim against `MUTATION_ALIASES`/`STAGE_CUES` if not already a
known key, and the known-gene list itself was expanded to include TBK1,
OPTN, VCP, UBQLN2, ANG, PFN1, NEK1, MATR3. Verified live: `--mutation TBK1`
now runs without error.

### 7. Citations to `reviews/trial_matcher.md` pointed at a file that didn't exist -- fixed
The README and code comments cited this file three times before it was
written, meaning the only record of Stage 1's testing was the README's own
"Known limitations" section -- self-certification, not the independent
record the citations implied. This file is that fix.

### 8. Defensive null-safety on nested `.get()` chains -- fixed, unverified live
`design.get("enrollmentInfo", {}).get("count")` (and similarly for date
structs and sponsor) would raise `AttributeError` if the API ever returned
an explicit JSON `null` for a present-but-empty nested object, rather than
omitting the key. Not observed across ~230 real records checked, but not
guaranteed by the API contract either.
**Fix:** changed to `(design.get("enrollmentInfo") or {}).get("count")`
pattern throughout `extract_fields`.

### 9. `atxn2` alias missing the common prose form "ataxin-2" -- fixed
No live trial exercised this in the review, but it's a clear gap.
**Fix:** added `"ataxin-2"` / `"ataxin 2"` to the `atxn2` alias list.

## Verdict

**SURVIVES**, in the sense that the tool ships -- but only after the
independent review process did exactly what it exists to do: catch a
severe, real error (the positive control itself was badly mis-scored) that
self-testing missed. The residual limitation (findings 4/5's failure shape)
is fundamental to keyword-based matching, not a bug queue item, and is why
the tool's design was changed from "confident classifier" to "objective
scoring for structured data + verbatim quotes for text-derived signals" --
per the project's honesty rules (`README.md` operating rule 3: "a method
that can't rediscover known truth isn't trusted to find new truth"), a
scoring approach that got its own positive control backwards should not
ship making confident claims about anything else, either.
