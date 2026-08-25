#!/usr/bin/env python
"""ALS trial matcher: ClinicalTrials.gov API v2 wrapper (Stream C, AGENTS.md).

Given a country/region, an optional known genetic mutation (SOD1, C9orf72,
TARDBP/TDP-43, FUS, TBK1, OPTN, VCP, ...), and an optional disease-stage
hint, fetches currently-recruiting ALS trials and ranks them.

This is a triage tool, not a medical decision tool. Eligibility is only ever
confirmed by the trial site itself. See "What this tool does NOT do" in
tools/trial_matcher/README.md and reviews/trial_matcher.md before trusting
any ranking -- an independent adversarial review found real cases where
keyword matching gets the direction backwards, including on this tool's own
originally-cited positive control. That is why mutation/stage matches are
shown as excerpted source sentences, not asserted verdicts, and why there
is no --strict / hide-non-matches mode: hiding was proven capable of
hiding the single best-matching trial for a real scenario.

Data source: https://clinicaltrials.gov/api/v2/studies (NIH/NLM, public, no
API key required). Endpoint and field list are the reproducibility contract
required by AGENTS.md rule 6 -- see FIELDS/API_BASE below and README.md.

Run:
  python trial_matcher.py --country "United States" --mutation C9orf72
  python trial_matcher.py --country France --stage early --json
"""
import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

API_BASE = "https://clinicaltrials.gov/api/v2/studies"

FIELDS = [
    "NCTId", "BriefTitle", "OverallStatus", "Phase", "StudyType",
    "LeadSponsorName", "EnrollmentCount", "StartDate",
    "PrimaryCompletionDate", "EligibilityCriteria", "MinimumAge",
    "MaximumAge", "LocationCountry", "LocationCity", "LocationFacility",
]

# Known ALS-associated genes, with common synonyms/aliases seen in
# eligibility text (research/01_biology_and_therapeutics.md, GWAS Catalog
# EFO_0000253). Not a closed list: any --mutation value not found here is
# used verbatim as its own single alias, so an uncommon gene still works.
MUTATION_ALIASES = {
    "sod1": ["sod1", "sod-1"],
    "c9orf72": ["c9orf72", "c9-als", "c9 als", "hexanucleotide repeat"],
    "tardbp": ["tardbp", "tdp-43", "tdp43"],
    "fus": ["fus gene", "fus mutation", "fus-als"],
    "atxn2": ["atxn2", "ataxin-2", "ataxin 2"],
    "tbk1": ["tbk1"],
    "optn": ["optn"],
    "vcp": ["vcp"],
    "ubqln2": ["ubqln2"],
    "ang": ["angiogenin", " ang "],
    "pfn1": ["pfn1", "profilin"],
    "nek1": ["nek1"],
    "matr3": ["matr3"],
    "sporadic": ["sporadic als", "sals"],
}

# Heuristic-only stage cues -- ClinicalTrials.gov has no structured "disease
# stage" field for ALS. Same "not a closed list" behavior as mutations
# above: an unrecognized --stage value is used verbatim as its own cue.
STAGE_CUES = {
    "early": [
        "early stage", "early-stage", "recently diagnosed", "within 24 months",
        "within 12 months", "disease duration", "diagnosis within",
        "king's stage 1", "kings stage 1",
    ],
    "late": [
        "advanced", "late stage", "late-stage", "king's stage 3",
        "king's stage 4", "kings stage 3", "kings stage 4", "non-invasive ventilation",
        "gastrostomy", "peg tube",
    ],
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


class FetchError(RuntimeError):
    """Raised when the ClinicalTrials.gov API could not be reached or
    returned an error. Deliberately distinct from "the API returned zero
    matching studies" -- an independent review confirmed that conflating
    the two lets a network blip get reported to the user as "no ALS trials
    exist anywhere", the worst possible false negative for this tool.
    See reviews/trial_matcher.md finding 1."""


# ----------------------------------------------------------------------------
# 1. Fetch
# ----------------------------------------------------------------------------

def fetch_trials(condition="Amyotrophic Lateral Sclerosis", country=None,
                  status=("RECRUITING",), page_size=100, max_pages=10,
                  timeout=20):
    """Page through ClinicalTrials.gov API v2 studies endpoint.

    Returns the raw list of `study` objects (protocolSection-wrapped dicts)
    exactly as returned by the API -- extract_fields() flattens them.
    Raises FetchError on any request failure; never silently returns a
    truncated/partial result as if it were complete.
    """
    studies = []
    page_token = None
    for page in range(max_pages):
        params = {
            "query.cond": condition,
            "filter.overallStatus": "|".join(status),
            "pageSize": str(page_size),
            "fields": ",".join(FIELDS),
        }
        if country:
            params["query.locn"] = country
        if page_token:
            params["pageToken"] = page_token

        url = f"{API_BASE}?{urllib.parse.urlencode(params)}"
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            raise FetchError(
                f"ClinicalTrials.gov request failed on page {page} "
                f"(after {len(studies)} studies already fetched this call): {exc}"
            ) from exc

        batch = data.get("studies", [])
        studies.extend(batch)
        page_token = data.get("nextPageToken")
        log(f"page {page}: +{len(batch)} studies (total {len(studies)})")
        if not page_token or not batch:
            break
    return studies


# ----------------------------------------------------------------------------
# 2. Flatten
# ----------------------------------------------------------------------------

def extract_fields(study):
    """Flatten one ClinicalTrials.gov v2 study object into a plain dict.

    Every nested .get() is guarded against an explicit JSON null (not just
    a missing key) with `... or {}` -- a live check across ~230 real
    records never hit this, but the API contract does not promise it can't
    happen (reviews/trial_matcher.md finding 8)."""
    ps = study.get("protocolSection") or {}
    ident = ps.get("identificationModule") or {}
    status_mod = ps.get("statusModule") or {}
    design = ps.get("designModule") or {}
    sponsor = (ps.get("sponsorCollaboratorsModule") or {}).get("leadSponsor") or {}
    elig = ps.get("eligibilityModule") or {}
    locs = (ps.get("contactsLocationsModule") or {}).get("locations") or []

    nct_id = ident.get("nctId")
    countries = sorted({loc.get("country") for loc in locs if loc.get("country")})

    return {
        "nct_id": nct_id,
        "title": ident.get("briefTitle"),
        "status": status_mod.get("overallStatus"),
        "phase": design.get("phases") or [],
        "study_type": design.get("studyType"),
        "sponsor": sponsor.get("name"),
        "enrollment": (design.get("enrollmentInfo") or {}).get("count"),
        "start_date": (status_mod.get("startDateStruct") or {}).get("date"),
        "primary_completion_date": (status_mod.get("primaryCompletionDateStruct") or {}).get("date"),
        "eligibility_criteria": elig.get("eligibilityCriteria") or "",
        "min_age": elig.get("minimumAge"),
        "max_age": elig.get("maximumAge"),
        "countries": countries,
        "locations_count": len(locs),
        "url": f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else None,
    }


def trial_matches_country(trial, country):
    return bool(trial.get("countries")) and any(
        country.lower() in c.lower() for c in trial["countries"]
    )


# ----------------------------------------------------------------------------
# 3. Score / rank
# ----------------------------------------------------------------------------

def split_eligibility(text):
    """Split free-text eligibility criteria into (inclusion, exclusion)
    halves on the first "Exclusion Criteria" heading. If no such heading is
    found, the whole text is treated as inclusion-only."""
    text = text or ""
    m = re.search(r"exclusion criteria", text, re.IGNORECASE)
    if m:
        return text[:m.start()], text[m.start():]
    return text, ""


NEGATION_CUES = [
    "without", "negative for", "non-", "non ", "absence of", "not have",
    "not carry", "not carrying", "excluding", "except", "lack of", "no ",
]


def _find_sentence(text, alias):
    """Return the bullet/sentence containing `alias`, for the user to read
    verbatim -- see module docstring on why this replaced a confident
    include/exclude verdict."""
    idx = text.find(alias)
    if idx == -1:
        return None
    start = max(text.rfind("\n", 0, idx), text.rfind(". ", 0, idx))
    start = start + 1 if start != -1 else 0
    end_candidates = [e for e in (text.find("\n", idx), text.find(". ", idx)) if e != -1]
    end = min(end_candidates) if end_candidates else len(text)
    return text[start:end].strip(" \t*-\n")


def _find_with_negation(aliases, text, window=40):
    """Find the first alias in text; report whether a negation cue appears
    in the window immediately preceding it, and the surrounding sentence.
    Returns (alias_or_None, negated, sentence_or_None). This is a plain
    substring/window heuristic, not real NLP negation-scope detection --
    an independent adversarial review confirmed it still gets the
    direction wrong on real trials in both directions (reviews/
    trial_matcher.md findings 2, 4, 5). That is why callers must treat the
    boolean as a weak hint and always surface `sentence` for the human to
    read, rather than asserting it as a verdict."""
    for a in aliases:
        idx = text.find(a)
        if idx == -1:
            continue
        preceding = text[max(0, idx - window):idx]
        negated = any(cue in preceding for cue in NEGATION_CUES)
        return a, negated, _find_sentence(text, a)
    return None, False, None


def score_trial(trial, country=None, mutation=None, stage=None):
    """Additive scoring with small weights for text-derived (mutation/
    stage) signals and larger weights for structured-data signals
    (country, phase, study type) that don't rely on free-text
    interpretation. Returns (score, notes).

    Mutation/stage weights are deliberately small (+-1) rather than the
    confident +4/-5 an earlier version used: an adversarial review found
    that version's own positive-control trial (a tofersen safety study
    explicitly FOR SOD1-ALS patients) scored the single most negative,
    "likely EXCLUDES" mutation score of all 182 fetched trials, because
    its exclusion criteria mention "SOD1-ALS" while excluding *concurrent
    enrollment in another trial*, not the mutation. See reviews/
    trial_matcher.md finding 2. Small weights mean a wrong call there
    nudges rank instead of burying a genuinely strong match."""
    score = 0
    notes = []
    inc_text, exc_text = split_eligibility(trial.get("eligibility_criteria"))
    inc_text, exc_text = inc_text.lower(), exc_text.lower()
    title = (trial.get("title") or "").lower()

    if country:
        if trial_matches_country(trial, country):
            score += 3
            notes.append(f"site located in requested country ({country})")
        else:
            notes.append(f"no site found in {country} -- verify travel distance")

    if mutation:
        aliases = MUTATION_ALIASES.get(mutation.lower(), [mutation.lower()])
        inc_hit, inc_negated, inc_sentence = _find_with_negation(aliases, inc_text + " " + title)
        exc_hit, exc_negated, exc_sentence = _find_with_negation(aliases, exc_text)
        if exc_hit:
            score += -1 if not exc_negated else 0
            sentence = exc_sentence or ""
            notes.append(
                f"'{exc_hit}' mentioned in EXCLUSION criteria -- READ THIS YOURSELF, "
                f"do not trust a +/- score for it: \"{sentence[:220]}\""
            )
        elif inc_hit:
            score += -1 if inc_negated else 1
            sentence = inc_sentence or ""
            notes.append(
                f"'{inc_hit}' mentioned in inclusion criteria -- READ THIS YOURSELF, "
                f"do not trust a +/- score for it: \"{sentence[:220]}\""
            )
        else:
            notes.append(f"no explicit mention of '{mutation}' -- may still be eligible (broad ALS trial) or may exclude this genotype; confirm with site")

    if stage:
        cues = STAGE_CUES.get(stage.lower(), [stage.lower()])
        exc_hit, _, exc_sentence = _find_with_negation(cues, exc_text, window=0)
        inc_hit, _, inc_sentence = _find_with_negation(cues, inc_text, window=0)
        if exc_hit:
            score -= 1
            notes.append(f"stage cue '{exc_hit}' appears in exclusion criteria: \"{(exc_sentence or '')[:220]}\" -- confirm with site")
        elif inc_hit:
            score += 1
            notes.append(f"stage cue '{inc_hit}' matches requested stage ({stage}): \"{(inc_sentence or '')[:220]}\"")
        else:
            notes.append(f"no explicit {stage}-stage cue found -- stage matching is heuristic, not authoritative")

    if trial.get("study_type") == "INTERVENTIONAL":
        score += 1
        notes.append("interventional trial (tests a treatment, not purely observational)")

    phases = trial.get("phase") or []
    if any(p in ("PHASE2", "PHASE3") for p in phases):
        score += 1
        notes.append(f"phase {'/'.join(phases)}")

    return score, notes


def rank_trials(studies, country=None, mutation=None, stage=None):
    flat = [extract_fields(s) for s in studies]
    ranked = []
    for t in flat:
        score, notes = score_trial(t, country=country, mutation=mutation, stage=stage)
        t["match_score"] = score
        t["match_notes"] = notes
        ranked.append(t)
    ranked.sort(key=lambda t: t["match_score"], reverse=True)
    return ranked


# ----------------------------------------------------------------------------
# 4. Output
# ----------------------------------------------------------------------------

def format_table(ranked, limit=20):
    """Text report. Country lists are never truncated -- a hidden country was
    the source of a real confusion bug during validation (see reviews/)."""
    lines = []
    for t in ranked[:limit]:
        phase = "/".join(t.get("phase") or []) or "N/A"
        countries = ", ".join(t.get("countries") or []) or "location TBD"
        lines.append(f"[{t['match_score']:>2}] {t['nct_id']}  {phase:<8}  {t['study_type'] or '?':<13}  {t['title']}")
        lines.append(f"       sites: {countries}")
        for note in t["match_notes"]:
            lines.append(f"       - {note}")
        lines.append(f"       {t['url']}")
    return "\n".join(lines)


DISCLAIMER = (
    "This is a triage aid, NOT medical advice and NOT an eligibility "
    "determination. Rankings are based on public trial metadata and simple "
    "keyword matching against free-text eligibility criteria, which cannot "
    "reliably tell inclusion from exclusion or parse negation -- confirmed "
    "on real trials during independent review (see reviews/trial_matcher.md). "
    "Mutation and stage notes quote the actual source sentence for a reason: "
    "read it yourself rather than trusting the +/- score. NEVER rule a "
    "trial in or out based on this tool alone -- confirm directly with the "
    "trial site (contact info at the URL below) or through your "
    "neurologist / ALS clinic."
)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--condition", default="Amyotrophic Lateral Sclerosis")
    ap.add_argument("--country", default=None, help="e.g. 'United States', 'France'")
    ap.add_argument("--mutation", default=None, help=f"known causal gene, if any. Common ones: {', '.join(sorted(MUTATION_ALIASES))} -- but any gene name works, matched verbatim if not in this list")
    ap.add_argument("--stage", default=None, help=f"rough disease stage, if known. Common ones: {', '.join(sorted(STAGE_CUES))} -- any free-text cue works, matched verbatim if not in this list")
    ap.add_argument("--status", nargs="+", default=["RECRUITING"])
    ap.add_argument("--max-pages", type=int, default=10)
    ap.add_argument("--limit", type=int, default=20, help="max rows to print")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of a text table")
    args = ap.parse_args(argv)

    try:
        studies = fetch_trials(condition=args.condition, country=args.country,
                                status=args.status, max_pages=args.max_pages)
    except FetchError as exc:
        print(f"ERROR: could not reach ClinicalTrials.gov: {exc}", file=sys.stderr)
        print("This is a NETWORK/API failure, not a result -- it does NOT mean zero ALS trials exist. Try again, or check https://clinicaltrials.gov directly.", file=sys.stderr)
        return 1

    ranked = rank_trials(studies, country=args.country, mutation=args.mutation, stage=args.stage)

    fallback_note = None
    country_hits = sum(1 for t in ranked if args.country and trial_matches_country(t, args.country))
    if args.country and country_hits == 0:
        # Either zero studies were fetched, or `query.locn` fuzzy-matched on
        # something else entirely (confirmed live: --country "Georgia" the
        # country returns only Georgia-USA-state trials -- reviews/
        # trial_matcher.md finding 3). Either way, zero genuine matches to
        # the requested country must never be presented as "these are your
        # options" without explanation, and must never look like "no ALS
        # trials exist".
        log(f"0 trials with a confirmed site in {args.country} ({len(studies)} fetched) -- widening to all countries")
        try:
            studies = fetch_trials(condition=args.condition, country=None,
                                    status=args.status, max_pages=args.max_pages)
        except FetchError as exc:
            print(f"ERROR: could not reach ClinicalTrials.gov while widening search: {exc}", file=sys.stderr)
            return 1
        ranked = rank_trials(studies, country=args.country, mutation=args.mutation, stage=args.stage)
        fallback_note = (
            f"No recruiting trial has a confirmed site in {args.country}. "
            f"Showing all {len(ranked)} recruiting trials worldwide instead -- "
            f"some may accept remote/international participants or have sites "
            f"reachable by travel. Check each trial's contact info directly."
        )

    if args.json:
        print(json.dumps({"disclaimer": DISCLAIMER, "fallback_note": fallback_note, "results": ranked}, indent=2))
    else:
        print(DISCLAIMER + "\n", file=sys.stderr)
        if fallback_note:
            print(fallback_note + "\n")
        print(f"{len(ranked)}/{len(studies)} recruiting trials matched (of those fetched); showing top {min(args.limit, len(ranked))}:\n")
        print(format_table(ranked, limit=args.limit))
    return 0


if __name__ == "__main__":
    sys.exit(main())
