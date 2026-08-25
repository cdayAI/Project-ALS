#!/usr/bin/env python
"""ALS trial matcher: ClinicalTrials.gov API v2 wrapper (Stream C, AGENTS.md).

Given a country/region, an optional known genetic mutation (SOD1, C9orf72,
TARDBP/TDP-43, FUS, ...), and an optional disease-stage hint, fetches
currently-recruiting ALS trials and ranks them by how well they match.

This is a triage tool, not a medical decision tool. Eligibility is only ever
confirmed by the trial site itself. See "What this tool does NOT do" in
tools/trial_matcher/README.md before trusting any ranking.

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

# Known ALS-associated genes, with common synonyms/aliases seen in eligibility
# text. Extend this as new causal genes are validated (research/01, GWAS
# Catalog EFO_0000253).
MUTATION_ALIASES = {
    "sod1": ["sod1"],
    "c9orf72": ["c9orf72", "c9-als", "c9 als", "hexanucleotide repeat"],
    "tardbp": ["tardbp", "tdp-43", "tdp43"],
    "fus": ["fus gene", "fus mutation", "fus-als"],
    "atxn2": ["atxn2"],
    "sporadic": ["sporadic als", "sals"],
}

# Heuristic-only stage cues. These are text-matched against free-form
# eligibility criteria; ClinicalTrials.gov has no structured "disease stage"
# field for ALS. Treated as a soft signal (annotate), never a hard filter,
# unless --strict is passed. See README caveats.
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


# ----------------------------------------------------------------------------
# 1. Fetch
# ----------------------------------------------------------------------------

def fetch_trials(condition="Amyotrophic Lateral Sclerosis", country=None,
                  status=("RECRUITING",), page_size=100, max_pages=10,
                  timeout=20):
    """Page through ClinicalTrials.gov API v2 studies endpoint.

    Returns the raw list of `study` objects (protocolSection-wrapped dicts)
    exactly as returned by the API -- extract_fields() flattens them.
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
        except urllib.error.URLError as exc:
            log(f"ClinicalTrials.gov request failed on page {page}: {exc}")
            break

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
    """Flatten one ClinicalTrials.gov v2 study object into a plain dict."""
    ps = study.get("protocolSection", {})
    ident = ps.get("identificationModule", {})
    status_mod = ps.get("statusModule", {})
    design = ps.get("designModule", {})
    sponsor = ps.get("sponsorCollaboratorsModule", {}).get("leadSponsor", {})
    elig = ps.get("eligibilityModule", {})
    locs = ps.get("contactsLocationsModule", {}).get("locations", []) or []

    nct_id = ident.get("nctId")
    countries = sorted({loc.get("country") for loc in locs if loc.get("country")})

    return {
        "nct_id": nct_id,
        "title": ident.get("briefTitle"),
        "status": status_mod.get("overallStatus"),
        "phase": design.get("phases") or [],
        "study_type": design.get("studyType"),
        "sponsor": sponsor.get("name"),
        "enrollment": design.get("enrollmentInfo", {}).get("count"),
        "start_date": status_mod.get("startDateStruct", {}).get("date"),
        "primary_completion_date": status_mod.get("primaryCompletionDateStruct", {}).get("date"),
        "eligibility_criteria": elig.get("eligibilityCriteria", ""),
        "min_age": elig.get("minimumAge"),
        "max_age": elig.get("maximumAge"),
        "countries": countries,
        "locations_count": len(locs),
        "url": f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else None,
    }


# ----------------------------------------------------------------------------
# 3. Score / rank
# ----------------------------------------------------------------------------

def split_eligibility(text):
    """Split free-text eligibility criteria into (inclusion, exclusion) halves
    on the first "Exclusion Criteria" heading. If no such heading is found,
    the whole text is treated as inclusion-only (exclusion="") rather than
    guessed at.

    This exists because a naive keyword search across the *whole* eligibility
    blob cannot tell "must have SOD1 mutation" from "must NOT have SOD1
    mutation" -- both contain the string "SOD1". Confirmed with a real trial
    during validation (NCT07322003 excludes SOD1/FUS/C9orf72 carriers) --
    see reviews/trial_matcher.md for the full case."""
    text = text or ""
    m = re.search(r"exclusion criteria", text, re.IGNORECASE)
    if m:
        return text[:m.start()], text[m.start():]
    return text, ""


# Negation cues checked in the ~40 chars immediately before a matched
# keyword. Exists because "Exclusion Criteria" headings alone miss inline
# negation inside the Inclusion section itself -- e.g. "Prior confirmed
# genetic testing negative for SOD1" is an INCLUSION-criteria sentence that
# still excludes SOD1 carriers. Confirmed with a real trial during
# validation (NCT07294144, "Tofersen in Non-SOD1 ALS"); see reviews/.
# This is a plain substring heuristic, not real NLP negation scope
# detection -- it will still miss negation phrased in unexpected ways.
NEGATION_CUES = [
    "without", "negative for", "non-", "non ", "absence of", "not have",
    "not carry", "not carrying", "excluding", "except", "lack of", "no ",
]


def _find_with_negation(aliases, text, window=40):
    """Find the first alias in text; report whether a negation cue appears
    in the window immediately preceding it. Returns (alias_or_None, negated)."""
    for a in aliases:
        idx = text.find(a)
        if idx == -1:
            continue
        preceding = text[max(0, idx - window):idx]
        negated = any(cue in preceding for cue in NEGATION_CUES)
        return a, negated
    return None, False


def score_trial(trial, country=None, mutation=None, stage=None):
    """Transparent, additive scoring. Returns (score, notes) -- notes explain
    every point awarded so the ranking is auditable, not a black box."""
    score = 0
    notes = []
    inc_text, exc_text = split_eligibility(trial.get("eligibility_criteria"))
    inc_text, exc_text = inc_text.lower(), exc_text.lower()
    title = (trial.get("title") or "").lower()

    if country:
        if trial.get("countries") and any(country.lower() in c.lower() for c in trial["countries"]):
            score += 3
            notes.append(f"site located in requested country ({country})")
        else:
            notes.append(f"no site found in {country} -- verify travel distance")

    if mutation:
        aliases = MUTATION_ALIASES.get(mutation.lower(), [mutation.lower()])
        inc_hit, inc_negated = _find_with_negation(aliases, inc_text + " " + title)
        exc_hit, exc_negated = _find_with_negation(aliases, exc_text)
        if exc_hit and not exc_negated:
            score -= 5
            notes.append(f"WARNING: '{exc_hit}' appears in EXCLUSION criteria -- this trial likely EXCLUDES this genotype, do not assume eligibility")
        elif exc_hit and exc_negated:
            score -= 1
            notes.append(f"ambiguous: '{exc_hit}' appears near a negation inside exclusion criteria -- verify manually, do not trust this ranking")
        elif inc_hit and inc_negated:
            score -= 5
            notes.append(f"WARNING: inclusion criteria require testing NEGATIVE for '{inc_hit}' -- this trial likely EXCLUDES confirmed carriers, do not assume eligibility")
        elif inc_hit:
            score += 4
            notes.append(f"mutation keyword '{inc_hit}' found in eligibility text/title")
        else:
            notes.append(f"no explicit mention of '{mutation}' -- may still be eligible (broad ALS trial) or may exclude this genotype; confirm with site")

    if stage:
        cues = STAGE_CUES.get(stage.lower(), [])
        exc_hit = next((c for c in cues if c in exc_text), None)
        inc_hit = next((c for c in cues if c in inc_text), None)
        if exc_hit:
            score -= 2
            notes.append(f"WARNING: stage cue '{exc_hit}' appears in EXCLUSION criteria -- may exclude this stage, confirm with site")
        elif inc_hit:
            score += 2
            notes.append(f"stage cue '{inc_hit}' matches requested stage ({stage})")
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


def rank_trials(studies, country=None, mutation=None, stage=None, strict=False):
    flat = [extract_fields(s) for s in studies]
    ranked = []
    for t in flat:
        score, notes = score_trial(t, country=country, mutation=mutation, stage=stage)
        if strict:
            # Reuse the exact same negation-aware matching as score_trial so
            # --strict can never disagree with the notes shown to the user.
            inc_text, exc_text = split_eligibility(t.get("eligibility_criteria"))
            inc_text, exc_text, title = inc_text.lower(), exc_text.lower(), (t.get("title") or "").lower()
            if mutation:
                aliases = MUTATION_ALIASES.get(mutation.lower(), [mutation.lower()])
                inc_hit, inc_negated = _find_with_negation(aliases, inc_text + " " + title)
                exc_hit, exc_negated = _find_with_negation(aliases, exc_text)
                genuinely_included = inc_hit and not inc_negated
                excluded = (exc_hit and not exc_negated) or (inc_hit and inc_negated) or (exc_hit and exc_negated)
                if excluded or not genuinely_included:
                    continue
            if country and not any(country.lower() in c.lower() for c in (t.get("countries") or [])):
                continue
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
    "keyword matching against free-text eligibility criteria. This matching "
    "can fail in BOTH directions -- confirmed during validation (see "
    "reviews/trial_matcher.md): a trial can score HIGH while actually "
    "excluding your mutation (the gene name appears in inclusion text but "
    "negated, e.g. 'tested negative for SOD1'), and a trial can score LOW "
    "or vanish under --strict while actually being a strong match (the gene "
    "name appears near 'exclusion' for an unrelated reason, e.g. excluding "
    "concurrent enrollment in another trial, not excluding the mutation "
    "itself). NEVER rule a trial in or out based on this tool alone -- read "
    "the linked eligibility criteria yourself, or better, always confirm "
    "directly with the trial site (contact info at the URL below) or "
    "through your neurologist / ALS clinic."
)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--condition", default="Amyotrophic Lateral Sclerosis")
    ap.add_argument("--country", default=None, help="e.g. 'United States', 'France'")
    ap.add_argument("--mutation", default=None, choices=sorted(MUTATION_ALIASES), help="known causal gene, if any")
    ap.add_argument("--stage", default=None, choices=sorted(STAGE_CUES), help="rough disease stage, if known")
    ap.add_argument("--status", nargs="+", default=["RECRUITING"])
    ap.add_argument("--max-pages", type=int, default=10)
    ap.add_argument("--limit", type=int, default=20, help="max rows to print")
    ap.add_argument("--strict", action="store_true", help="HIDE trials that don't match country/mutation instead of just annotating them. "
                     "DANGER: proven during validation to be able to hide a genuinely strong match "
                     "(false-negative keyword matching) -- prefer the default (non-strict) mode, which "
                     "shows everything with notes, and read the notes yourself.")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of a text table")
    args = ap.parse_args(argv)

    studies = fetch_trials(condition=args.condition, country=args.country,
                            status=args.status, max_pages=args.max_pages)
    ranked = rank_trials(studies, country=args.country, mutation=args.mutation,
                          stage=args.stage, strict=args.strict)

    fallback_note = None
    if not studies and args.country:
        # Zero recruiting trials with a site in this country is a real,
        # meaningful result (small/no ALS trial presence there) -- but it
        # must never read as "no ALS trials exist". Widen automatically and
        # say so explicitly, rather than leaving a blank table.
        log(f"0 recruiting trials with a site in {args.country} -- widening to all countries")
        studies = fetch_trials(condition=args.condition, country=None,
                                status=args.status, max_pages=args.max_pages)
        ranked = rank_trials(studies, country=args.country, mutation=args.mutation,
                              stage=args.stage, strict=False)
        fallback_note = (
            f"No recruiting trial currently lists a site in {args.country}. "
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


if __name__ == "__main__":
    main()
