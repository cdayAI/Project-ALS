# Project ALS

**An open, AI-accelerated research project targeting amyotrophic lateral sclerosis (ALS).**

Human judgment + AI agents working in public. Everything here — hypotheses, code,
results, and failed ideas — is open by default, because patients don't have time
for closed pipelines and duplicated effort.

## Why this exists

ALS has no cure. Four approved drugs extend life by roughly 2-6 months combined.
Meanwhile AI can now compress the discovery loop: mining literature, screening
public multi-omics data for repurposable drugs, designing candidate molecules,
and killing bad hypotheses *before* humans waste years on them.

This project runs that loop in the open:

> HYPOTHESIZE -> EXPERIMENT (in silico) -> ADVERSARIAL REVIEW -> TRIAGE -> repeat

Human attention goes only where judgment is irreplaceable. Wet-lab validation
happens through partners, purchased assays, and collaboration.

## Repository layout

| Path | Purpose |
|---|---|
| `research/` | Literature-grounded research briefs with citations |
| `hypotheses/` | The hypothesis ledger - every idea is structured & falsifiable |
| `experiments/` | One folder per run: config, code, outputs, verdict |
| `pipelines/` | Reusable analysis pipelines (differential expression, signature reversal, ...) |
| `reviews/` | Adversarial agent reviews of every result |
| `digests/` | Periodic human-readable summaries of what survived review |

## Operating rules

1. No result counts until an adversarial review exists.
2. Every hypothesis declares falsification criteria BEFORE testing.
3. Positive controls are mandatory: a method that can't rediscover known truth isn't trusted to find new truth.
4. All state lives in this repo, not in anyone's chat history.

## Current status

- [x] Research factory design (`docs_plan.md`)
- [x] Disease & therapeutics landscape brief (`research/01_biology_and_therapeutics.md`)
- [x] Public data resource survey with verified access (`research/03_data_resources.md`)
- [x] AI methods brief with honest capability ledger (`research/02_ai_methods.md`)
- [x] Cross-brief synthesis & stream prioritization (`research/00_synthesis.md`)
- [ ] Pilot Sprint #1: transcriptome-based drug repurposing loop (exp001 running)
- [x] Stream C: patient-facing trial matcher (`tools/trial_matcher/`, `reviews/trial_matcher.md`)

## Who this is for

- **Patients & families**: digests will summarize what's real vs hype; trial-matching tooling is planned.
- **Researchers**: a transparent, reproducible pipeline you can audit, reuse, or challenge.
- **AI builders**: a real domain where careful falsification beats demo-ware.

## Contributing

We especially need:
- Domain scientists willing to attack our conclusions (that's the point)
- Access to wet-lab validation (iPSC motor-neuron assays)
- Data-wrangling and bioinformatics help

Open an issue or a PR. Adversarial input is the most valuable contribution of all.

## Disclaimer

This repository is research tooling, not medical advice. Clinical decisions belong
to patients and their physicians. Nothing here has been validated in the lab or clinic.

## License

MIT
