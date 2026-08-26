# AnswerALS portal access + metadata pull log (H-011/H-012 gate locking)

date_probed: 2026-08-26 (UTC)
agent: c9orf72-factory
authorization: parent agentmsg_e86a3f85 (free-tier registration authorized)
account_email_used: NONE - REGISTRATION BLOCKED (see below)
gate_status: UNRESOLVED AT METADATA STAGE - documented blocker, stopped per instruction
("if gates fail at metadata stage, document and stop - that IS a result")

## What was attempted

1. Portal endpoints probed:
   - https://dataportal.answerals.org/home -> 200 (SPA shell, "Neuromine Data Portal")
   - /env.js reveals backend APIs (Azure app services):
     search-api, ext-int-api, user-management-api, billing-api, doi-management-api,
     index-management-api, statistics-api, wordpress (all *.azurewebsites.net/api)
2. API probing WITHOUT credentials: all tested routes (/search, /datasets, /metadata,
   /files, user-management /register|/users|/roles) return 404 - routes are versioned/
   auth-gated. No swagger/openapi exposed. No blob-storage URLs embedded in SPA bundle.
   No public GitHub org repos; no sitemap. Conclusion: metadata pull REQUIRES a
   registered, email-verified portal account.
3. Programmatic registration is not possible for this agent: registration requires an
   email address we control for verification. None exists in this environment.

## Verified publication-level facts usable meanwhile (PMID 35115730, PMC8825283 via E-utilities)

- >1,000 ALS/control participants enrolled; ~850+ iPS cell lines banked (one clone/donor)
- 830 participants WGS'd; 41 ALS patients + 4 pre-fALS subjects carry C9orf72 expansions
  >26 repeats (Fig 4f / Suppl Table 15)
- 217 control/ALS iPS lines had spinal-neuron cultures evaluated for cell-type markers
- RNA-seq/proteomics/ATAC generated from iPS neurons at 32-day differentiation (diMN protocol)
- Published C9-vs-control RNA-seq splicing analysis exists ("male C9 samples vs male controls",
  SE 52% / RI 35% enrichment) - confirms per-donor RNA-seq with genotype labels EXISTS on the
  portal, but exact per-genotype donor counts are NOT stated in accessible text (supplementary
  tables only)

## Gate implications for H-011/H-012

Known upper bound: >=41 C9 carriers exist in the population, so the n>=15 C9 gate is
PLAUSIBLE but UNVERIFIED for the iPS-neuron RNA-seq subset specifically. Controls count
unknown. Gates remain UNLOCKED until first authenticated metadata pull.

## Unblock path (needs human)

1. Human registers at https://dataportal.answerals.org (free tier) with their email;
   record the email here.
2. Hand agent the session/API token or download the metadata tables (curated expression +
   sample metadata are GB-scale) into data/gse_answerals/ (gitignored).
3. Agent then locks gates per H-011/H-012 pre-registration BEFORE any module testing.
