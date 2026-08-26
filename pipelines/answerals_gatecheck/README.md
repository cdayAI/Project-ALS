# pipelines/answerals_gatecheck/

Go/no-go gate checker for hypotheses H-011 and H-012 (AnswerALS iPS-neuron RNA-seq).

Runs identically on MOCK metadata (shaped per PMID 35115730 / PMC8825283 supplementary
descriptions) and on REAL portal metadata, so gate-checking executes the moment real
counts arrive. The mock NEVER feeds hypothesis testing - it exists only to validate the
code path (per _LESSONS.md L4/L8: gates run before any module test).

Schema (sample_metadata.csv - one row per RNA-seq sample):
    donor_id            unique donor
    sample_id           RNA-seq sample/library id
    genotype_group      one of: C9_carrier | nonC9_ALS | control   (control = unaffected donor)
    sex                 M / F / unknown
    c9_repeat_size      modal repeat count; >26 for carriers per paper definition; NA otherwise
    culture_batch       differentiation batch label
    rna_qc_pass         TRUE/FALSE library QC flag from portal

Usage:
    python make_mock_metadata.py [--n-per-group N]      # writes outputs/mock_sample_metadata.csv
    python gate_check.py --metadata <csv>               # applies H-011/H-012 pre-registered gates
