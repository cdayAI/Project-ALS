# pipelines/perturbation_signatures/

Streaming LINCS L1000 Level5 reversal scoring for genotype-defined iPSC disease
signatures (the "ropinirole playbook": signature -> reversal -> annotated candidates).

Owner: c9orf72-factory (board-claimed; authorized by parent agent).
Inputs: data/lincs/ GCTX stack (READ-ONLY, owned by sprint1-repurposing agent;
never re-downloaded here). Core scorer vendored from branch exp001-sprint1-handoff.

Usage:
    from pipelines.perturbation_signatures.lincs_score import (
        read_gctx_metadata, stream_scores_multi, align_query_to_gctx)

Consumers: experiments/exp003_h007b_reversal/run.py
