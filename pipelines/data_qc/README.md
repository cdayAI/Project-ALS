# data_qc - pre-flight count-matrix contamination checks

## te_contamination_check.py
Mandatory pre-flight QC for any raw count matrix before normalization.

It detects non-gene feature classes (transposable-element/repeat consensus rows,
pseudogene-style leftovers, mixed-ID rows) against a gene-reference whitelist,
quantifies their share of total raw library counts per sample, and emits a
GO / CAUTION / NO-GO normalization recommendation:

- non-gene share >= 5% -> **NO-GO** (exit 2 with `--fail-on-no-go`)
- non-gene share >= 1% -> CAUTION
- otherwise GO

Why: in GSE124439, ~1k transposable-element rows carry ~25-26% of all raw
counts. Normalizing without removing them dilutes every real gene signal by a
global factor; the first exp001 ALS signature failed silently because of this.

### Usage
```bash
# directory of GEO per-sample files (worked example, see demo_output_GSE124439/)
python pipelines/data_qc/te_contamination_check.py \
    --input-dir data/gse124439 --pattern '*_counts.txt.gz' \
    --gene-whitelist data/ref/human.gene_symbols.txt \
    --outdir pipelines/data_qc/demo_output_GSE124439 --fail-on-no-go

# single matrix
python pipelines/data_qc/te_contamination_check.py --matrix counts.tsv \
    --gene-whitelist symbols.txt --outdir /tmp/qc

# self-test
python pipelines/data_qc/te_contamination_check.py --self-test
```

### Whitelist
One gene symbol per line, or an NCBI `gene_info`-style TSV - if `Synonyms` /
`Symbol_from_nomenclature_authority` columns are present they are added too,
because older matrices use outdated symbols (SEPT7 vs SEPTIN7). Without a
whitelist the tool runs heuristic-only mode (':' pattern rule) and says so in
the report.

Reference used for the worked example: NCBI Mammalia/Homo_sapiens.gene_info
(ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz).

### Outputs
`qc_report.json` (machine-readable, includes verdict + exit code) and
`qc_report.md` (human-readable: class table, top offending features,
recommendation).
