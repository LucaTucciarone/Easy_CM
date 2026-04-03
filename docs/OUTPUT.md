# Output Files

Results can be written to any directory by setting `outputs.results_dir` in
your config file. The structure within that directory is:

```
results/
├── pipeline_summary.pdf                    step 06 -- DEG bar chart across cell types
├── pipeline_summary.csv                    step 06 -- one row per cell type
└── {celltype}/
    ├── intermediates/                           working files
    │   ├── coldata.csv                          step 02 -- combined sample metadata
    │   ├── counts_filtered.csv                  step 02 -- combined count matrix (interest + other)
    │   ├── counts_interest.csv                  step 02 -- interest-only counts
    │   ├── counts_other.csv                     step 02 -- other-only counts
    │   └── model_info.csv                       step 03 -- formula, gene/sample counts, DEG stats
    │
    ├── finals/                                  clean result tables
    │   ├── results_deseq.tsv                    step 03 -- DESeq2 results (all genes)
    │   ├── fgsea_all.tsv                        step 04 -- all tested pathways
    │   └── fgsea_significant.tsv                step 04 -- significant pathways only
    │
    ├── plots/
    │   ├── volcano.pdf                          step 03 -- volcano plot
    │   └── fgsea_pathways.pdf                   step 04 -- top pathway bar chart
    │
    └── fgsea.done                               step 04 -- Snakemake touch file
│
└── summary/                                     cross-celltype aggregation (step 05)
    ├── celltype_summary.csv                     one row per cell type, all counts + status
    ├── all_markers_deseq.tsv                    merged DESeq2 results (all cell types)
    ├── all_fgsea.tsv                            merged fGSEA results (all cell types)
    ├── status_overview.pdf                      status heatmap (cell type x status + DEG counts)
    └── log_summary.tsv                          all flagged log lines (WARN/ERROR/SKIP)
```

---

## DESeq2 Result Columns

`results_deseq.tsv` contains one row per gene:

| Column | Description |
|--------|-------------|
| `gene` | Gene symbol |
| `baseMean` | Mean normalized count across all samples (interest + other) |
| `log2FoldChange` | DESeq2 MLE estimate. Positive = upregulated in cell type of interest |
| `lfcSE` | Standard error of log2FoldChange |
| `stat` | Wald statistic (default fGSEA ranking metric) |
| `pvalue` | Raw p-value |
| `padj` | BH-adjusted p-value |

### Merged file

`all_markers_deseq.tsv` contains the same columns as above, plus:

| Column | Description |
|--------|-------------|
| `celltype` | Which cell type produced this row |

This file is the recommended input for cross-celltype analyses -- load it
and filter by `celltype` and `padj` as needed.

---

## fGSEA Result Columns

`fgsea_all.tsv` and `fgsea_significant.tsv` share this schema:

| Column | Description |
|--------|-------------|
| `pathway` | Pathway name (from GMT file) |
| `pval` | Raw p-value |
| `padj` | BH-adjusted p-value |
| `log2err` | Log2 fold enrichment error |
| `ES` | Enrichment score |
| `NES` | Normalized enrichment score. Positive = enriched in cell type markers |
| `size` | Number of pathway genes found in the ranked list |
| `leadingEdge` | Space-separated list of leading-edge genes |

---

## Model Info Columns

`model_info.csv` records metadata about the DESeq2 run:

| Column | Description |
|--------|-------------|
| `celltype` | Cell type name |
| `formula` | Design formula used (e.g. `~ subtype + sample`) |
| `n_genes` | Number of genes in the count matrix |
| `n_columns` | Number of columns (2 x n_samples) |
| `n_degs` | Total DEGs at configured FDR threshold |
| `n_up` | Upregulated markers |
| `n_down` | Downregulated markers |
| `fdr_threshold` | FDR threshold used |

---

## Coldata Columns

`coldata.csv` describes each column of the combined count matrix:

| Column | Description |
|--------|-------------|
| `test` | Column name in the count matrix (e.g. `interest_SAMN001`) |
| `subtype` | Factor: `"interest"` or `"other"` |
| `sample` | Donor/sample ID (shared between interest and other) |

---

## Pipeline Status Values

The `status` column in `celltype_summary.csv` uses a fine-grained taxonomy:

| Status | Meaning | When |
|--------|---------|------|
| `success` | DESeq2 completed, markers found | DEGs > 0 at FDR threshold |
| `success_no_significant` | DESeq2 ran, no DEGs | All genes have padj >= threshold |
| `skipped_no_samples` | Zero samples | After metadata filtering or cell type not found in data |
| `skipped_min_cells` | Too few samples | Fewer than 3 samples after deduplication |
| `skipped_preflight` | Gene filter failed | Fewer than 10 genes passed filtering |
| `skipped` | Other skip reason | Unclassified -- check `error_message` for details |
| `error` | Analysis crashed | DESeq2 or upstream step threw an error |
| `not_run` | No data found | No log entries or output files for this cell type |

The `error_message` column is populated for non-success statuses,
extracted from structured log entries.

---

## Log Files

Per-celltype logs are written to `logs/{celltype}.log` with structured
entries:

```
[2026-04-03 19:43:38] [INFO ] celltype=Beta | step=start | Preparing pseudobulk for celltype=Beta
[2026-04-03 19:43:38] [INFO ] celltype=Beta | step=load_metadata | n_samples=191, n_cols=86
[2026-04-03 19:43:38] [SKIP ] celltype=X | step=no_samples | reason=no_samples | ...
```

Error-only logs are written to `logs/{celltype}.errors.log` for quick
triage. Step 05 collects all flagged lines into
`results/summary/log_summary.tsv`.
