# Methods

## Pseudobulk design (rule 02)

For each cell type, rule 02 builds a **combined count matrix** in which every retained
sample appears twice:

- **interest**: pseudobulk counts of cells of this type in that sample
- **other**: pseudobulk counts of all remaining cells in that sample

This paired layout is what makes the design formula work. With the same donors present in
both groups, `~ subtype + sample` isolates the cell-type effect (`subtype`: interest vs other)
while absorbing all donor-level variance into the `sample` term. It is the negative-binomial
analog of a paired t-test, applied gene by gene.

### Sample filters (PankBase pancreas run)

Applied to `sample_metadata.csv` before matrix construction:

| Filter | Value |
|--------|-------|
| `derived_diabetes_status` | `Normal` |
| `treatments` | `no_treatment` |

Only samples present in both the interest and other matrices are kept.

### Donor deduplication

When multiple samples share a `donor_accession`, one is kept at random per donor with
`set.seed(123)`. This prevents pseudoreplication without introducing selection bias across runs.

### Gene filter

Applied to the **interest matrix only**, then the surviving gene set is applied to the other
matrix:

- A gene must have **≥ 5 raw counts in ≥ 25 % of samples**.

Genes that pass in the interest matrix but are zero in the other matrix are retained
(biologically meaningful: cell-type-specific expression).

### Preflight checks

Before DESeq2, rule 02 requires:

- ≥ 3 shared samples (DESeq2 dispersion minimum)
- ≥ 10 genes surviving the filter
- No all-zero columns

Failure writes a skip sentinel and propagates a `skipped_*` status downstream.

---

## DESeq2 (rule 03)

Standard pseudobulk workflow on the combined matrix:

1. **Size factors**: DESeq2 median-of-ratios.
2. **Dispersion**: gene-wise estimates shrunk to the parametric fitted curve.
3. **Wald test** on the `subtype` coefficient.

Results are extracted with `contrast = c("subtype", "interest", "other")`:

- `log2FoldChange > 0` → upregulated in the cell type of interest (positive marker)
- `log2FoldChange < 0` → higher in other cells (negative marker)

### apeglm LFC shrinkage

The raw `log2FoldChange` is noisy for low-count genes. After the Wald test, rule 03 calls
`lfcShrink(..., type = "apeglm")` on the `subtype_interest_vs_other` coefficient. The shrunk
estimate is reported as `log2FoldChange_shrunk` alongside the raw value.

Use the shrunk estimate for ranking, visualisation, and downstream comparison. The raw value
is retained for transparency and for any test that assumes unshrunk effect sizes.

### DEG counting

At FDR 0.05 and `lfc_threshold = 0` (absolute shrunk LFC), genes are counted as up- or
downregulated and reported separately in the per-cell-type summary and the volcano plot.

---

## fGSEA (rule 04)

`fgseaMultilevel(..., eps = 0.0)` — multilevel splitting for exact small p-values.

### Ranking metric

**Wald statistic** (`stat` column from DESeq2). Combines effect size and precision in one
signed value, and is the standard recommended GSEA input for DESeq2 output.

### Gene exclusions

Removed before ranking:

- Ribosomal proteins: RPL, RPS families (`rpl_genes.csv`, `rps_genes.csv`)
- Mitochondrial genes: `MT-` prefix and mt-encoded genes (`mtr_genes.csv`)
- Any regex in `fgsea.exclude_gene_patterns` (e.g. `^NA`)

These genes are highly variable across samples for reasons unrelated to cell-type biology
and would otherwise dominate pathway rankings.

### Gene sets

Default: merged Reactome + KEGG from MSigDB
(`resources/gsea_files/reactome_kegg.gmt.txt`). Only pathways with **10–500 genes** are
tested.

### Thresholds

- Pathway FDR: **0.10** (`fgsea_significant.tsv`)
- All tested pathways (any p-value) go to `fgsea_all.tsv`

---

## Outputs

### `results_deseq.tsv` (rule 03)

| Column | Meaning |
|--------|---------|
| `gene` | Gene symbol |
| `baseMean` | Mean normalised count across samples |
| `log2FoldChange` | Raw DESeq2 LFC, interest vs other |
| `log2FoldChange_shrunk` | apeglm-shrunk LFC — **use this for ranking/plots** |
| `lfcSE` | Standard error of raw LFC |
| `stat` | Wald statistic (used as fGSEA ranking metric) |
| `pvalue` | Wald p-value |
| `padj` | BH-adjusted p-value |

### `fgsea_significant.tsv` (rule 04)

Pathways with `padj ≤ 0.10`:

| Column | Meaning |
|--------|---------|
| `pathway` | Pathway name |
| `pval` / `padj` | Multilevel p-value, BH-adjusted |
| `ES` / `NES` | Enrichment score and normalised enrichment score |
| `size` | Genes in pathway after exclusions |
| `leadingEdge` | Comma-separated core enrichment genes |

---

## Aggregation (rules 05–06)

Rule 05 walks `results/{celltype}/` and produces:

- `celltype_summary.csv` — status, sample/DEG/pathway counts, error messages per cell type
- `all_markers_deseq.tsv`, `all_fgsea.tsv` — concatenated results with `celltype` column
- `log_summary.tsv` — WARN/ERROR/SKIP lines from all logs
- `status_overview.pdf` — heatmap: status colour + DEG counts

Rule 06 copies the summary to `results/pipeline_summary.csv` and renders
`pipeline_summary.pdf` (DEG bar chart, success-only, sorted descending).

---

## Skip sentinel system

Snakemake requires every declared output to exist. When a cell type is skipped or errors,
rule 02/03 writes sentinel TSVs containing `skipped=TRUE`. Downstream steps detect these
via `is_skip_sentinel()` and propagate the skip without producing misleading files.
