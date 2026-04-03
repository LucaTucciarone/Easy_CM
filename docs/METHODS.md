# Methods

## Pseudobulk Design (step 02)

EasyCM discovers cell-type marker genes by comparing each cell type against
all other cells in the dataset. For each cell type, the pipeline constructs
a **combined count matrix** where every sample appears twice:

- **interest**: pseudobulk counts from cells of this type
- **other**: pseudobulk counts from all remaining cells

This paired design is critical. By including the same samples in both
groups, the DESeq2 `~ subtype + sample` formula can isolate the cell-type
effect (`subtype`) while perfectly controlling for donor-to-donor variation
(`sample`). This is equivalent to a paired t-test generalized to the
negative binomial framework.

### Sample filtering

Before building the combined matrix:

1. **Metadata filters** are applied (e.g. `diabetes_status: "Control"`),
   restricting analysis to a homogeneous cohort.
2. **Sample matching**: only samples present in both the interest and other
   matrices are retained.
3. **Donor deduplication**: when multiple samples share a donor (via
   `donor_col`), one is randomly kept per donor (seed 123 for
   reproducibility). This prevents pseudoreplication.

### Gene filtering

Genes are filtered on the **interest matrix only** (following the OG
pipeline convention):

- A gene must have at least `min_gene_counts` raw counts in at least
  `min_gene_proportion` of samples.
- Default: 5 counts in 25% of samples.
- After filtering, the same gene set is applied to both matrices.

### Preflight checks

Before proceeding to DESeq2, the pipeline verifies:

- At least 3 shared samples exist (minimum for DESeq2 dispersion estimation)
- At least 10 genes pass filtering
- Zero-count columns are removed from the combined matrix

Cell types that fail these checks are **skipped** with a specific status
(`skipped_no_samples`, `skipped_min_cells`, `skipped_preflight`) and
sentinel files are written to satisfy Snakemake's output contract.

---

## DESeq2 (step 03)

Standard DESeq2 pseudobulk analysis with the design formula
`~ subtype + sample`:

1. **Size factor estimation** via DESeq2's median-of-ratios method
2. **Dispersion estimation** -- gene-wise, then shrunk to the fitted curve
3. **Wald test** for the `subtype` coefficient

Results are extracted with `contrast = c("subtype", "interest", "other")`:

- **Positive log2FoldChange** = gene is upregulated in the cell type of interest
  (i.e. a positive marker)
- **Negative log2FoldChange** = gene is downregulated in the cell type
  (i.e. more highly expressed in other cell types)

### DEG counting

DEGs are counted at the configured thresholds:

- `deseq2.fdr_threshold` (default 0.05) for adjusted p-value
- `deseq2.lfc_threshold` (default 0) for minimum absolute fold change

Both upregulated and downregulated genes are counted and reported separately.

### Volcano plot

A volcano plot is generated per cell type showing:
- Upregulated markers in red
- Downregulated genes in blue
- Horizontal dashed line at the FDR threshold

---

## fGSEA (step 04)

Fast Gene Set Enrichment Analysis using `fgseaMultilevel()` with `eps = 0.0`
for exact p-value computation.

### Gene ranking

Before ranking, genes are filtered to remove:

1. **Ribosomal genes**: RPL and RPS family (from `rpl_genes.csv` and
   `rps_genes.csv`)
2. **Mitochondrial genes**: MT- prefix and mitochondrially-encoded genes
   (from `mtr_genes.csv`)
3. **Additional patterns**: regex patterns from `exclude_gene_patterns`
   (e.g. `^NA` for non-annotated genes)

These highly variable housekeeping genes would otherwise dominate pathway
results without providing biological insight.

**Ranking metric** (configurable via `fgsea.ranking_stat`):

| Option | Formula | When to use |
|--------|---------|-------------|
| `stat` (default) | DESeq2 Wald statistic | Recommended -- combines effect size and precision |
| `lfc_pvalue` | `(-log10(pvalue)) * log2FoldChange` | OG pipeline method -- directional, p-value weighted |

### Gene sets

Default: merged Reactome + KEGG pathways from MSigDB
(`resources/gsea_files/reactome_kegg.gmt.txt`). To use a different
collection, set `fgsea.gene_sets_file` in your config.

### Pathway size filtering

Only pathways with gene counts between `min_gene_set_size` (default 10)
and `max_gene_set_size` (default 500) are tested. This removes pathways
that are too small (unreliable enrichment scores) or too broad (low
specificity).

### Output

- **`fgsea_all.tsv`**: all tested pathways with NES, p-value, and adjusted
  p-value
- **`fgsea_significant.tsv`**: pathways passing `fgsea.fdr_threshold`
  (default 0.10)
- **`fgsea_pathways.pdf`**: bar chart of top N pathways by NES, colored by
  enrichment direction

---

## Aggregation (step 05)

Step 05 walks all cell type result directories and builds:

1. **`celltype_summary.csv`**: one row per cell type with status,
   sample counts, DEG counts, pathway counts, and error messages.
   Status is determined by checking output files and log entries.

2. **`all_markers_deseq.tsv`**: concatenated DESeq2 results across all
   cell types (with added `celltype` column).

3. **`all_fgsea.tsv`**: concatenated fGSEA results across all cell types.

4. **`log_summary.tsv`**: all WARN, ERROR, and SKIP log lines extracted
   from per-celltype logs.

5. **`status_overview.pdf`**: heatmap showing status (color-coded) and
   DEG counts (numeric labels) for every cell type.

---

## Pipeline Summary (step 06)

Generates the final outputs:

- **`pipeline_summary.csv`**: copy of `celltype_summary.csv` at the results
  root for convenience.
- **`pipeline_summary.pdf`**: bar chart of DEG counts per cell type,
  ordered from most to fewest markers. Only cell types with `success`
  status are included.

---

## Skip Sentinel System

Snakemake requires all declared output files to exist. When a cell type is
skipped (e.g. too few samples) or errors out, the pipeline writes
placeholder "sentinel" files containing `skipped=TRUE` for each expected
output. This satisfies Snakemake's file contract without generating
misleading result files. Downstream steps detect sentinels via
`is_skip_sentinel()` and propagate the skip cleanly.
