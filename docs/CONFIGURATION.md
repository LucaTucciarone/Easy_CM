# Configuration

EasyCM uses a single YAML config file for all analysis parameters. No
contrasts file is needed -- cell types are automatically discovered from
the count matrices.

| File | Purpose |
|------|---------|
| `config/config.yaml` | Analysis parameters: inputs, filtering, thresholds, step toggles |
| `profiles/local/config.yaml` | Local execution: core count, keep-going, latency wait |
| `profiles/slurm/config.yaml` | SLURM execution: job count, sbatch template, account/partition |

---

## Config File Reference

### `inputs`

```yaml
inputs:
  seurat_object: ""                          # Path to Seurat .rds (only if make_pseudobulk = true)
  counts_dir: "data/counts"                  # Directory with pseudobulk count matrices
  sample_metadata: "data/sample_metadata.csv"
  celltype_column: "coarse_annot"            # Cell type column in metadata/Seurat
  sample_column: "samples"                   # Sample ID column
  celltypes: []                              # Empty = all discovered; or list specific ones
```

- **`seurat_object`**: Only used when `steps.make_pseudobulk: true`. The
  pipeline will subset cells per cell type and aggregate via
  `AggregateExpression()`.
- **`counts_dir`**: Used when `steps.make_pseudobulk: false`. See
  [Count Matrix Format](#count-matrix-format) below.
- **`celltypes`**: Leave empty to auto-discover all cell types from
  `counts_dir`. Specify a list to restrict analysis to specific types.

### `filtering`

```yaml
filtering:
  sample_filters:                            # Key-value metadata filters
    description_of_diabetes_status: "Control Without Diabetes"
    treatments: "no_treatment"
  min_cells_per_sample: 20                   # Minimum cells per sample (Seurat mode)
  min_gene_counts: 5                         # Minimum raw counts for gene filtering
  min_gene_proportion: 0.25                  # Gene must pass min_gene_counts in this fraction of samples
  donor_col: "rrid"                          # Column for donor deduplication
```

- **`sample_filters`**: Applied before any analysis. Each key is a metadata
  column name; the value is the required value. All filters are AND-combined.
- **`donor_col`**: When multiple samples share a donor, one is kept per
  donor (deterministic random selection with seed 123).
- **`min_gene_counts`** / **`min_gene_proportion`**: A gene must have at
  least `min_gene_counts` raw counts in at least `min_gene_proportion` of
  samples. Applied to the interest matrix only, following the OG pipeline.

### `deseq2`

```yaml
deseq2:
  design: "~ subtype + sample"
  fdr_threshold: 0.05
  lfc_threshold: 0
```

- **`design`**: The DESeq2 formula. Default `~ subtype + sample` tests
  cell-type effect while blocking on donor. In most cases this should not
  be changed.
- **`lfc_threshold`**: Minimum absolute log2 fold change. Set to 0 to
  report all DEGs (default).

### `fgsea`

```yaml
fgsea:
  fdr_threshold: 0.10
  min_gene_set_size: 10
  max_gene_set_size: 500
  gene_sets_file: "resources/gsea_files/reactome_kegg.gmt.txt"
  exclude_gene_lists:
    - "resources/gsea_files/rpl_genes.csv"
    - "resources/gsea_files/rps_genes.csv"
    - "resources/gsea_files/mtr_genes.csv"
  exclude_gene_patterns:
    - "^MT-"
    - "^NA"
  ranking_stat: "stat"
  n_top_pathways_plot: 20
```

- **`ranking_stat`**: `"stat"` uses the DESeq2 Wald statistic directly.
  `"lfc_pvalue"` uses `(-log10(pvalue)) * log2FoldChange`.
- **`exclude_gene_lists`**: CSV files containing gene symbols to remove
  before ranking. Expected format: HGNC download with "Approved symbol"
  column, or single-column gene list.
- **`exclude_gene_patterns`**: Regex patterns applied to gene names.

### `steps`

```yaml
steps:
  make_pseudobulk: false                     # true = build from Seurat; false = use counts_dir
  fgsea: true                                # Run pathway enrichment after DESeq2
```

### `outputs` and `logging`

```yaml
outputs:
  results_dir: "results"
  logs_dir: "logs"

logging:
  level: "INFO"                              # DEBUG | INFO | WARN | ERROR
```

---

## Count Matrix Format

EasyCM requires **paired** count matrices per cell type -- one for the
cell type of interest and one for all other cells. Two directory layouts
are supported:

### Layout 1: OG Pipeline (recommended for existing data)

```
data/counts/
├── cell_mtx/
│   ├── Alpha_persample_RNA_counts.tsv
│   ├── Beta_persample_RNA_counts.tsv
│   └── ...
└── allbut_mtx/
    ├── Alpha_persample_RNA_counts.tsv
    ├── Beta_persample_RNA_counts.tsv
    └── ...
```

Files are tab-separated. First column is gene name, remaining columns are
sample IDs.

### Layout 2: EasyCM Paired Naming

```
data/counts/
├── Alpha_interest.counts.csv
├── Alpha_other.counts.csv
├── Beta_interest.counts.csv
├── Beta_other.counts.csv
└── ...
```

Files are comma-separated with the same structure.

### Layout 3: Seurat Object (automated)

Set `steps.make_pseudobulk: true` and provide `inputs.seurat_object`.
The pipeline builds both matrices internally using `AggregateExpression()`.
List cell types in `inputs.celltypes`.

### Matrix structure

```
gene    SAMN001  SAMN002  SAMN003  ...
INS     142      0        88       ...
GCG     0        201      0        ...
```

Genes as rows, samples as columns. First column is the gene symbol.

---

## Sample Metadata Format

`data/sample_metadata.csv` -- one row per sample, must contain:

- The `sample_column` specified in config (e.g. `samples`)
- The `donor_col` for deduplication (e.g. `rrid`)
- All columns referenced in `sample_filters`

---

## Cell Type Discovery

When `inputs.celltypes` is empty, EasyCM auto-discovers cell types from
the count directory:

1. First checks for `{CellType}.counts.csv` files directly in `counts_dir`
2. Falls back to `cell_mtx/{CellType}_persample_RNA_counts.tsv` files
3. Cell type names are extracted from filenames

To analyze only specific cell types, list them explicitly:

```yaml
inputs:
  celltypes:
    - Beta
    - Alpha
    - Delta
```
