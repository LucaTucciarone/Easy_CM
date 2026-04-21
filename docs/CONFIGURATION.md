# Configuration

One YAML config drives the analysis. Cell types are auto-discovered from the count
directory unless `inputs.celltypes` is set.

| File | Purpose |
|------|---------|
| `config/config.yaml` | Analysis parameters |
| `profiles/local/config.yaml` | Local execution: cores, keep-going, latency wait |
| `profiles/slurm/config.yaml` | SLURM execution (not tested on ophelia) |

---

## Config reference

### `inputs`

```yaml
inputs:
  counts_dir: "data/counts"
  sample_metadata: "data/sample_metadata.csv"
  celltype_column: "coarse_annot"
  sample_column: "samples"
  celltypes: []              # empty = auto-discover
```

| Key | Meaning |
|-----|---------|
| `counts_dir` | Directory holding the pseudobulk matrices from EasyPseudobulk |
| `sample_metadata` | CSV with one row per sample |
| `celltype_column` | Metadata column containing cell-type labels (used for validation) |
| `sample_column` | Column carrying sample IDs (must match count-matrix column headers) |
| `celltypes` | Restrict to a list, or leave empty to run every discovered cell type |

### `filtering`

```yaml
filtering:
  sample_filters:
    derived_diabetes_status: "Normal"
    treatments: "no_treatment"
  min_gene_counts: 5
  min_gene_proportion: 0.25
  donor_col: "donor_accession"
```

| Key | Meaning |
|-----|---------|
| `sample_filters` | Key/value metadata filters, AND-combined. Applied before matrix construction |
| `min_gene_counts` | Minimum raw counts per sample for gene to pass (interest matrix) |
| `min_gene_proportion` | Gene must hit `min_gene_counts` in this fraction of samples |
| `donor_col` | Column used to keep one sample per donor (seed 123) |

### `deseq2`

```yaml
deseq2:
  design: "~ subtype + sample"
  fdr_threshold: 0.05
  lfc_threshold: 0
```

| Key | Meaning |
|-----|---------|
| `design` | DESeq2 formula. Default pairs on donor — do not change unless you know why |
| `fdr_threshold` | `padj` cutoff for DEG counts |
| `lfc_threshold` | Absolute shrunk-LFC cutoff for DEG counts (0 = all significant genes) |

Rule 03 always runs `lfcShrink(..., type = "apeglm")` and writes `log2FoldChange_shrunk`.

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

| Key | Meaning |
|-----|---------|
| `ranking_stat` | `stat` = Wald statistic (recommended). `lfc_pvalue` = `-log10(p) * LFC` |
| `exclude_gene_lists` | CSVs of gene symbols to drop before ranking |
| `exclude_gene_patterns` | Regex filters applied to gene names |
| Pathway size | Only pathways with 10–500 genes tested |

### `steps`

```yaml
steps:
  fgsea: true
```

Toggle fGSEA off to stop after DESeq2.

### `outputs` and `logging`

```yaml
outputs:
  results_dir: "results"
  logs_dir: "logs"

logging:
  level: "INFO"          # DEBUG | INFO | WARN | ERROR
```

---

## Count matrix format

EasyCM consumes the paired TSVs produced by the **EasyPseudobulk** notebook:

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

Tab-separated. First column = gene symbol, remaining columns = sample IDs matching
`inputs.sample_column` values in the metadata CSV.

```
gene    SAMN001  SAMN002  SAMN003
INS     142      0        88
GCG     0        201      0
```

### Cell-type discovery

When `inputs.celltypes` is empty, rule 01 scans `counts_dir/cell_mtx/` and extracts cell
type names from `{CellType}_persample_RNA_counts.tsv`. A matching file must exist in
`counts_dir/allbut_mtx/` or the cell type is rejected.

---

## Sample metadata format

`sample_metadata.csv`, one row per sample. Must contain:

- `inputs.sample_column` (sample IDs matching matrix columns)
- `filtering.donor_col` (donor key for deduplication)
- Every column referenced in `filtering.sample_filters`

---

## Legacy note

Earlier versions of EasyCM supported building pseudobulk from a Seurat `.rds` via
`steps.make_pseudobulk: true` and `inputs.seurat_object`. Both keys are removed.
Pseudobulk construction now lives entirely in the separate **EasyPseudobulk** notebook;
EasyCM only consumes its output.
