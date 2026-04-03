# EasyCM

A modular R pipeline for cell-type marker gene discovery via pseudobulk
DESeq2 and fGSEA pathway enrichment. Designed for scRNA-seq data
(PanKbase and similar resources) but works with any pseudobulk count matrices.

Orchestrated by Snakemake for local or SLURM cluster execution.

---

## Pipeline Overview

Six steps. Steps 02--04 run in parallel across every **cell type**:

```
01  Validate config + inputs             (once per run -- soft gate)
02  Prepare pseudobulk matrices          (per celltype)
03  DESeq2 marker gene detection         (per celltype)
04  fGSEA pathway enrichment             (per celltype)  <- optional
05  Aggregate results + status plot      (once after all celltypes)
06  Pipeline summary                     (once -- final)
```

For each cell type, the pipeline builds a combined count matrix where every
sample appears twice -- once as "interest" (cells of this type) and once as
"other" (all remaining cells). DESeq2 tests `~ subtype + sample`, where
`subtype` captures the cell-type effect and `sample` blocks on donor.
Positive log2FoldChange = upregulated in the cell type of interest.

Step 02 includes **preflight checks** -- minimum sample counts, gene
filtering, and donor deduplication. Rejected cell types are skipped cleanly
through the rest of the pipeline.

---

## Quick Start

```bash
# 1. Install and activate
micromamba env create -f installation/EasyCM_install.yml
micromamba activate EasyCM

# 2. Place your data
#    - Count matrices in data/counts/ (cell_mtx/ + allbut_mtx/ or {CellType}.counts.csv)
#    - Sample metadata at data/sample_metadata.csv

# 3. Edit config
#    cp config/config.yaml config/my_config.yaml
#    # Edit filtering, thresholds, paths as needed

# 4. Run the full pipeline
snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml

# 5. Check results
cat results/pipeline_summary.csv | column -t -s,
```

---

## Project Layout

```
Easy_CM/
├── Snakefile                             <- Snakemake workflow (do not edit)
├── config/
│   └── config.yaml                       <- EDIT THIS
├── profiles/                             <- Snakemake execution profiles
│   ├── local/config.yaml
│   └── slurm/config.yaml
├── data/                                 <- YOUR DATA
│   ├── counts/                           <- pseudobulk count matrices
│   │   ├── cell_mtx/                     <- per-celltype interest matrices
│   │   └── allbut_mtx/                   <- per-celltype "all others" matrices
│   └── sample_metadata.csv
├── resources/
│   └── gsea_files/                       <- pathway GMT files + gene exclusion lists
│       ├── reactome_kegg.gmt.txt
│       ├── rpl_genes.csv
│       ├── rps_genes.csv
│       └── mtr_genes.csv
├── workflow/scripts/                     <- pipeline code (do not edit)
│   ├── 01_validate.R                     <- input validation
│   ├── 02_prepare_pseudobulk.R           <- matrix preparation + filtering
│   ├── 03_run_deseq.R                    <- DESeq2 marker detection
│   ├── 04_run_fgsea.R                    <- pathway enrichment
│   ├── 05_aggregate_results.R            <- cross-celltype aggregation
│   ├── 06_pipeline_summary.R             <- final summary + plots
│   └── utils/                            <- shared R utilities
│       ├── io_utils.R
│       ├── logging_utils.R
│       ├── filter_utils.R
│       └── validation_utils.R
├── installation/
│   └── EasyCM_install.yml
├── notebooks/
│   └── explore_results.ipynb
└── docs/                                 <- detailed documentation
```

| What to edit | Where |
|-------------|-------|
| Analysis parameters | `config/config.yaml` |
| Execution settings | `profiles/*/config.yaml` |
| Data files | `data/counts/` + `data/sample_metadata.csv` |

---

## Documentation

| Guide | Contents |
|-------|----------|
| **[Installation](docs/INSTALLATION.md)** | Environment setup, verification |
| **[Configuration](docs/CONFIGURATION.md)** | Config file, data formats, cell type discovery |
| **[Running](docs/RUNNING.md)** | Snakemake commands, SLURM setup, step-by-step manual run |
| **[Methods](docs/METHODS.md)** | Pseudobulk design, DESeq2, fGSEA, gene filtering |
| **[Output](docs/OUTPUT.md)** | File tree, column descriptions, status taxonomy |
| **[Troubleshooting](docs/TROUBLESHOOTING.md)** | Common errors, log reading, resource files |

---

## Pipeline Status Values

The `celltype_summary.csv` produced by step 05 categorizes each cell type
with a fine-grained status taxonomy:

| Status | Meaning |
|--------|---------|
| `success` | DESeq2 completed, marker genes found |
| `success_no_significant` | DESeq2 ran but no genes passed FDR threshold |
| `skipped_no_samples` | Zero samples after filtering or cell type not found |
| `skipped_min_cells` | Too few samples after deduplication |
| `skipped_preflight` | Too few genes passed filtering |
| `skipped` | Other skip reason -- check `error_message` column |
| `error` | DESeq2 or upstream step crashed |
| `not_run` | No log or output data found |

Step 05 generates a **status overview heatmap** (`status_overview.pdf`)
showing status and DEG counts across all cell types at a glance.
Step 06 generates a **DEG bar chart** (`pipeline_summary.pdf`) ranking
cell types by marker gene count.

