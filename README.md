# EasyCM

Snakemake pipeline for cell-type marker gene discovery from scRNA-seq pseudobulk data.
Per cell type: paired DESeq2 (`~ subtype + sample`) vs all other cells, apeglm LFC shrinkage,
then fGSEA against Reactome+KEGG. Consumes pre-made pseudobulk matrices from
[EasyPseudobulk](https://github.com/).

All R scripts must run inside the `EasyDE` micromamba environment (provides DESeq2, apeglm, fGSEA).

---

## Pipeline Overview

Six rules, single `{celltype}` wildcard. Rules 02–04 run in parallel across all cell types.

```
01  Validate config + inputs           (once)
02  Prepare pseudobulk matrices        (per celltype)
03  DESeq2 + apeglm LFC shrinkage      (per celltype)
04  fGSEA pathway enrichment           (per celltype, optional)
05  Aggregate results + status plot    (once, after all celltypes)
06  Pipeline summary + DEG bar chart   (once, final)
```

For each cell type, rule 02 builds a combined count matrix in which every sample appears
twice — once as **interest** (cells of that type) and once as **other** (all remaining cells).
DESeq2 then tests `~ subtype + sample` to isolate the cell-type effect while blocking on donor.
Positive `log2FoldChange` = upregulated in the cell type of interest.

---

## Quick Start

```bash
# 1. Install
micromamba env create -f installation/EasyCM_install.yml

# 2. Place inputs
#    data/counts/cell_mtx/{CellType}_persample_RNA_counts.tsv
#    data/counts/allbut_mtx/{CellType}_persample_RNA_counts.tsv
#    data/sample_metadata.csv

# 3. Edit config
cp config/config.yaml config/my_config.yaml

# 4. Dry run, then run
micromamba run -n EasyDE snakemake -n --config pipeline_config=config/my_config.yaml
micromamba run -n EasyDE snakemake --profile profiles/local \
    --config pipeline_config=config/my_config.yaml

# 5. Inspect
column -t -s, results/pipeline_summary.csv
```

---

## Project Layout

```
EasyCM/
├── Snakefile
├── config/config.yaml                <- edit this
├── profiles/{local,slurm}/config.yaml
├── data/
│   ├── counts/{cell_mtx,allbut_mtx}/ <- pseudobulk TSVs from EasyPseudobulk
│   └── sample_metadata.csv
├── resources/gsea_files/             <- GMT + gene exclusion lists
├── workflow/scripts/
│   ├── 01_validate.R
│   ├── 02_prepare_pseudobulk.R
│   ├── 03_run_deseq.R
│   ├── 04_run_fgsea.R
│   ├── 05_aggregate_results.R
│   ├── 06_pipeline_summary.R
│   └── utils/{io,logging,filter,validation}_utils.R
├── installation/EasyCM_install.yml
└── docs/
```

| Edit | Where |
|------|-------|
| Analysis parameters | `config/config.yaml` |
| Execution settings | `profiles/*/config.yaml` |
| Inputs | `data/counts/`, `data/sample_metadata.csv` |

---

## Documentation

| Guide | Contents |
|-------|----------|
| [Installation](docs/INSTALLATION.md) | Environment setup |
| [Configuration](docs/CONFIGURATION.md) | Config reference, input formats |
| [Running](docs/RUNNING.md) | Snakemake invocation, manual step-by-step |
| [Methods](docs/METHODS.md) | Paired design, DESeq2, apeglm, fGSEA |
| [Output](docs/OUTPUT.md) | File tree, column definitions, status taxonomy |
| [Troubleshooting](docs/TROUBLESHOOTING.md) | Common errors, logs |

---

## Pipeline Status Values

Rule 05 writes `celltype_summary.csv` with a status per cell type:

| Status | Meaning |
|--------|---------|
| `success` | DESeq2 completed, marker genes found |
| `success_no_significant` | DESeq2 ran, no genes passed FDR |
| `skipped_no_samples` | Cell type absent or zero samples after filtering |
| `skipped_min_cells` | Too few samples after donor deduplication |
| `skipped_preflight` | Too few genes passed filtering |
| `skipped` | Other skip — see `error_message` |
| `error` | DESeq2 or upstream crashed |
| `not_run` | No log or output found |

Rule 05 emits `status_overview.pdf` (status + DEG counts heatmap). Rule 06 emits
`pipeline_summary.pdf` (DEG bar chart, ordered by marker count).
