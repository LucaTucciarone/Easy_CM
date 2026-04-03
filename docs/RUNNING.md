# Running the Pipeline

## Running with Snakemake

Snakemake is the recommended way to run EasyCM. It handles the full DAG
of dependencies, parallelizes steps 02--04 across all cell types, and
automatically retries failed jobs.

### Run all cell types

```bash
micromamba activate EasyCM

snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml
```

This runs **all** cell types discovered from the count matrices. For 12
cell types, that is ~40 jobs total (12 x 3 per-celltype steps + 3 global
steps).

### Dry run (see what would execute)

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml -n
```

### Run specific cell types

Restrict to specific cell types by listing them in the config:

```yaml
inputs:
  celltypes:
    - Beta
    - Alpha
```

Or request specific output files:

```bash
snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml \
    results/Beta/finals/results_deseq.tsv
```

### Clean and re-run

```bash
# Remove results for one cell type
rm -rf results/Beta/

# Remove all results
rm -rf results/ logs/

# Then re-run
snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml
```

---

## Running on SLURM

For cluster execution, use the SLURM profile. Edit
`profiles/slurm/config.yaml` to set your account, partition, and QOS:

```yaml
# profiles/slurm/config.yaml (edit these)
default-resources:
  slurm_account: your_account
  slurm_partition: your_partition
  slurm_qos: your_qos
  mem_mb: 8000
  time: 240      # minutes
  cpus: 1
```

Then run:

```bash
snakemake --profile $PWD/profiles/slurm \
    --config pipeline_config=config/config.yaml
```

> **Note:** Use `$PWD/profiles/slurm` (absolute path) for the SLURM profile.
> SLURM jobs inherit the working directory, so all paths resolve correctly.

---

## Running Step by Step

For debugging or educational purposes, you can run each script manually.
All commands assume you are in the pipeline root directory.

```bash
micromamba activate EasyCM

CONFIG="config/config.yaml"
CELLTYPE="Beta"

# Step 01 -- Validate (once per run)
Rscript workflow/scripts/01_validate.R --config "$CONFIG"

# Step 02 -- Prepare pseudobulk matrices
Rscript workflow/scripts/02_prepare_pseudobulk.R \
    --config "$CONFIG" --celltype "$CELLTYPE"

# Step 03 -- DESeq2 marker detection
Rscript workflow/scripts/03_run_deseq.R \
    --config "$CONFIG" --celltype "$CELLTYPE"

# Step 04 -- fGSEA pathway enrichment (optional)
Rscript workflow/scripts/04_run_fgsea.R \
    --config "$CONFIG" --celltype "$CELLTYPE"

# Step 05 -- Aggregate results (once, after all cell types complete)
Rscript workflow/scripts/05_aggregate_results.R --config "$CONFIG"

# Step 06 -- Pipeline summary (once, after aggregation)
Rscript workflow/scripts/06_pipeline_summary.R --config "$CONFIG"
```

### Timing expectations

On a local machine (64 cores):

| Step | Time per cell type | Notes |
|------|-------------------|-------|
| 01 Validate | ~2 seconds | Once per run |
| 02 Pseudobulk | ~1 second | Pre-made matrices; longer with Seurat |
| 03 DESeq2 | 1--5 minutes | Depends on gene count and sample size |
| 04 fGSEA | ~5 seconds | After DESeq2 results are available |
| 05 Aggregate | ~3 seconds | Once after all cell types |
| 06 Summary | ~2 seconds | Once after aggregation |

A full run on 12 cell types typically completes in 10--20 minutes
sequentially, or 3--5 minutes with Snakemake parallelization.
