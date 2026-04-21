# Running the Pipeline

All invocations go through the `EasyDE` micromamba environment, which provides DESeq2,
apeglm, fGSEA, and the rest of the R stack.

```bash
MAMBA=/home/luca/bin/micromamba
ENV=EasyDE
```

---

## Run with Snakemake

### Full run (all cell types)

```bash
$MAMBA run -n $ENV snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml
```

Rules 02–04 run in parallel across cell types. For the PankBase 12-cell-type run this is
~40 jobs total (12 × 3 per-celltype + validate + aggregate + summary).

### Dry run

```bash
$MAMBA run -n $ENV snakemake -n --config pipeline_config=config/config.yaml
```

### Restrict to specific cell types

Either list them in config:

```yaml
inputs:
  celltypes: [Beta, Alpha]
```

…or request specific outputs:

```bash
$MAMBA run -n $ENV snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml \
    results/Beta/finals/results_deseq.tsv
```

### Clean and re-run

```bash
rm -rf results/Beta/              # one cell type
rm -rf results/ logs/             # everything

$MAMBA run -n $ENV snakemake --profile profiles/local \
    --config pipeline_config=config/config.yaml
```

---

## Run step by step (manual / debugging)

```bash
CONFIG=config/config.yaml
CELLTYPE=Beta
RUN="$MAMBA run -n $ENV Rscript"

$RUN workflow/scripts/01_validate.R            --config $CONFIG
$RUN workflow/scripts/02_prepare_pseudobulk.R  --config $CONFIG --celltype $CELLTYPE
$RUN workflow/scripts/03_run_deseq.R           --config $CONFIG --celltype $CELLTYPE
$RUN workflow/scripts/04_run_fgsea.R           --config $CONFIG --celltype $CELLTYPE
$RUN workflow/scripts/05_aggregate_results.R   --config $CONFIG    # after all celltypes
$RUN workflow/scripts/06_pipeline_summary.R    --config $CONFIG
```

---

## Timing (PankBase, 191 samples, 12 cell types, ophelia)

| Rule | Per cell type | Notes |
|------|---------------|-------|
| 02 Pseudobulk | ~1 s | pre-made matrices |
| 03 DESeq2 + apeglm | 1–5 min | gene-count dependent |
| 04 fGSEA | ~5 s | |
| 05 + 06 | ~5 s total | once |

Sequential: ~15 min. Parallel (local profile, 12 jobs): ~5 min.

---

## SLURM

`profiles/slurm/config.yaml` is provided but **not tested on ophelia** (no `sbatch`).
To use on a SLURM cluster, edit the profile to set your account/partition/QOS and invoke
with `--profile $PWD/profiles/slurm` (absolute path).
