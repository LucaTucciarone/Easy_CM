# Installation

## Option 1: Conda/micromamba (recommended)

```bash
cd /path/to/Easy_CM
micromamba env create -f installation/EasyCM_install.yml
micromamba activate EasyCM
```

This installs R, DESeq2, fgsea, Snakemake, and all dependencies.
Takes ~5-10 minutes on first run.

> Using conda/mamba? Replace `micromamba` with `conda` or `mamba`.

### Verify

```bash
Rscript -e "cat('R ok\n')"        # should print "R ok"
snakemake --version                 # should print 7.x
```

---

## Option 2: Install micromamba first (HPC / fresh system)

```bash
# Install micromamba (user-level, no root needed)
wget -O ~/micromamba.tar.bz2 https://micro.mamba.pm/api/micromamba/linux-64/latest
tar -xvjf ~/micromamba.tar.bz2 -C ~/
rm ~/micromamba.tar.bz2

# Initialize shell integration
~/bin/micromamba shell init -s bash -r ~/micromamba
source ~/.bashrc

# Configure channels
micromamba config append channels bioconda
micromamba config append channels conda-forge
micromamba config set channel_priority strict

# Then proceed with Option 1
micromamba env create -f installation/EasyCM_install.yml
```

---

## Optional: Seurat (for pseudobulk from Seurat objects)

If you want EasyCM to build pseudobulk matrices directly from a Seurat
`.rds` object (i.e. `steps.make_pseudobulk: true`), install Seurat:

```bash
micromamba install -n EasyCM -c conda-forge r-seurat
```

Or uncomment the `r-seurat` line in `installation/EasyCM_install.yml`
before creating the environment. Seurat is **not required** if you provide
pre-made count matrices.

---

## Key Dependencies

| Package | Purpose |
|---------|---------|
| R (r-base) | Runtime for all analysis scripts |
| DESeq2 | Differential expression (negative binomial GLM) |
| fgsea | Fast gene set enrichment analysis |
| Snakemake | Pipeline orchestration |
| r-yaml, r-optparse | Config parsing, CLI arguments |
| r-data.table | Fast file I/O |
| r-tidyverse | Data wrangling |
| r-ggplot2, r-ggrepel | Visualization |
| r-seurat (optional) | Pseudobulk from Seurat objects |
