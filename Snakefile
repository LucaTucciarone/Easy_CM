# =============================================================================
# EasyCM -- Snakefile
# Cell-type marker gene discovery via pseudobulk DESeq2 + fGSEA
# =============================================================================

import yaml
import os
import re


# -----------------------------------------------------------------------------
# Config and helpers
# -----------------------------------------------------------------------------

PIPELINE_CONFIG = config.get("pipeline_config", "config/config.yaml")

with open(PIPELINE_CONFIG) as f:
    cfg = yaml.safe_load(f)

RESULTS_DIR = cfg["outputs"]["results_dir"]
LOGS_DIR    = cfg["outputs"]["logs_dir"]
SCRIPTS_DIR = "workflow/scripts"


def sanitize(name):
    """Replace non-alphanumeric chars with underscores for safe file paths."""
    return re.sub(r"[^A-Za-z0-9]", "_", name)


# --- Discover cell types ---
# If celltypes list is specified in config, use that.
# Otherwise, discover from count matrix files in counts_dir.
if cfg["inputs"].get("celltypes"):
    CELLTYPES = cfg["inputs"]["celltypes"]
else:
    counts_dir = cfg["inputs"]["counts_dir"]
    if cfg["steps"].get("make_pseudobulk", False):
        # When building from Seurat, celltypes must be specified in config
        raise ValueError(
            "steps.make_pseudobulk = true requires inputs.celltypes to be set. "
            "List the cell types to analyze in config.yaml."
        )
    else:
        # Discover from {CellType}.counts.csv files
        count_files = [f for f in os.listdir(counts_dir)
                       if f.endswith(".counts.csv") or f.endswith(".counts.csv.gz")]
        CELLTYPES = sorted(set(
            re.sub(r"\.counts\.csv(\.gz)?$", "", f) for f in count_files
        ))

        # Also check for OG pipeline layout (cell_mtx/ subdirectory)
        if not CELLTYPES:
            cell_mtx_dir = os.path.join(counts_dir, "cell_mtx")
            if os.path.isdir(cell_mtx_dir):
                tsv_files = [f for f in os.listdir(cell_mtx_dir)
                             if f.endswith("_persample_RNA_counts.tsv")]
                CELLTYPES = sorted(set(
                    f.replace("_persample_RNA_counts.tsv", "") for f in tsv_files
                ))

if not CELLTYPES:
    raise ValueError("No cell types found. Set inputs.celltypes in config or provide count matrices.")

# Map sanitized names back to originals
SAFE_MAP = {sanitize(ct): ct for ct in CELLTYPES}
CELLTYPES_SAFE = list(SAFE_MAP.keys())


def original_celltype(wildcards):
    """Recover original cell type name from sanitized wildcard."""
    return SAFE_MAP[wildcards.celltype]


def inter(celltype, filename):
    return os.path.join(RESULTS_DIR, celltype, "intermediates", filename)


def final(celltype, filename):
    return os.path.join(RESULTS_DIR, celltype, "finals", filename)


# -----------------------------------------------------------------------------
# Target rule
# -----------------------------------------------------------------------------

rule all:
    input:
        os.path.join(RESULTS_DIR, "pipeline_summary.pdf")


# -----------------------------------------------------------------------------
# Rule 01 -- Validate inputs (once)
# -----------------------------------------------------------------------------

rule validate:
    input:
        config = PIPELINE_CONFIG,
        meta   = cfg["inputs"]["sample_metadata"],
    output:
        touch(os.path.join(LOGS_DIR, "01_validate.done"))
    log:
        os.path.join(LOGS_DIR, "01_validate.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/01_validate.R \
            --config {input.config} \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 02 -- Prepare pseudobulk matrices (per cell type)
# -----------------------------------------------------------------------------

rule prepare_pseudobulk:
    input:
        done   = os.path.join(LOGS_DIR, "01_validate.done"),
        config = PIPELINE_CONFIG,
        meta   = cfg["inputs"]["sample_metadata"],
    output:
        coldata = inter("{celltype}", "coldata.csv"),
        counts  = inter("{celltype}", "counts_filtered.csv"),
    params:
        celltype = original_celltype
    log:
        os.path.join(LOGS_DIR, "{celltype}_02.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/02_prepare_pseudobulk.R \
            --config {input.config} \
            --celltype "{params.celltype}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 03 -- DESeq2 marker gene detection (per cell type)
# -----------------------------------------------------------------------------

rule run_deseq:
    input:
        config  = PIPELINE_CONFIG,
        coldata = inter("{celltype}", "coldata.csv"),
        counts  = inter("{celltype}", "counts_filtered.csv"),
    output:
        results = final("{celltype}", "results_deseq.tsv"),
    params:
        celltype = original_celltype
    log:
        os.path.join(LOGS_DIR, "{celltype}_03.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/03_run_deseq.R \
            --config {input.config} \
            --celltype "{params.celltype}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 04 -- fGSEA pathway enrichment (per cell type, optional)
# -----------------------------------------------------------------------------

rule run_fgsea:
    input:
        config  = PIPELINE_CONFIG,
        results = final("{celltype}", "results_deseq.tsv"),
    output:
        touch(os.path.join(RESULTS_DIR, "{celltype}", "fgsea.done"))
    params:
        celltype = original_celltype
    log:
        os.path.join(LOGS_DIR, "{celltype}_04.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/04_run_fgsea.R \
            --config {input.config} \
            --celltype "{params.celltype}" \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 05 -- Aggregate results (once)
# -----------------------------------------------------------------------------

rule aggregate:
    input:
        config = PIPELINE_CONFIG,
        done   = expand(
            os.path.join(RESULTS_DIR, "{celltype}", "fgsea.done"),
            celltype=CELLTYPES_SAFE
        ),
    output:
        summary = os.path.join(RESULTS_DIR, "summary", "celltype_summary.csv"),
    log:
        os.path.join(LOGS_DIR, "05_aggregate.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/05_aggregate_results.R \
            --config {input.config} \
            > {log} 2>&1
        """


# -----------------------------------------------------------------------------
# Rule 06 -- Pipeline summary (once, final)
# -----------------------------------------------------------------------------

rule pipeline_summary:
    input:
        config  = PIPELINE_CONFIG,
        summary = os.path.join(RESULTS_DIR, "summary", "celltype_summary.csv"),
    output:
        os.path.join(RESULTS_DIR, "pipeline_summary.pdf")
    log:
        os.path.join(LOGS_DIR, "06_pipeline_summary.log")
    shell:
        """
        Rscript {SCRIPTS_DIR}/06_pipeline_summary.R \
            --config {input.config} \
            > {log} 2>&1
        """
