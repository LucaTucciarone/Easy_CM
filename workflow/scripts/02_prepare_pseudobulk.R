# =============================================================================
# 02_prepare_pseudobulk.R
# EasyCM -- Pseudobulk Matrix Preparation
#
# PURPOSE:
#   For a given cell type, prepare the combined count matrix and coldata
#   for DESeq2. Either builds pseudobulk from a Seurat object or loads
#   pre-made count matrices.
#
#   The key design: each sample appears TWICE in the combined matrix --
#   once as "interest" (cells of this type) and once as "other" (all other
#   cells). DESeq2 then tests ~ subtype + sample, where subtype = interest/other
#   and sample controls for donor effects.
#
# USAGE:
#   Rscript workflow/scripts/02_prepare_pseudobulk.R \
#       --config config/config.yaml --celltype Beta
#
# INPUTS:
#   - Config YAML
#   - Sample metadata CSV
#   - Either: Seurat RDS object (make_pseudobulk = true)
#   - Or: Pre-made count matrices in counts_dir (make_pseudobulk = false)
#
# OUTPUTS:
#   results/{celltype}/intermediates/coldata.csv
#   results/{celltype}/intermediates/counts_filtered.csv
#   results/{celltype}/intermediates/counts_interest.csv
#   results/{celltype}/intermediates/counts_other.csv
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(optparse)
})

# --- Resolve utils ---
.script_dir <- tryCatch({
    args     <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("--file=", "", file_arg[1])))
    } else {
        dirname(normalizePath(sys.frame(1)$ofile))
    }
}, error = function(e) getwd())

source(file.path(.script_dir, "utils", "io_utils.R"))
source(file.path(.script_dir, "utils", "logging_utils.R"))
source(file.path(.script_dir, "utils", "filter_utils.R"))


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, celltype) {

    cfg <- load_config(config_path)

    # --- Set up paths and logging ---
    results_dir <- cfg$outputs$results_dir
    inter_dir   <- file.path(results_dir, celltype, "intermediates")
    dir.create(inter_dir, recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(cfg$outputs$logs_dir, paste0(celltype, ".log"))
    logger   <- make_logger(celltype = celltype, log_file = log_file,
                            level = cfg$logging$level %||% "INFO")

    # Register skip sentinels for expected outputs
    register_skip_sentinel(inter_dir,
        c("coldata.csv", "counts_filtered.csv"), env = environment())

    logger$info("start", sprintf("Preparing pseudobulk for celltype=%s", celltype))

    # --- Load sample metadata ---
    sample_meta <- fread(cfg$inputs$sample_metadata) %>% as.data.frame()
    logger$info("load_metadata", sprintf("n_samples=%d, n_cols=%d",
                nrow(sample_meta), ncol(sample_meta)))

    # --- Apply sample filters from config ---
    filters <- cfg$filtering$sample_filters
    if (!is.null(filters) && length(filters) > 0) {
        sample_meta <- apply_subset_filters(sample_meta, filters)
        logger$info("sample_filters", sprintf("n_samples_after_filter=%d", nrow(sample_meta)))
    }

    if (nrow(sample_meta) == 0) {
        logger$skip("no_samples", "Zero samples after applying sample filters")
        return(invisible(NULL))
    }

    # --- Load count matrices ---
    if (isTRUE(cfg$steps$make_pseudobulk)) {
        # Build from Seurat object
        interest_and_other <- build_pseudobulk_from_seurat(
            cfg, celltype, sample_meta, logger)
        if (is.null(interest_and_other)) return(invisible(NULL))
        interest_counts <- interest_and_other$interest
        other_counts    <- interest_and_other$other
    } else {
        # Load pre-made matrices
        loaded <- load_premade_matrices(cfg, celltype, logger)
        if (is.null(loaded)) return(invisible(NULL))
        interest_counts <- loaded$interest
        other_counts    <- loaded$other
    }

    # --- Donor filtering: keep only samples with enough cells ---
    sample_col <- cfg$inputs$sample_column %||% "samples"
    donor_col  <- cfg$filtering$donor_col

    # Standardize sample names (replace hyphens with dots to match matrix colnames)
    sample_meta[[sample_col]] <- gsub("-", ".", sample_meta[[sample_col]])

    # Filter to samples present in the interest matrix
    available_samples <- colnames(interest_counts)
    keep_samples <- intersect(available_samples, sample_meta[[sample_col]])

    if (length(keep_samples) == 0) {
        logger$skip("no_samples", "No samples in interest matrix match metadata")
        return(invisible(NULL))
    }

    interest_counts <- interest_counts[, keep_samples, drop = FALSE]
    other_counts    <- other_counts[, intersect(colnames(other_counts), keep_samples), drop = FALSE]

    # Ensure both matrices have same samples
    shared_samples  <- intersect(colnames(interest_counts), colnames(other_counts))
    interest_counts <- interest_counts[, shared_samples, drop = FALSE]
    other_counts    <- other_counts[, shared_samples, drop = FALSE]

    logger$info("sample_match", sprintf("n_shared_samples=%d", length(shared_samples)))

    if (length(shared_samples) < 3) {
        logger$skip("too_few_samples",
            sprintf("Only %d samples with data -- need at least 3", length(shared_samples)))
        return(invisible(NULL))
    }

    # --- Donor deduplication ---
    meta_filtered <- sample_meta[sample_meta[[sample_col]] %in% shared_samples, , drop = FALSE]
    if (!is.null(donor_col) && donor_col %in% colnames(meta_filtered)) {
        meta_deduped <- dedup_donors(meta_filtered, donor_col)
        keep_after_dedup <- meta_deduped[[sample_col]]
        interest_counts <- interest_counts[, keep_after_dedup, drop = FALSE]
        other_counts    <- other_counts[, keep_after_dedup, drop = FALSE]
        logger$info("donor_dedup", sprintf("n_samples=%d->%d",
                    length(shared_samples), length(keep_after_dedup)))
        shared_samples <- keep_after_dedup
    }

    if (length(shared_samples) < 3) {
        logger$skip("too_few_samples_after_dedup",
            sprintf("Only %d samples after dedup", length(shared_samples)))
        return(invisible(NULL))
    }

    # --- Gene filtering (on interest counts only, as in OG pipeline) ---
    min_counts <- cfg$filtering$min_gene_counts %||% 5
    min_prop   <- cfg$filtering$min_gene_proportion %||% 0.25

    interest_filtered <- filter_genes(interest_counts, min_reads = min_counts, min_prop = min_prop)
    genes_to_keep     <- rownames(interest_filtered)

    logger$info("gene_filter", sprintf("n_genes=%d->%d (min_counts=%d, min_prop=%.2f)",
                nrow(interest_counts), length(genes_to_keep), min_counts, min_prop))

    if (length(genes_to_keep) < 10) {
        logger$skip("too_few_genes",
            sprintf("Only %d genes passed filter", length(genes_to_keep)))
        return(invisible(NULL))
    }

    # Subset both matrices to passing genes
    common_genes    <- intersect(genes_to_keep, rownames(other_counts))
    interest_counts <- interest_counts[common_genes, , drop = FALSE]
    other_counts    <- other_counts[common_genes, , drop = FALSE]

    # --- Build combined matrix and coldata ---
    # Prefix column names: interest_SampleA, other_SampleA
    colnames(interest_counts) <- paste0("interest_", colnames(interest_counts))
    colnames(other_counts)    <- paste0("other_", colnames(other_counts))

    combined_counts <- cbind(interest_counts, other_counts)

    # Remove samples with zero total counts
    nonzero_cols    <- colSums(combined_counts) > 0
    combined_counts <- combined_counts[, nonzero_cols, drop = FALSE]

    # Build coldata
    coldata <- data.frame(
        test    = colnames(combined_counts),
        subtype = ifelse(grepl("^interest_", colnames(combined_counts)), "interest", "other"),
        sample  = gsub("^(interest_|other_)", "", colnames(combined_counts)),
        stringsAsFactors = FALSE
    )
    rownames(coldata) <- coldata$test

    # Ensure subtype is a factor with interest as the second level
    # (so DESeq2 contrast interest vs other gives positive LFC for markers)
    coldata$subtype <- factor(coldata$subtype, levels = c("other", "interest"))
    coldata$sample  <- factor(coldata$sample)

    logger$info("combined_matrix", sprintf("n_genes=%d, n_columns=%d (interest=%d, other=%d)",
                nrow(combined_counts), ncol(combined_counts),
                sum(coldata$subtype == "interest"), sum(coldata$subtype == "other")))

    # --- Write outputs ---
    write_result(coldata, file.path(inter_dir, "coldata.csv"), sep = ",")
    write_result(combined_counts, file.path(inter_dir, "counts_filtered.csv"), sep = ",",
                 row_names = TRUE)
    write_result(interest_counts, file.path(inter_dir, "counts_interest.csv"), sep = ",",
                 row_names = TRUE)
    write_result(other_counts, file.path(inter_dir, "counts_other.csv"), sep = ",",
                 row_names = TRUE)

    logger$info("done", "Pseudobulk preparation complete")
    return(invisible(NULL))
}


# =============================================================================
# Pseudobulk from Seurat
# =============================================================================

build_pseudobulk_from_seurat <- function(cfg, celltype, sample_meta, logger) {

    suppressMessages(library(Seurat))

    seurat_path  <- cfg$inputs$seurat_object
    ct_col       <- cfg$inputs$celltype_column %||% "coarse_annot"
    sample_col   <- cfg$inputs$sample_column %||% "samples"

    logger$info("load_seurat", sprintf("path=%s", seurat_path))

    data <- try_logged(
        readRDS(seurat_path),
        logger, "load_seurat", "Seurat object loaded"
    )
    if (is.null(data)) return(NULL)

    # Clean cell type annotations (remove spaces, as in OG pipeline)
    data@meta.data[[ct_col]] <- gsub(" ", "", data@meta.data[[ct_col]])
    Idents(data) <- data@meta.data[[ct_col]]

    available_types <- unique(data@meta.data[[ct_col]])
    if (!celltype %in% available_types) {
        logger$skip("celltype_not_found",
            sprintf("'%s' not in Seurat celltypes: %s",
                    celltype, paste(available_types, collapse = ", ")))
        return(NULL)
    }

    # Subset and aggregate
    interest_seurat <- subset(data, idents = celltype)
    other_seurat    <- subset(data, idents = celltype, invert = TRUE)

    interest_bulk <- try_logged(
        AggregateExpression(interest_seurat, group.by = sample_col,
                           return.seurat = TRUE)[["RNA"]]$counts %>% as.data.frame(),
        logger, "aggregate_interest", "Interest pseudobulk created"
    )
    if (is.null(interest_bulk)) return(NULL)

    other_bulk <- try_logged(
        AggregateExpression(other_seurat, group.by = sample_col,
                           return.seurat = TRUE)[["RNA"]]$counts %>% as.data.frame(),
        logger, "aggregate_other", "Other pseudobulk created"
    )
    if (is.null(other_bulk)) return(NULL)

    logger$info("pseudobulk", sprintf("interest: %d genes x %d samples, other: %d genes x %d samples",
                nrow(interest_bulk), ncol(interest_bulk),
                nrow(other_bulk), ncol(other_bulk)))

    list(interest = interest_bulk, other = other_bulk)
}


# =============================================================================
# Load pre-made matrices
# =============================================================================

load_premade_matrices <- function(cfg, celltype, logger) {

    counts_dir <- cfg$inputs$counts_dir
    safe_name  <- gsub("[^[:alnum:]]", "_", celltype)

    # Look for interest (cell_mtx) and other (allbut_mtx) matrices
    # Support both EasyCM naming ({CellType}.counts.csv) and OG naming
    interest_path <- NULL
    other_path    <- NULL

    # Try EasyCM naming: {CellType}.counts.csv in counts_dir
    easycm_path <- file.path(counts_dir, paste0(safe_name, ".counts.csv"))
    if (file.exists(easycm_path)) {
        # Single matrix per celltype -- need to handle differently
        # For now, we need paired matrices (interest + other)
        logger$info("load_counts", sprintf("Found single matrix: %s", easycm_path))
        logger$error("load_counts",
            "EasyCM requires paired matrices (interest + other). Provide cell_mtx/ and allbut_mtx/ subdirectories, or use make_pseudobulk=true with a Seurat object.")
        return(NULL)
    }

    # Try OG naming: cell_mtx/{CellType}_persample_RNA_counts.tsv
    og_interest <- file.path(counts_dir, "cell_mtx",
                             paste0(celltype, "_persample_RNA_counts.tsv"))
    og_other    <- file.path(counts_dir, "allbut_mtx",
                             paste0(celltype, "_persample_RNA_counts.tsv"))

    # Also try with safe_name
    if (!file.exists(og_interest)) {
        og_interest <- file.path(counts_dir, "cell_mtx",
                                 paste0(safe_name, "_persample_RNA_counts.tsv"))
    }
    if (!file.exists(og_other)) {
        og_other <- file.path(counts_dir, "allbut_mtx",
                              paste0(safe_name, "_persample_RNA_counts.tsv"))
    }

    # Try paired naming: {CellType}_interest.counts.csv / {CellType}_other.counts.csv
    paired_interest <- file.path(counts_dir, paste0(safe_name, "_interest.counts.csv"))
    paired_other    <- file.path(counts_dir, paste0(safe_name, "_other.counts.csv"))

    if (file.exists(og_interest) && file.exists(og_other)) {
        interest_path <- og_interest
        other_path    <- og_other
    } else if (file.exists(paired_interest) && file.exists(paired_other)) {
        interest_path <- paired_interest
        other_path    <- paired_other
    } else {
        logger$skip("no_matrices",
            sprintf("Could not find paired count matrices for '%s' in %s", celltype, counts_dir))
        return(NULL)
    }

    logger$info("load_counts", sprintf("interest=%s, other=%s", interest_path, other_path))

    interest <- load_count_matrix(interest_path)
    other    <- load_count_matrix(other_path)

    if (is.null(interest) || is.null(other)) {
        logger$error("load_counts", "Failed to load one or both count matrices")
        return(NULL)
    }

    # First column is gene names -- set as rownames
    if (!is.numeric(interest[[1]])) {
        rownames(interest) <- interest[[1]]
        interest <- interest[, -1, drop = FALSE]
    }
    if (!is.numeric(other[[1]])) {
        rownames(other) <- other[[1]]
        other <- other[, -1, drop = FALSE]
    }

    # Standardize column names
    colnames(interest) <- gsub("-", ".", colnames(interest))
    colnames(other)    <- gsub("-", ".", colnames(other))

    logger$info("load_counts", sprintf("interest: %d genes x %d samples, other: %d genes x %d samples",
                nrow(interest), ncol(interest), nrow(other), ncol(other)))

    list(interest = interest, other = other)
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive()) {

    opt_list <- list(
        make_option("--config", type = "character",
                    default = "config/config.yaml",
                    help = "Path to config.yaml"),
        make_option("--celltype", type = "character",
                    help = "Cell type to process")
    )

    opts <- parse_args(OptionParser(option_list = opt_list))

    if (is.null(opts$celltype)) {
        stop("--celltype is required", call. = FALSE)
    }

    tryCatch({
        main(opts$config, opts$celltype)
        quit(status = 0)
    }, error = function(e) {
        message("ERROR: ", conditionMessage(e))
        quit(status = 1)
    })
}
