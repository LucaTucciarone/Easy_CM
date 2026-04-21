# =============================================================================
# 02_prepare_pseudobulk.R
# EasyCM -- Pseudobulk Matrix Preparation (consumes EasyPseudobulk output)
#
# PURPOSE:
#   For a given cell type, assemble the combined count matrix and coldata
#   for DESeq2 from pre-made pseudobulk matrices produced by EasyPseudobulk.
#
#   Each donor appears TWICE in the combined matrix -- once as "interest"
#   (cells of this type) and once as "other" (all other cells). DESeq2 then
#   tests ~ subtype + sample, where subtype = interest/other and sample
#   controls for donor effects.
#
#   NOTE: this script used to build pseudobulks directly from a Seurat RDS.
#   That path was extracted into the standalone EasyPseudobulk notebook so
#   every tool in the Easy suite can share one upstream pseudobulk step.
#   EasyCM now ONLY consumes the counts/easycm layout EasyPseudobulk writes.
#
# USAGE:
#   Rscript workflow/scripts/02_prepare_pseudobulk.R \
#       --config config/config.yaml --celltype Beta
#
# INPUTS:
#   - Config YAML (inputs.counts_dir must point at EasyPseudobulk's
#     {output_dir}/counts/easycm/ directory)
#   - Sample metadata tsv/csv (EasyPseudobulk's donor_metadata.tsv works
#     out of the box once the sample_column is set to the donor column)
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

    # --- Derive missing filter columns from sample names ---
    # Some datasets encode treatments in sample names (e.g., "HP.21070__DHT_10nM").
    # If config filters reference a column that doesn't exist in metadata but sample
    # names contain "__", derive the column: text after last "__" = treatment value,
    # no "__" = "no_treatment".
    filters <- cfg$filtering$sample_filters
    sample_col <- cfg$inputs$sample_column %||% "donor_accession"
    if (!is.null(filters) && length(filters) > 0) {
        for (filt_col in names(filters)) {
            if (!filt_col %in% colnames(sample_meta) && any(grepl("__", sample_meta[[sample_col]]))) {
                sample_meta[[filt_col]] <- ifelse(
                    grepl("__", sample_meta[[sample_col]]),
                    sub(".*__", "", sample_meta[[sample_col]]),
                    "no_treatment"
                )
                logger$info("derive_column",
                    sprintf("Derived '%s' from sample names (%d unique values)",
                            filt_col, length(unique(sample_meta[[filt_col]]))))
            }
        }
        sample_meta <- apply_subset_filters(sample_meta, filters)
        logger$info("sample_filters", sprintf("n_samples_after_filter=%d", nrow(sample_meta)))
    }

    if (nrow(sample_meta) == 0) {
        logger$skip("no_samples", "Zero samples after applying sample filters")
        return(invisible(NULL))
    }

    # --- Load pre-made count matrices from EasyPseudobulk ---
    loaded <- load_premade_matrices(cfg, celltype, logger)
    if (is.null(loaded)) return(invisible(NULL))
    interest_counts <- loaded$interest
    other_counts    <- loaded$other

    # --- Sample matching and optional donor dedup ---
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

    # --- Min-cells-per-sample filter ---
    # EasyPseudobulk writes cell_counts.tsv alongside sample_metadata.tsv.
    # Use it to drop samples with too few cells of this type.
    min_cells <- cfg$filtering$min_cells_per_sample %||% 0
    if (min_cells > 0) {
        cc_path <- file.path(dirname(cfg$inputs$sample_metadata), "cell_counts.tsv")
        if (file.exists(cc_path)) {
            cell_counts <- fread(cc_path) %>% as.data.frame()
            # Standardize sample names
            cc_sample_col <- colnames(cell_counts)[1]
            cell_counts[[cc_sample_col]] <- gsub("-", ".", cell_counts[[cc_sample_col]])

            # Find column matching this cell type (try exact, then unsanitized)
            ct_col <- NULL
            for (cand in c(celltype, gsub("_", " ", celltype))) {
                if (cand %in% colnames(cell_counts)) { ct_col <- cand; break }
            }

            if (!is.null(ct_col)) {
                passing <- cell_counts[[cc_sample_col]][
                    !is.na(cell_counts[[ct_col]]) & cell_counts[[ct_col]] >= min_cells
                ]
                before_n <- length(shared_samples)
                shared_samples  <- intersect(shared_samples, passing)
                interest_counts <- interest_counts[, shared_samples, drop = FALSE]
                other_counts    <- other_counts[, shared_samples, drop = FALSE]
                logger$info("min_cells_filter",
                    sprintf("min_cells>=%d for '%s': %d->%d samples",
                            min_cells, ct_col, before_n, length(shared_samples)))
            } else {
                logger$warn("min_cells_filter",
                    sprintf("Cell count column for '%s' not found in %s -- skipping",
                            celltype, cc_path))
            }
        } else {
            logger$warn("min_cells_filter",
                sprintf("cell_counts.tsv not found at %s -- skipping min-cells filter", cc_path))
        }
    }

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
# Load pre-made matrices (EasyPseudobulk counts/easycm layout)
# =============================================================================
#
# Expects OG-style layout produced by EasyPseudobulk:
#   {counts_dir}/cell_mtx/{celltype}_persample_RNA_counts.tsv   # interest
#   {counts_dir}/allbut_mtx/{celltype}_persample_RNA_counts.tsv # other
#
# Cell-type names are filesystem-safe (gsub("[^A-Za-z0-9]+", "_", label)).
# This function tries the literal name first, then the safe name.
# =============================================================================

load_premade_matrices <- function(cfg, celltype, logger) {

    counts_dir <- cfg$inputs$counts_dir
    safe_name  <- gsub("[^A-Za-z0-9]+", "_", celltype)

    candidates <- list(
        list(interest = file.path(counts_dir, "cell_mtx",
                                  paste0(celltype,  "_persample_RNA_counts.tsv")),
             other    = file.path(counts_dir, "allbut_mtx",
                                  paste0(celltype,  "_persample_RNA_counts.tsv"))),
        list(interest = file.path(counts_dir, "cell_mtx",
                                  paste0(safe_name, "_persample_RNA_counts.tsv")),
             other    = file.path(counts_dir, "allbut_mtx",
                                  paste0(safe_name, "_persample_RNA_counts.tsv")))
    )

    interest_path <- NULL
    other_path    <- NULL
    for (cand in candidates) {
        if (file.exists(cand$interest) && file.exists(cand$other)) {
            interest_path <- cand$interest
            other_path    <- cand$other
            break
        }
    }

    if (is.null(interest_path)) {
        logger$skip("no_matrices",
            sprintf("No paired count matrices for '%s' under %s/{cell_mtx,allbut_mtx}/",
                    celltype, counts_dir))
        return(NULL)
    }

    logger$info("load_counts", sprintf("interest=%s, other=%s",
                basename(interest_path), basename(other_path)))

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

    logger$info("load_counts",
        sprintf("interest: %d genes x %d samples, other: %d genes x %d samples",
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
