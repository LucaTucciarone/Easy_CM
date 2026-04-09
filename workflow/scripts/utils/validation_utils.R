# =============================================================================
# validation_utils.R
# EasyCM -- Config and input validation
#
# Pure validation functions with no side effects. Each returns a named list
# with $errors and $warnings character vectors.
# =============================================================================


#' Quick structural validation of a count matrix file
#'
#' Reads only the first few rows to check structure.
#'
#' @param path          Path to the count matrix file
#' @param celltype_name Name of the cell type (for error messages)
#' @return Named list: $errors, $warnings, $n_samples
validate_count_matrix <- function(path, celltype_name = "") {

    errors   <- character(0)
    warnings <- character(0)
    n_samples <- NA_integer_

    mat <- tryCatch(
        suppressWarnings(data.table::fread(path, nrows = 5, check.names = FALSE)),
        error = function(e) NULL
    )

    if (is.null(mat)) {
        errors <- c(errors, sprintf("celltype '%s': cannot read count matrix: %s",
                                    celltype_name, path))
        return(list(errors = errors, warnings = warnings, n_samples = n_samples))
    }

    if (ncol(mat) < 2) {
        errors <- c(errors, sprintf("celltype '%s': count matrix has fewer than 2 columns",
                                    celltype_name))
        return(list(errors = errors, warnings = warnings, n_samples = n_samples))
    }

    if (nrow(mat) == 0) {
        errors <- c(errors, sprintf("celltype '%s': count matrix has 0 rows (no genes)",
                                    celltype_name))
        return(list(errors = errors, warnings = warnings, n_samples = n_samples))
    }

    sample_cols <- colnames(mat)[-1]
    n_samples   <- length(sample_cols)
    dup_cols    <- sample_cols[duplicated(sample_cols)]
    if (length(dup_cols) > 0) {
        errors <- c(errors, sprintf("celltype '%s': duplicate column names: %s",
                                    celltype_name, paste(unique(dup_cols), collapse = ", ")))
    }

    if (is.numeric(mat[[1]])) {
        warnings <- c(warnings, sprintf("celltype '%s': first column is numeric (expected gene names)",
                                        celltype_name))
    }

    non_numeric <- character(0)
    for (j in 2:ncol(mat)) {
        if (!is.numeric(mat[[j]])) {
            non_numeric <- c(non_numeric, colnames(mat)[j])
        }
    }
    if (length(non_numeric) > 0) {
        errors <- c(errors, sprintf("celltype '%s': non-numeric data columns: %s",
                                    celltype_name, paste(head(non_numeric, 5), collapse = ", ")))
    }

    list(errors = errors, warnings = warnings, n_samples = n_samples)
}


#' Validate EasyCM config.yaml structure, paths, and value ranges
#'
#' @param cfg  Parsed config list from yaml::read_yaml()
#' @return Named list: $errors, $warnings
validate_config <- function(cfg) {

    errors   <- character(0)
    warnings <- character(0)

    cat("\n-- Validating config.yaml ------------------------------------------\n")

    # --- Required top-level sections ---
    required_sections <- c("inputs", "outputs", "filtering", "deseq2", "logging")
    for (section in required_sections) {
        if (is.null(cfg[[section]])) {
            errors <- c(errors, sprintf("  [ERROR] Missing required section: '%s'", section))
        }
    }

    # --- Input file paths ---
    if (!is.null(cfg$inputs)) {

        # counts_dir points at EasyPseudobulk's counts/easycm/ output
        has_counts <- !is.null(cfg$inputs$counts_dir) &&
                      dir.exists(cfg$inputs$counts_dir)

        if (!has_counts) {
            errors <- c(errors, sprintf("  [ERROR] counts_dir not found: '%s'",
                                        cfg$inputs$counts_dir %||% "<unset>"))
        } else {
            cat(sprintf("  [OK] counts_dir: %s\n", cfg$inputs$counts_dir))
        }

        if (!is.null(cfg$inputs$sample_metadata) && nzchar(cfg$inputs$sample_metadata)) {
            if (!file.exists(cfg$inputs$sample_metadata)) {
                errors <- c(errors, sprintf("  [ERROR] sample_metadata not found: '%s'",
                                            cfg$inputs$sample_metadata))
            } else {
                cat(sprintf("  [OK] sample_metadata: %s\n", cfg$inputs$sample_metadata))
            }
        } else {
            errors <- c(errors, "  [ERROR] inputs.sample_metadata is required")
        }
    }

    # --- Steps flags must be logical ---
    if (!is.null(cfg$steps)) {
        for (flag in names(cfg$steps)) {
            if (!is.logical(cfg$steps[[flag]])) {
                errors <- c(errors, sprintf("  [ERROR] steps.%s must be true or false", flag))
            }
        }
    }

    # --- DESeq2 thresholds ---
    if (!is.null(cfg$deseq2)) {
        fdr <- cfg$deseq2$fdr_threshold
        if (!is.null(fdr) && (!is.numeric(fdr) || fdr <= 0 || fdr >= 1)) {
            errors <- c(errors, sprintf("  [ERROR] deseq2.fdr_threshold must be between 0 and 1, got: %s", fdr))
        }
        lfc <- cfg$deseq2$lfc_threshold
        if (!is.null(lfc) && (!is.numeric(lfc) || lfc < 0)) {
            errors <- c(errors, sprintf("  [ERROR] deseq2.lfc_threshold must be >= 0, got: %s", lfc))
        }
    }

    # --- fGSEA settings ---
    if (!is.null(cfg$fgsea) && isTRUE(cfg$steps$fgsea)) {

        if (!is.null(cfg$fgsea$gene_sets_file)) {
            if (!file.exists(cfg$fgsea$gene_sets_file)) {
                errors <- c(errors, sprintf("  [ERROR] fgsea.gene_sets_file not found: '%s'",
                                            cfg$fgsea$gene_sets_file))
            } else {
                cat(sprintf("  [OK] gene_sets_file: %s\n", cfg$fgsea$gene_sets_file))
            }
        }

        if (!is.null(cfg$fgsea$exclude_gene_lists)) {
            for (gene_list in cfg$fgsea$exclude_gene_lists) {
                if (!file.exists(gene_list)) {
                    errors <- c(errors, sprintf("  [ERROR] exclude_gene_list not found: '%s'", gene_list))
                } else {
                    cat(sprintf("  [OK] exclude_gene_list: %s\n", gene_list))
                }
            }
        }

        min_s <- cfg$fgsea$min_gene_set_size
        max_s <- cfg$fgsea$max_gene_set_size
        if (!is.null(min_s) && !is.null(max_s) && min_s >= max_s) {
            errors <- c(errors, sprintf(
                "  [ERROR] fgsea.min_gene_set_size (%d) must be < max_gene_set_size (%d)",
                min_s, max_s))
        }

        valid_stats <- c("stat", "lfc_pvalue")
        if (!is.null(cfg$fgsea$ranking_stat) && !cfg$fgsea$ranking_stat %in% valid_stats) {
            errors <- c(errors, sprintf(
                "  [ERROR] fgsea.ranking_stat must be one of: %s. Got: '%s'",
                paste(valid_stats, collapse = ", "), cfg$fgsea$ranking_stat))
        }
    }

    # --- Logging level ---
    if (!is.null(cfg$logging$level)) {
        valid_levels <- c("DEBUG", "INFO", "WARN", "ERROR")
        if (!toupper(cfg$logging$level) %in% valid_levels) {
            errors <- c(errors, sprintf(
                "  [ERROR] logging.level must be one of: %s. Got: '%s'",
                paste(valid_levels, collapse = ", "), cfg$logging$level))
        }
    }

    list(errors = errors, warnings = warnings)
}


#' Discover cell types from count matrix files in a directory
#'
#' Supports two layouts:
#'   1. EasyCM: {CellType}.counts.csv directly in counts_dir
#'   2. OG pipeline: cell_mtx/{CellType}_persample_RNA_counts.tsv
#'
#' @param counts_dir  Path to directory with count matrices
#' @return Character vector of cell type names
discover_celltypes <- function(counts_dir) {

    # Try EasyCM naming first
    files <- list.files(counts_dir, pattern = "\\.counts\\.csv(\\.gz)?$", full.names = FALSE)
    celltypes <- sub("\\.counts\\.csv(\\.gz)?$", "", files)

    # Fall back to OG pipeline naming (cell_mtx/ subdirectory)
    if (length(celltypes) == 0) {
        cell_mtx_dir <- file.path(counts_dir, "cell_mtx")
        if (dir.exists(cell_mtx_dir)) {
            files <- list.files(cell_mtx_dir, pattern = "_persample_RNA_counts\\.tsv$",
                                full.names = FALSE)
            celltypes <- sub("_persample_RNA_counts\\.tsv$", "", files)
        }
    }

    sort(unique(celltypes))
}
