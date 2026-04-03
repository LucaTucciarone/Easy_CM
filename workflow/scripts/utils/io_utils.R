# =============================================================================
# io_utils.R
# Input/output helpers for EasyCM
#
# General-purpose functions for reading and writing data.
# No project-specific logic lives here.
# =============================================================================

suppressMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
})


#' Load a count matrix from a local file, return NULL on failure
#'
#' @param path  Path to the count matrix file
#' @return data.frame or NULL if loading failed
load_count_matrix <- function(path) {

    result <- tryCatch(
        withCallingHandlers(
            fread(path, check.names = FALSE) %>% as.data.frame(),
            warning = function(w) {
                if (grepl("column names but the data has", conditionMessage(w))) {
                    invokeRestart("muffleWarning")
                }
            }
        ),
        error = function(e) {
            message(sprintf("  [skip] Could not load matrix from: %s\n  Reason: %s",
                            path, e$message))
            return(NULL)
        }
    )

    return(result)
}


#' Build a celltype_summary row for one cell type (success or failure)
#'
#' @param celltype        Cell type being analyzed
#' @param status          One of: "success", "success_no_significant", "skipped_*", "error", "not_run"
#' @param error_message   Human-readable reason for skip/error (or NA)
#' @param n_interest      Number of "interest" samples (or NA)
#' @param n_other         Number of "other" samples (or NA)
#' @param n_samples_total Total samples (or NA)
#' @param formula_used    DESeq2 formula string (or NA)
#' @param n_degs          Number of DEGs (or NA)
#' @param n_pathways      Number of significant fGSEA pathways (or NA)
#' @return Single-row data.frame
make_summary_row <- function(celltype,
                             status,
                             error_message   = NA,
                             n_interest      = NA,
                             n_other         = NA,
                             n_samples_total = NA,
                             formula_used    = NA,
                             n_degs          = NA,
                             n_pathways      = NA) {
    data.frame(
        celltype        = celltype,
        status          = status,
        n_interest      = n_interest,
        n_other         = n_other,
        n_samples_total = n_samples_total,
        formula_used    = formula_used,
        n_degs          = n_degs,
        n_pathways      = n_pathways,
        error_message   = error_message,
        stringsAsFactors = FALSE
    )
}


#' Write a TSV or CSV result file, creating the output directory if needed
#'
#' @param df        data.frame to write
#' @param path      Output file path
#' @param sep       Separator: "\t" for TSV (default), "," for CSV
#' @param row_names Whether to include row names (default FALSE)
write_result <- function(df, path, sep = "\t", row_names = FALSE) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    write.table(df, path, sep = sep, quote = TRUE, row.names = row_names)
    message(sprintf("  [write] %s", path))
}


#' Load EasyCM pipeline config
#'
#' @param config_path  Path to config.yaml
#' @return Parsed config list
load_config <- function(config_path) {
    cfg <- suppressWarnings(yaml::read_yaml(config_path))
    return(cfg)
}


# =============================================================================
# Snakemake skip sentinel helpers
# =============================================================================

#' Register an on.exit handler that writes skip sentinels for missing outputs
#'
#' @param output_dir  Directory where outputs are written
#' @param filenames   Character vector of expected output filenames
#' @param env         Environment to register on.exit in (default: parent frame)
register_skip_sentinel <- function(output_dir, filenames, env = parent.frame()) {
    do.call(on.exit, list(substitute({
        for (.sentinel_f in .sentinel_files) {
            .sentinel_p <- file.path(.sentinel_dir, .sentinel_f)
            if (!file.exists(.sentinel_p)) {
                dir.create(dirname(.sentinel_p), recursive = TRUE, showWarnings = FALSE)
                writeLines("skipped=TRUE", .sentinel_p)
            }
        }
    }, list(.sentinel_dir = output_dir, .sentinel_files = filenames)),
    add = TRUE), envir = env)
}


#' Check whether a file is a skip sentinel
#'
#' @param path  Path to file to check
#' @return TRUE if the file is a sentinel placeholder, FALSE otherwise
is_skip_sentinel <- function(path) {
    if (!file.exists(path)) return(FALSE)
    tryCatch({
        first_line <- readLines(path, n = 1, warn = FALSE)
        identical(trimws(first_line), "skipped=TRUE")
    }, error = function(e) FALSE)
}


#' Null-coalescing operator
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
