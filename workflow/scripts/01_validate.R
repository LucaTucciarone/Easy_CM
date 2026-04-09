# =============================================================================
# 01_validate.R
# EasyCM -- Input Validation
#
# PURPOSE:
#   Validates all inputs before any compute runs. Checks config.yaml,
#   verifies input files exist, validates count matrix structure,
#   and discovers cell types to analyze.
#
# USAGE:
#   Rscript workflow/scripts/01_validate.R --config config/config.yaml
#
# OUTPUT:
#   Prints validation report to stdout.
#   Writes logs/validation.log
#   On success: exits 0
#   On failure: exits 1 with consolidated error list.
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(optparse)
})

# --- Resolve utils directory relative to this script ---
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
source(file.path(.script_dir, "utils", "validation_utils.R"))


# =============================================================================
# Main
# =============================================================================

main <- function(config_path) {

    all_errors   <- character(0)
    all_warnings <- character(0)

    cat("==================================================================\n")
    cat("           EasyCM -- Input Validation                              \n")
    cat("==================================================================\n")
    cat(sprintf("Config: %s\n", config_path))
    cat(sprintf("Time:   %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

    # --- Load config ---
    if (!file.exists(config_path)) {
        stop(sprintf("Cannot find config file: '%s'", config_path))
    }
    cfg <- suppressWarnings(yaml::read_yaml(config_path))
    cat("  [OK] config.yaml loaded\n")

    # --- Validate config structure ---
    config_result <- validate_config(cfg)
    all_errors    <- c(all_errors,   config_result$errors)
    all_warnings  <- c(all_warnings, config_result$warnings)

    if (length(all_errors) > 0) {
        cat("\n-- Validation failed -- fix config errors before continuing -----\n")
        cat(paste(all_errors, collapse = "\n"), "\n")
        stop("Validation failed. See errors above.", call. = FALSE)
    }

    # --- Load and validate metadata ---
    cat("\n-- Loading input files ---------------------------------------------\n")

    sample_meta <- tryCatch(
        fread(cfg$inputs$sample_metadata) %>% as.data.frame(),
        error = function(e) {
            stop(sprintf("Could not read sample_metadata: %s", e$message), call. = FALSE)
        }
    )
    cat(sprintf("  [OK] sample_metadata: %d rows x %d columns\n",
                nrow(sample_meta), ncol(sample_meta)))

    # Check required metadata columns
    required_meta_cols <- c(cfg$filtering$donor_col)
    missing_meta <- setdiff(required_meta_cols, colnames(sample_meta))
    if (length(missing_meta) > 0) {
        all_errors <- c(all_errors, sprintf(
            "  [ERROR] Required columns not in sample_metadata: %s",
            paste(missing_meta, collapse = ", ")))
    }

    # --- Discover and validate cell types ---
    cat("\n-- Discovering cell types ------------------------------------------\n")

    # Discover from EasyPseudobulk counts/easycm/ layout
    discovered <- discover_celltypes(cfg$inputs$counts_dir)
    cat(sprintf("  Found %d cell type count matrices: %s\n",
                length(discovered), paste(discovered, collapse = ", ")))

    if (length(cfg$inputs$celltypes) > 0) {
        celltypes <- cfg$inputs$celltypes
        missing_ct <- setdiff(celltypes, discovered)
        if (length(missing_ct) > 0) {
            all_errors <- c(all_errors, sprintf(
                "  [ERROR] Requested cell types not found in counts_dir: %s",
                paste(missing_ct, collapse = ", ")))
        }
    } else {
        celltypes <- discovered
    }

    if (length(celltypes) == 0) {
        all_errors <- c(all_errors,
            "  [ERROR] No cell types found. Check counts_dir or celltypes config.")
    }

    # --- Final report ---
    cat("\n-- Validation Report -----------------------------------------------\n")

    if (length(all_warnings) > 0) {
        cat(sprintf("\n%d warning(s):\n", length(all_warnings)))
        cat(paste(all_warnings, collapse = "\n"), "\n")
    }

    if (length(all_errors) > 0) {
        cat(sprintf("\n%d error(s) -- pipeline cannot run:\n", length(all_errors)))
        cat(paste(all_errors, collapse = "\n"), "\n")
        stop("Validation failed.", call. = FALSE)
    }

    # --- Create output directory stubs ---
    dir.create(cfg$outputs$results_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(cfg$outputs$logs_dir,    recursive = TRUE, showWarnings = FALSE)

    # Write validation log
    validation_log <- file.path(cfg$outputs$logs_dir, "validation.log")
    writeLines(c(
        sprintf("EasyCM validation: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        sprintf("Config:     %s", config_path),
        sprintf("Cell types: %s", paste(celltypes, collapse = ", ")),
        sprintf("Warnings:   %d", length(all_warnings)),
        if (length(all_warnings) > 0) paste(all_warnings, collapse = "\n") else ""
    ), validation_log)

    cat(sprintf("\n  %d cell type(s) ready to analyze: %s\n",
                length(celltypes), paste(celltypes, collapse = ", ")))
    cat(sprintf("  Validation log written: %s\n", validation_log))
    cat("  Ready to run.\n\n")

    invisible(list(cfg = cfg, celltypes = celltypes, sample_meta = sample_meta))
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive()) {

    opt_list <- list(
        make_option("--config", type = "character",
                    default = "config/config.yaml",
                    help = "Path to config.yaml [default: config/config.yaml]")
    )

    opts <- parse_args(OptionParser(option_list = opt_list))
    tryCatch({
        main(opts$config)
        quit(status = 0)
    }, error = function(e) {
        message("ERROR: ", conditionMessage(e))
        quit(status = 1)
    })
}
