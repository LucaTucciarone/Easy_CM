# =============================================================================
# 05_aggregate_results.R
# EasyCM -- Aggregate Results Across Cell Types
#
# PURPOSE:
#   Collects per-celltype DESeq2 and fGSEA results into consolidated files.
#   Builds a celltype_summary.csv with status and DEG counts per cell type.
#   Generates a status overview heatmap.
#
# USAGE:
#   Rscript workflow/scripts/05_aggregate_results.R --config config/config.yaml
#
# OUTPUTS:
#   results/summary/all_markers_deseq.tsv
#   results/summary/all_fgsea.tsv
#   results/summary/celltype_summary.csv
#   results/summary/status_overview.pdf
#   results/summary/log_summary.tsv
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(ggplot2)
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
source(file.path(.script_dir, "utils", "validation_utils.R"))


# =============================================================================
# Helpers
# =============================================================================

#' Safely read a TSV file, returning NULL if missing or a skip sentinel
safe_read <- function(path) {
    if (!file.exists(path)) return(NULL)
    if (is_skip_sentinel(path)) return(NULL)
    tryCatch(fread(path) %>% as.data.frame(), error = function(e) NULL)
}


#' Determine the status for a cell type based on its outputs
determine_status <- function(celltype, results_dir, cfg) {

    inter_dir  <- file.path(results_dir, celltype, "intermediates")
    finals_dir <- file.path(results_dir, celltype, "finals")

    # Check coldata sentinel
    coldata_path <- file.path(inter_dir, "coldata.csv")
    if (!file.exists(coldata_path)) return("not_run")
    if (is_skip_sentinel(coldata_path)) {
        # Try to determine skip reason from log
        log_msg <- extract_skip_reason(celltype, cfg$outputs$logs_dir)
        return(log_msg$status)
    }

    # Check DESeq2 results
    deseq_path <- file.path(finals_dir, "results_deseq.tsv")
    if (!file.exists(deseq_path)) return("not_run")
    if (is_skip_sentinel(deseq_path)) return("error")

    # DESeq2 ran -- check for significant results
    res <- safe_read(deseq_path)
    if (is.null(res)) return("error")

    fdr <- cfg$deseq2$fdr_threshold %||% 0.05
    n_degs <- sum(!is.na(res$padj) & res$padj < fdr)

    if (n_degs == 0) return("success_no_significant")
    return("success")
}


#' Extract skip/error reason from log file
extract_skip_reason <- function(celltype, logs_dir) {
    log_file <- file.path(logs_dir, paste0(celltype, ".log"))
    if (!file.exists(log_file)) {
        return(list(status = "not_run", message = NA))
    }

    lines <- readLines(log_file, warn = FALSE)
    skip_lines  <- grep("\\[SKIP \\]", lines, value = TRUE)
    error_lines <- grep("\\[ERROR\\]", lines, value = TRUE)

    if (length(skip_lines) > 0) {
        last_skip <- tail(skip_lines, 1)
        reason <- sub(".*reason=([^ |]+).*", "\\1", last_skip)
        status <- switch(reason,
            "no_samples"                = "skipped_no_samples",
            "too_few_samples"           = "skipped_min_cells",
            "too_few_samples_after_dedup" = "skipped_min_cells",
            "too_few_genes"             = "skipped_preflight",
            "no_matrices"               = "skipped_no_samples",
            "celltype_not_found"        = "skipped_no_samples",
            "upstream_skip"             = "skipped",
            paste0("skipped_", reason)
        )
        return(list(status = status, message = last_skip))
    }

    if (length(error_lines) > 0) {
        return(list(status = "error", message = tail(error_lines, 1)))
    }

    return(list(status = "not_run", message = NA))
}


# =============================================================================
# Main
# =============================================================================

main <- function(config_path) {

    cfg <- load_config(config_path)

    results_dir <- cfg$outputs$results_dir
    summary_dir <- file.path(results_dir, "summary")
    dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

    cat("==================================================================\n")
    cat("           EasyCM -- Aggregate Results                             \n")
    cat("==================================================================\n")

    # --- Discover cell types from results directories ---
    all_dirs <- list.dirs(results_dir, full.names = FALSE, recursive = FALSE)
    celltypes <- all_dirs[all_dirs != "summary" & all_dirs != ""]

    if (length(celltypes) == 0) {
        stop("No cell type result directories found in: ", results_dir, call. = FALSE)
    }

    cat(sprintf("  Found %d cell type directories: %s\n",
                length(celltypes), paste(celltypes, collapse = ", ")))

    # --- Build summary rows ---
    summary_rows <- list()
    all_deseq    <- list()
    all_fgsea    <- list()

    for (ct in celltypes) {

        finals_dir <- file.path(results_dir, ct, "finals")
        inter_dir  <- file.path(results_dir, ct, "intermediates")
        status     <- determine_status(ct, results_dir, cfg)
        skip_info  <- extract_skip_reason(ct, cfg$outputs$logs_dir)

        # Read model info if available
        model_info_path <- file.path(inter_dir, "model_info.csv")
        model_info <- safe_read(model_info_path)

        # Read DESeq2 results
        deseq_path <- file.path(finals_dir, "results_deseq.tsv")
        deseq_res  <- safe_read(deseq_path)

        n_degs <- NA
        formula_used <- NA
        if (!is.null(deseq_res)) {
            fdr <- cfg$deseq2$fdr_threshold %||% 0.05
            n_degs <- sum(!is.na(deseq_res$padj) & deseq_res$padj < fdr)
            deseq_res$celltype <- ct
            all_deseq[[ct]] <- deseq_res
        }
        if (!is.null(model_info)) {
            formula_used <- model_info$formula[1]
        }

        # Read fGSEA results
        fgsea_path <- file.path(finals_dir, "fgsea_all.tsv")
        fgsea_res  <- safe_read(fgsea_path)

        n_pathways <- NA
        if (!is.null(fgsea_res)) {
            fgsea_fdr  <- cfg$fgsea$fdr_threshold %||% 0.10
            n_pathways <- sum(!is.na(fgsea_res$padj) & fgsea_res$padj < fgsea_fdr)
            fgsea_res$celltype <- ct
            all_fgsea[[ct]] <- fgsea_res
        }

        # Read coldata for sample counts
        coldata_path <- file.path(inter_dir, "coldata.csv")
        coldata <- safe_read(coldata_path)

        n_interest <- NA
        n_other    <- NA
        n_total    <- NA
        if (!is.null(coldata) && "subtype" %in% colnames(coldata)) {
            n_interest <- sum(coldata$subtype == "interest")
            n_other    <- sum(coldata$subtype == "other")
            n_total    <- nrow(coldata)
        }

        summary_rows[[ct]] <- make_summary_row(
            celltype        = ct,
            status          = status,
            error_message   = if (status %in% c("error", "not_run") ||
                                  grepl("^skipped", status)) skip_info$message else NA,
            n_interest      = n_interest,
            n_other         = n_other,
            n_samples_total = n_total,
            formula_used    = formula_used,
            n_degs          = n_degs,
            n_pathways      = n_pathways
        )

        cat(sprintf("  [%s] %s: %s (DEGs=%s, pathways=%s)\n",
                    status, ct,
                    ifelse(is.na(skip_info$message), "", skip_info$message),
                    ifelse(is.na(n_degs), "-", as.character(n_degs)),
                    ifelse(is.na(n_pathways), "-", as.character(n_pathways))))
    }

    # --- Write consolidated files ---
    summary_df <- do.call(rbind, summary_rows)
    write_result(summary_df, file.path(summary_dir, "celltype_summary.csv"), sep = ",")

    if (length(all_deseq) > 0) {
        combined_deseq <- do.call(rbind, all_deseq)
        write_result(combined_deseq, file.path(summary_dir, "all_markers_deseq.tsv"))
    }

    if (length(all_fgsea) > 0) {
        combined_fgsea <- do.call(rbind, all_fgsea)
        fwrite(combined_fgsea, file.path(summary_dir, "all_fgsea.tsv"),
               sep = "\t", sep2 = c("", " ", ""), quote = FALSE)
        message(sprintf("  [write] %s", file.path(summary_dir, "all_fgsea.tsv")))
    }

    # --- Log summary ---
    summarize_logs(cfg$outputs$logs_dir, file.path(summary_dir, "log_summary.tsv"))

    # --- Status heatmap ---
    tryCatch({
        plot_status_heatmap(summary_df, file.path(summary_dir, "status_overview.pdf"))
        cat("  [write] status_overview.pdf\n")
    }, error = function(e) {
        message(sprintf("  [warn] Could not generate status heatmap: %s", e$message))
    })

    # --- Final report ---
    n_success <- sum(grepl("^success", summary_df$status))
    n_skipped <- sum(grepl("^skipped", summary_df$status))
    n_error   <- sum(summary_df$status == "error")
    n_not_run <- sum(summary_df$status == "not_run")

    cat(sprintf("\n  Summary: %d success, %d skipped, %d error, %d not_run\n",
                n_success, n_skipped, n_error, n_not_run))
    cat("  Aggregation complete.\n\n")

    invisible(summary_df)
}


# =============================================================================
# Status heatmap
# =============================================================================

plot_status_heatmap <- function(summary_df, output_path) {

    status_colors <- c(
        "success"               = "#4CAF50",
        "success_no_significant" = "#81C784",
        "skipped_no_samples"    = "#FFB74D",
        "skipped_min_cells"     = "#FFA726",
        "skipped_preflight"     = "#FF9800",
        "skipped"               = "#FFE0B2",
        "error"                 = "#E53935",
        "not_run"               = "#BDBDBD"
    )

    summary_df$status <- factor(summary_df$status,
        levels = names(status_colors))

    # Add DEG count labels
    summary_df$label <- ifelse(is.na(summary_df$n_degs), "",
                               as.character(summary_df$n_degs))

    p <- ggplot(summary_df, aes(x = celltype, y = 1, fill = status)) +
        geom_tile(color = "white", linewidth = 1) +
        geom_text(aes(label = label), size = 3, color = "black") +
        scale_fill_manual(values = status_colors, drop = FALSE) +
        labs(
            title = "EasyCM -- Cell Type Marker Status",
            x = "Cell Type",
            y = NULL,
            fill = "Status"
        ) +
        theme_minimal(base_size = 12) +
        theme(
            axis.text.x  = element_text(angle = 45, hjust = 1),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid   = element_blank(),
            legend.position = "bottom"
        )

    ggsave(output_path, p, width = max(8, length(unique(summary_df$celltype)) * 0.8),
           height = 4)
}


# =============================================================================
# CLI entry point
# =============================================================================

if (!interactive()) {

    opt_list <- list(
        make_option("--config", type = "character",
                    default = "config/config.yaml",
                    help = "Path to config.yaml")
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
