# =============================================================================
# 06_pipeline_summary.R
# EasyCM -- Pipeline-Level Summary
#
# PURPOSE:
#   Generates the final pipeline summary: a CSV with one row per cell type
#   and a summary PDF showing the marker gene landscape across all cell types.
#
# USAGE:
#   Rscript workflow/scripts/06_pipeline_summary.R --config config/config.yaml
#
# INPUTS:
#   results/summary/celltype_summary.csv
#   results/summary/all_markers_deseq.tsv
#
# OUTPUTS:
#   results/pipeline_summary.csv
#   results/pipeline_summary.pdf
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


# =============================================================================
# Main
# =============================================================================

main <- function(config_path) {

    cfg <- load_config(config_path)
    results_dir <- cfg$outputs$results_dir
    summary_dir <- file.path(results_dir, "summary")

    cat("==================================================================\n")
    cat("           EasyCM -- Pipeline Summary                              \n")
    cat("==================================================================\n")

    # --- Read celltype summary ---
    summary_path <- file.path(summary_dir, "celltype_summary.csv")
    if (!file.exists(summary_path)) {
        stop("celltype_summary.csv not found. Run step 05 first.", call. = FALSE)
    }

    summary_df <- read.csv(summary_path, stringsAsFactors = FALSE)

    # Copy to pipeline root
    write_result(summary_df, file.path(results_dir, "pipeline_summary.csv"), sep = ",")

    cat(sprintf("  %d cell types processed\n", nrow(summary_df)))
    cat(sprintf("  Success: %d\n", sum(grepl("^success", summary_df$status))))
    cat(sprintf("  Skipped: %d\n", sum(grepl("^skipped", summary_df$status))))
    cat(sprintf("  Error:   %d\n", sum(summary_df$status == "error")))

    # --- Summary bar chart: DEGs per cell type ---
    tryCatch({
        successful <- summary_df[grepl("^success", summary_df$status) & !is.na(summary_df$n_degs), ]

        if (nrow(successful) > 0) {
            p <- ggplot(successful, aes(x = reorder(celltype, -n_degs), y = n_degs)) +
                geom_col(fill = "#4CAF50", alpha = 0.8) +
                geom_text(aes(label = n_degs), vjust = -0.3, size = 3) +
                labs(
                    title = "EasyCM -- Marker Genes per Cell Type",
                    subtitle = sprintf("FDR < %s", cfg$deseq2$fdr_threshold %||% 0.05),
                    x = "Cell Type",
                    y = "Number of Marker Genes (DEGs)"
                ) +
                theme_minimal(base_size = 12) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))

            output_pdf <- file.path(results_dir, "pipeline_summary.pdf")
            ggsave(output_pdf, p,
                   width = max(8, nrow(successful) * 0.8), height = 6)
            cat(sprintf("  [write] %s\n", output_pdf))
        } else {
            cat("  No successful cell types to plot.\n")
            # Write empty PDF sentinel
            pdf(file.path(results_dir, "pipeline_summary.pdf"), width = 1, height = 1)
            plot.new()
            text(0.5, 0.5, "No results")
            dev.off()
        }
    }, error = function(e) {
        message(sprintf("  [warn] Could not generate summary plot: %s", e$message))
        # Write empty PDF so Snakemake output contract is satisfied
        pdf(file.path(results_dir, "pipeline_summary.pdf"), width = 1, height = 1)
        plot.new()
        text(0.5, 0.5, "Plot failed")
        dev.off()
    })

    cat("  Pipeline summary complete.\n\n")
    invisible(summary_df)
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
