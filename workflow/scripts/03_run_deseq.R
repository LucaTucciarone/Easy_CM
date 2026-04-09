# =============================================================================
# 03_run_deseq.R
# EasyCM -- DESeq2 Marker Gene Detection
#
# PURPOSE:
#   Run DESeq2 on the combined interest/other count matrix to identify
#   genes differentially expressed in this cell type vs all others.
#   Design: ~ subtype + sample (cell type effect, controlling for donor).
#   Positive log2FoldChange = upregulated in the cell type of interest.
#
# USAGE:
#   Rscript workflow/scripts/03_run_deseq.R \
#       --config config/config.yaml --celltype Beta
#
# INPUTS:
#   results/{celltype}/intermediates/coldata.csv
#   results/{celltype}/intermediates/counts_filtered.csv
#
# OUTPUTS:
#   results/{celltype}/finals/results_deseq.tsv
#   results/{celltype}/plots/volcano.pdf
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(DESeq2)
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


# =============================================================================
# Main
# =============================================================================

main <- function(config_path, celltype) {

    cfg <- load_config(config_path)

    # --- Paths ---
    results_dir <- cfg$outputs$results_dir
    inter_dir   <- file.path(results_dir, celltype, "intermediates")
    finals_dir  <- file.path(results_dir, celltype, "finals")
    plots_dir   <- file.path(results_dir, celltype, "plots")
    dir.create(finals_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(plots_dir,  recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(cfg$outputs$logs_dir, paste0(celltype, ".log"))
    logger   <- make_logger(celltype = celltype, log_file = log_file,
                            level = cfg$logging$level %||% "INFO")

    # Register skip sentinels
    register_skip_sentinel(finals_dir, c("results_deseq.tsv"), env = environment())

    logger$info("start", sprintf("Running DESeq2 for celltype=%s", celltype))

    # --- Load inputs ---
    coldata_path <- file.path(inter_dir, "coldata.csv")
    counts_path  <- file.path(inter_dir, "counts_filtered.csv")

    if (is_skip_sentinel(coldata_path) || is_skip_sentinel(counts_path)) {
        logger$skip("upstream_skip", "Step 02 was skipped -- nothing to run")
        return(invisible(NULL))
    }

    coldata <- read.csv(coldata_path, row.names = 1, stringsAsFactors = FALSE)
    counts  <- read.csv(counts_path,  row.names = 1, check.names = FALSE)

    # Ensure factor levels
    coldata$subtype <- factor(coldata$subtype, levels = c("other", "interest"))
    coldata$sample  <- factor(coldata$sample)

    # Align column order
    counts <- counts[, rownames(coldata), drop = FALSE]

    logger$info("load_data", sprintf("n_genes=%d, n_columns=%d", nrow(counts), ncol(counts)))

    # --- Build DESeq2 dataset ---
    design_formula <- as.formula(cfg$deseq2$design %||% "~ subtype + sample")

    dds <- try_logged(
        DESeqDataSetFromMatrix(
            countData = round(as.matrix(counts)),
            colData   = coldata,
            design    = design_formula
        ),
        logger, "create_dds", "DESeqDataSet created"
    )
    if (is.null(dds)) return(invisible(NULL))

    # --- Run DESeq2 ---
    dds <- try_logged(
        {
            dds <- estimateSizeFactors(dds)
            dds <- estimateDispersions(dds)
            nbinomWaldTest(dds)
        },
        logger, "deseq2", "DESeq2 completed (sizeFactors + dispersions + Wald test)"
    )
    if (is.null(dds)) return(invisible(NULL))

    # --- Extract results ---
    # Contrast: interest vs other (positive LFC = upregulated in cell type)
    res <- try_logged(
        {
            res_raw <- results(dds, contrast = c("subtype", "interest", "other"))
            res_df  <- as.data.frame(res_raw)
            res_df  <- res_df[order(res_df$pvalue), ]
            res_df  <- cbind(gene = rownames(res_df), res_df)
            rownames(res_df) <- NULL
            res_df
        },
        logger, "extract_results", "Results extracted"
    )
    if (is.null(res)) return(invisible(NULL))

    # --- LFC shrinkage (apeglm) ---
    # Adds log2FoldChange_shrunk as a sibling column in the same tsv. The
    # regular log2FoldChange / stat / pvalue / padj columns stay as MLE
    # estimates. apeglm runs on the coef, not the contrast, so we need the
    # coefficient name DESeq2 assigned for "interest vs other". For the
    # design ~ subtype + sample with factor levels c("other","interest"),
    # that is "subtype_interest_vs_other".
    res$log2FoldChange_shrunk <- NA_real_

    shrunk_coef_name <- "subtype_interest_vs_other"
    available_coefs  <- resultsNames(dds)
    if (!shrunk_coef_name %in% available_coefs) {
        # Fall back to whichever coef name matches 'interest'
        fallback <- available_coefs[grepl("interest", available_coefs, fixed = TRUE)]
        if (length(fallback) == 1) shrunk_coef_name <- fallback
    }

    res_shrunk <- tryCatch(
        suppressMessages(DESeq2::lfcShrink(dds, coef = shrunk_coef_name,
                                           type = "apeglm", quiet = TRUE)),
        error = function(e) {
            logger$warn("lfc_shrinkage",
                sprintf("lfcShrink failed (coef=%s): %s",
                        shrunk_coef_name, e$message))
            NULL
        }
    )
    if (!is.null(res_shrunk)) {
        shrunk_df <- as.data.frame(res_shrunk)
        idx <- match(res$gene, rownames(shrunk_df))
        res$log2FoldChange_shrunk <- shrunk_df$log2FoldChange[idx]
        logger$info("lfc_shrinkage",
            sprintf("log2FoldChange_shrunk added (coef=%s)", shrunk_coef_name))
    }

    # --- Count DEGs ---
    fdr_threshold <- cfg$deseq2$fdr_threshold %||% 0.05
    lfc_threshold <- cfg$deseq2$lfc_threshold %||% 0
    n_degs <- sum(!is.na(res$padj) & res$padj < fdr_threshold &
                  abs(res$log2FoldChange) > lfc_threshold)
    n_up   <- sum(!is.na(res$padj) & res$padj < fdr_threshold &
                  res$log2FoldChange > lfc_threshold)
    n_down <- sum(!is.na(res$padj) & res$padj < fdr_threshold &
                  res$log2FoldChange < -lfc_threshold)

    logger$info("degs", sprintf("n_degs=%d (up=%d, down=%d) at FDR<%.2f",
                n_degs, n_up, n_down, fdr_threshold))

    # --- Write results ---
    write_result(res, file.path(finals_dir, "results_deseq.tsv"))

    # --- Save model info ---
    model_info <- data.frame(
        celltype     = celltype,
        formula      = deparse(design_formula),
        n_genes      = nrow(counts),
        n_columns    = ncol(counts),
        n_degs       = n_degs,
        n_up         = n_up,
        n_down       = n_down,
        fdr_threshold = fdr_threshold,
        stringsAsFactors = FALSE
    )
    write_result(model_info, file.path(inter_dir, "model_info.csv"), sep = ",")

    # --- Volcano plot ---
    tryCatch({
        plot_volcano(res, celltype, fdr_threshold, lfc_threshold,
                     file.path(plots_dir, "volcano.pdf"))
        logger$info("volcano_plot", "Saved")
    }, error = function(e) {
        logger$warn("volcano_plot", sprintf("Failed: %s", e$message))
    })

    logger$info("done", "DESeq2 step complete")
    return(invisible(NULL))
}


# =============================================================================
# Volcano plot
# =============================================================================

plot_volcano <- function(res, celltype, fdr_threshold, lfc_threshold, output_path) {

    res$significance <- "NS"
    res$significance[!is.na(res$padj) & res$padj < fdr_threshold &
                     res$log2FoldChange > lfc_threshold] <- "Up"
    res$significance[!is.na(res$padj) & res$padj < fdr_threshold &
                     res$log2FoldChange < -lfc_threshold] <- "Down"

    # Cap -log10(pvalue) for visualization
    res$neg_log10_pval <- -log10(pmax(res$pvalue, 1e-300))

    colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

    p <- ggplot(res, aes(x = log2FoldChange, y = neg_log10_pval, color = significance)) +
        geom_point(alpha = 0.5, size = 0.8) +
        scale_color_manual(values = colors) +
        geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "grey40") +
        labs(
            title = sprintf("%s -- Cell Type Markers", celltype),
            x = "log2 Fold Change (interest vs other)",
            y = "-log10(p-value)",
            color = "Significance"
        ) +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")

    ggsave(output_path, p, width = 8, height = 6)
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
