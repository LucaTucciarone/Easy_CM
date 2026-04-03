# =============================================================================
# 04_run_fgsea.R
# EasyCM -- fGSEA Pathway Enrichment on Marker Genes
#
# PURPOSE:
#   Run fast gene set enrichment analysis (fGSEA) on DESeq2 marker results.
#   Genes are ranked by (-log10(pvalue)) * log2FoldChange, after removing
#   ribosomal and mitochondrial genes.
#
# USAGE:
#   Rscript workflow/scripts/04_run_fgsea.R \
#       --config config/config.yaml --celltype Beta
#
# INPUTS:
#   results/{celltype}/finals/results_deseq.tsv
#
# OUTPUTS:
#   results/{celltype}/finals/fgsea_all.tsv
#   results/{celltype}/finals/fgsea_significant.tsv
#   results/{celltype}/plots/fgsea_pathways.pdf
# =============================================================================

suppressMessages({
    library(yaml)
    library(data.table)
    library(dplyr)
    library(tibble)
    library(fgsea)
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

    # Check if fGSEA is enabled
    if (!isTRUE(cfg$steps$fgsea)) {
        message("fGSEA step disabled in config -- skipping")
        return(invisible(NULL))
    }

    # --- Paths ---
    results_dir <- cfg$outputs$results_dir
    finals_dir  <- file.path(results_dir, celltype, "finals")
    plots_dir   <- file.path(results_dir, celltype, "plots")
    dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

    log_file <- file.path(cfg$outputs$logs_dir, paste0(celltype, ".log"))
    logger   <- make_logger(celltype = celltype, log_file = log_file,
                            level = cfg$logging$level %||% "INFO")

    # Register skip sentinels
    register_skip_sentinel(finals_dir, c("fgsea_all.tsv"), env = environment())

    logger$info("start", sprintf("Running fGSEA for celltype=%s", celltype))

    # --- Load DESeq2 results ---
    deseq_path <- file.path(finals_dir, "results_deseq.tsv")
    if (!file.exists(deseq_path) || is_skip_sentinel(deseq_path)) {
        logger$skip("upstream_skip", "DESeq2 results not available")
        return(invisible(NULL))
    }

    res <- fread(deseq_path) %>% as.data.frame()
    logger$info("load_results", sprintf("n_genes=%d", nrow(res)))

    # --- Load gene exclusion lists ---
    exclude_genes <- character(0)
    if (!is.null(cfg$fgsea$exclude_gene_lists)) {
        for (gene_list_path in cfg$fgsea$exclude_gene_lists) {
            if (file.exists(gene_list_path)) {
                gl <- fread(gene_list_path, fill = TRUE, header = TRUE)
                # Look for "Approved symbol" column (HGNC format) or first column
                if ("Approved symbol" %in% colnames(gl)) {
                    genes <- gl[["Approved symbol"]]
                } else {
                    genes <- gl[[1]]
                }
                genes <- genes[!is.na(genes) & genes != "" & genes != "Approved symbol"]
                exclude_genes <- c(exclude_genes, genes)
            }
        }
    }

    # --- Filter genes ---
    res_filtered <- res[!res$gene %in% exclude_genes, ]

    # Apply regex exclusion patterns
    if (!is.null(cfg$fgsea$exclude_gene_patterns)) {
        for (pattern in cfg$fgsea$exclude_gene_patterns) {
            res_filtered <- res_filtered[!grepl(pattern, res_filtered$gene), ]
        }
    }

    logger$info("gene_exclusion", sprintf("n_genes=%d->%d after excluding ribosomal/mito",
                nrow(res), nrow(res_filtered)))

    # --- Rank genes ---
    ranking_stat <- cfg$fgsea$ranking_stat %||% "stat"

    if (ranking_stat == "stat" && "stat" %in% colnames(res_filtered)) {
        # Use DESeq2 Wald statistic directly
        res_filtered$rank_value <- res_filtered$stat
    } else {
        # Use (-log10(pvalue)) * log2FoldChange (OG pipeline method)
        res_filtered$pvalue <- ifelse(res_filtered$pvalue == 0, 1e-306, res_filtered$pvalue)
        res_filtered$rank_value <- (-log10(res_filtered$pvalue)) * res_filtered$log2FoldChange
    }

    # Remove NA ranks
    res_filtered <- res_filtered[!is.na(res_filtered$rank_value), ]

    # Create named vector of ranks
    ranks <- setNames(res_filtered$rank_value, res_filtered$gene)
    ranks <- sort(ranks, decreasing = TRUE)

    logger$info("ranking", sprintf("n_ranked_genes=%d, stat=%s", length(ranks), ranking_stat))

    # --- Load gene sets ---
    gene_sets_file <- cfg$fgsea$gene_sets_file
    if (is.null(gene_sets_file) || !file.exists(gene_sets_file)) {
        logger$error("load_pathways", sprintf("Gene sets file not found: %s", gene_sets_file))
        return(invisible(NULL))
    }

    pathways <- gmtPathways(gene_sets_file)
    logger$info("load_pathways", sprintf("n_pathways=%d", length(pathways)))

    # --- Run fGSEA ---
    min_size <- cfg$fgsea$min_gene_set_size %||% 10
    max_size <- cfg$fgsea$max_gene_set_size %||% 500

    fgsea_res <- try_logged(
        fgseaMultilevel(
            pathways = pathways,
            stats    = ranks,
            eps      = 0.0,
            minSize  = min_size,
            maxSize  = max_size
        ),
        logger, "fgsea", "fGSEA completed"
    )
    if (is.null(fgsea_res)) return(invisible(NULL))

    fgsea_res <- fgsea_res[order(fgsea_res$pval), ]

    logger$info("fgsea_results", sprintf("n_total=%d pathways tested", nrow(fgsea_res)))

    # --- Filter significant ---
    fdr_threshold <- cfg$fgsea$fdr_threshold %||% 0.10
    fgsea_signif  <- fgsea_res[fgsea_res$padj < fdr_threshold, ]

    logger$info("fgsea_significant", sprintf("n_significant=%d at FDR<%.2f",
                nrow(fgsea_signif), fdr_threshold))

    # --- Write results ---
    fwrite(fgsea_res,    file.path(finals_dir, "fgsea_all.tsv"),
           sep = "\t", sep2 = c("", " ", ""), quote = FALSE)
    message(sprintf("  [write] %s", file.path(finals_dir, "fgsea_all.tsv")))

    fwrite(fgsea_signif, file.path(finals_dir, "fgsea_significant.tsv"),
           sep = "\t", sep2 = c("", " ", ""), quote = FALSE)
    message(sprintf("  [write] %s", file.path(finals_dir, "fgsea_significant.tsv")))

    # --- Pathway bar chart ---
    n_top <- cfg$fgsea$n_top_pathways_plot %||% 20
    tryCatch({
        plot_fgsea_bars(fgsea_res, celltype, n_top,
                        file.path(plots_dir, "fgsea_pathways.pdf"))
        logger$info("fgsea_plot", "Saved")
    }, error = function(e) {
        logger$warn("fgsea_plot", sprintf("Failed: %s", e$message))
    })

    logger$info("done", "fGSEA step complete")
    return(invisible(NULL))
}


# =============================================================================
# fGSEA pathway bar chart
# =============================================================================

plot_fgsea_bars <- function(fgsea_res, celltype, n_top, output_path) {

    # Take top N pathways by p-value
    top_pathways <- head(fgsea_res[order(fgsea_res$pval), ], n_top)

    if (nrow(top_pathways) == 0) return(invisible(NULL))

    top_pathways$pathway_short <- gsub("_", " ", top_pathways$pathway)
    # Truncate long names
    top_pathways$pathway_short <- ifelse(
        nchar(top_pathways$pathway_short) > 60,
        paste0(substr(top_pathways$pathway_short, 1, 57), "..."),
        top_pathways$pathway_short
    )

    top_pathways$direction <- ifelse(top_pathways$NES > 0, "Enriched", "Depleted")

    p <- ggplot(top_pathways, aes(x = reorder(pathway_short, NES), y = NES, fill = direction)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("Enriched" = "#E41A1C", "Depleted" = "#377EB8")) +
        labs(
            title = sprintf("%s -- Top %d Pathways (fGSEA)", celltype, n_top),
            x = NULL,
            y = "Normalized Enrichment Score (NES)",
            fill = NULL
        ) +
        theme_minimal(base_size = 10) +
        theme(legend.position = "bottom")

    ggsave(output_path, p, width = 10, height = max(6, n_top * 0.3))
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
