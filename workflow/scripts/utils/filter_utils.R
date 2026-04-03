# =============================================================================
# filter_utils.R
# Gene and sample filtering functions for EasyCM
#
# All functions are PURE: inputs in, outputs out, no side effects.
# =============================================================================


# -----------------------------------------------------------------------------
# Gene-level filters
# -----------------------------------------------------------------------------

#' Proportion-based gene filter
#'
#' Keeps genes with >= min_reads counts in >= min_prop fraction of samples.
#'
#' @param row       Numeric vector of counts for one gene across samples
#' @param min_reads Minimum count threshold (default 5)
#' @param min_prop  Minimum proportion of samples that must meet threshold (default 0.25)
#' @return Logical scalar -- TRUE means keep this gene
proportion_filter <- function(row, min_reads = 5, min_prop = 0.25) {
    mean(row >= min_reads) >= min_prop
}


#' Apply proportion filter to a count matrix
#'
#' Filters genes (rows) that pass the proportion threshold across samples.
#' Uses the "interest" group samples for filtering (consistent with OG pipeline).
#'
#' @param count_matrix  Genes x samples matrix (numeric, rownames = gene names)
#' @param min_reads     Minimum count threshold
#' @param min_prop      Minimum proportion of samples passing threshold
#' @return Filtered count matrix (subset of rows)
filter_genes <- function(count_matrix, min_reads = 5, min_prop = 0.25) {
    keep <- apply(count_matrix, 1, proportion_filter,
                  min_reads = min_reads, min_prop = min_prop)
    count_matrix[keep, , drop = FALSE]
}


# -----------------------------------------------------------------------------
# Sample-level filters
# -----------------------------------------------------------------------------

#' Apply an arbitrary list of subset filters to a metadata table
#'
#' @param meta    data.frame of sample metadata
#' @param filters Named list: column names -> allowed value(s)
#'                Empty list = no filtering.
#' @return Filtered data.frame
apply_subset_filters <- function(meta, filters) {

    if (length(filters) == 0) return(meta)

    for (col_name in names(filters)) {

        allowed_values <- filters[[col_name]]

        if (!col_name %in% colnames(meta)) {
            warning(sprintf(
                "Filter column '%s' not found in metadata -- skipping this filter.",
                col_name
            ))
            next
        }

        before <- nrow(meta)
        meta   <- meta[meta[[col_name]] %in% allowed_values, , drop = FALSE]
        after  <- nrow(meta)

        message(sprintf(
            "  Filter [%s = %s]: %d -> %d samples",
            col_name, paste(allowed_values, collapse = " | "), before, after
        ))
    }

    return(meta)
}


#' Deduplicate donors: keep one sample per donor
#'
#' @param meta       data.frame of sample metadata
#' @param donor_col  Column name holding donor/subject IDs
#' @return Deduplicated data.frame
dedup_donors <- function(meta, donor_col) {

    if (!donor_col %in% colnames(meta)) {
        stop(sprintf(
            "donor_col '%s' not found in metadata. Check filtering.donor_col in config.yaml.",
            donor_col
        ))
    }

    donor_counts <- table(meta[[donor_col]])

    # Donors with exactly 1 sample -- keep as-is
    single <- meta[meta[[donor_col]] %in% names(donor_counts[donor_counts == 1]), ]

    # Donors with >1 sample -- keep first deterministically
    multi_donors <- names(donor_counts[donor_counts > 1])
    set.seed(123)
    multi_rows <- lapply(multi_donors, function(d) {
        rows <- meta[meta[[donor_col]] == d, ]
        rows[sample(nrow(rows), 1), ]
    })

    result <- do.call(rbind, c(list(single), multi_rows))

    message(sprintf(
        "  Donor dedup: %d samples -> %d donors retained",
        nrow(meta), nrow(result)
    ))

    return(result)
}


#' Filter samples by minimum cell count for a given cell type
#'
#' Removes samples where the cell type of interest has fewer than
#' min_cells cells. Requires a cell-count column in metadata.
#'
#' @param meta          data.frame of sample metadata
#' @param celltype      Name of the cell type (used to find count column)
#' @param min_cells     Minimum cell count threshold
#' @param count_col     Column name with cell counts for this type (if NULL, tries celltype name)
#' @return Filtered data.frame
filter_by_cell_count <- function(meta, celltype, min_cells, count_col = NULL) {

    if (is.null(count_col)) count_col <- celltype

    if (!count_col %in% colnames(meta)) {
        message(sprintf("  Cell count column '%s' not found -- skipping min-cell filter", count_col))
        return(meta)
    }

    before <- nrow(meta)
    meta <- meta[!is.na(meta[[count_col]]) & meta[[count_col]] >= min_cells, , drop = FALSE]
    after <- nrow(meta)

    message(sprintf(
        "  Min cell filter [%s >= %d]: %d -> %d samples",
        count_col, min_cells, before, after
    ))

    return(meta)
}
