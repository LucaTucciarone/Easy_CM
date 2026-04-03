# =============================================================================
# logging_utils.R
# Structured logging for EasyCM
#
# Every log line is machine-readable AND human-readable:
#   [2026-04-02 14:32:01] [INFO]  celltype=Beta | step=filter_samples | n_samples=47->31
#   [2026-04-02 14:32:03] [ERROR] celltype=Acinar | step=deseq | msg=every gene contains at least one zero
#
# Two log files per cell type:
#   {celltype}.log        -- full structured log (all levels)
#   {celltype}.errors.log -- ONLY created if WARN/ERROR/SKIP lines exist
#
# Usage:
#   logger <- make_logger(celltype = "Beta",
#                         log_file = "logs/Beta.log", level = "INFO")
#   logger$info("filter_samples", "n_samples=47->31")
#   logger$warn("donor_dedup", "4 samples removed")
#   logger$error("deseq", "every gene contains at least one zero")
#   logger$skip("too_few_samples", "only 2 samples after filtering")
# =============================================================================


# Log level hierarchy -- messages below the configured level are suppressed
LOG_LEVELS <- c(DEBUG = 0, INFO = 1, WARN = 2, ERROR = 3, SKIP = 3)

#' Append a formatted error to a running error list
#'
#' Used by validation scripts that collect ALL errors before stopping.
#'
#' @param errors  Character vector of accumulated error strings
#' @param msg     Human-readable description of the error
#' @return Updated character vector with the new error appended
add_error <- function(errors, msg) {
    c(errors, paste0("  [ERROR] ", msg))
}


#' Append a formatted warning to a running warning list
#'
#' @param warnings  Character vector of accumulated warning strings
#' @param msg       Human-readable description of the warning
#' @return Updated character vector with the new warning appended
add_warning <- function(warnings, msg) {
    c(warnings, paste0("  [WARN]  ", msg))
}


#' Create a logger instance bound to a specific cell type
#'
#' Returns a named list of logging functions (info, warn, error, skip, debug).
#' Each function writes a structured line to the full log file and stderr.
#' WARN/ERROR/SKIP lines are also written to a separate .errors.log file.
#'
#' @param celltype   Cell type identifier string (e.g. "Beta")
#' @param log_file   Path to the full log file. Created if missing.
#' @param level      Minimum log level: "DEBUG" | "INFO" | "WARN" | "ERROR"
#' @return Named list with functions: $debug, $info, $warn, $error, $skip
make_logger <- function(celltype, log_file, level = "INFO") {

    min_level <- LOG_LEVELS[[toupper(level)]]
    if (is.null(min_level)) {
        stop(sprintf("Invalid log level: '%s'. Choose from: DEBUG, INFO, WARN, ERROR", level))
    }

    error_log_file <- sub("\\.log$", ".errors.log", log_file)
    dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)

    write_log <- function(log_level, step, message) {

        if (LOG_LEVELS[[log_level]] < min_level) return(invisible(NULL))

        timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

        line <- sprintf(
            "[%s] [%-5s] celltype=%s | step=%s | %s",
            timestamp,
            log_level,
            celltype,
            step,
            message
        )

        cat(line, "\n", file = log_file, append = TRUE)

        if (log_level %in% c("WARN", "ERROR", "SKIP")) {
            cat(line, "\n", file = error_log_file, append = TRUE)
        }

        if (log_level %in% c("WARN", "ERROR", "SKIP")) {
            message(line)
        }

        invisible(NULL)
    }

    list(
        debug = function(step, msg) write_log("DEBUG", step, msg),
        info  = function(step, msg) write_log("INFO",  step, msg),
        warn  = function(step, msg) write_log("WARN",  step, msg),
        error = function(step, msg) write_log("ERROR", step, msg),
        skip  = function(reason, msg) write_log("SKIP",
                    step    = "SKIPPED",
                    message = sprintf("reason=%s | %s", reason, msg))
    )
}


#' Summarize all log files into a single error/skip report
#'
#' @param logs_dir     Root logs directory
#' @param output_file  Path to write the summary report (TSV)
#' @return data.frame of flagged log lines (invisibly)
summarize_logs <- function(logs_dir, output_file) {

    log_files <- list.files(logs_dir, pattern = "\\.log$", recursive = TRUE, full.names = TRUE)

    if (length(log_files) == 0) {
        message("No log files found in: ", logs_dir)
        return(invisible(data.frame()))
    }

    flagged_lines <- lapply(log_files, function(f) {
        lines <- readLines(f, warn = FALSE)
        flagged <- lines[grepl("\\[(WARN |ERROR|SKIP )\\]", lines)]
        if (length(flagged) == 0) return(NULL)
        data.frame(log_file = f, line = flagged, stringsAsFactors = FALSE)
    })

    flagged_df <- do.call(rbind, Filter(Negate(is.null), flagged_lines))

    if (is.null(flagged_df) || nrow(flagged_df) == 0) {
        message("No warnings, errors or skips found across all logs.")
        return(invisible(data.frame()))
    }

    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    write.table(flagged_df, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    message(sprintf("[log summary] %d flagged lines written to: %s", nrow(flagged_df), output_file))

    return(invisible(flagged_df))
}


#' Wrap a code block in tryCatch with automatic structured logging
#'
#' On success: logs INFO and returns result.
#' On failure: logs ERROR and returns NULL (caller handles the skip).
#'
#' @param expr     Expression to evaluate
#' @param logger   Logger instance from make_logger()
#' @param step     Step name string
#' @param success  Message to log on success (default: "completed")
#' @return Result of expr on success, NULL on failure
try_logged <- function(expr, logger, step, success = "completed") {

    result <- tryCatch(
        {
            out <- expr
            logger$info(step, success)
            out
        },
        error = function(e) {
            msg <- conditionMessage(e)
            if (is.null(msg) || !nzchar(trimws(msg))) {
                msg <- paste0("unknown error (class: ", paste(class(e), collapse = "/"), ")")
            }
            logger$error(step, sprintf("msg=%s", msg))
            NULL
        }
    )

    return(result)
}
