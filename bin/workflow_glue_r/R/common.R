workflow_glue_r_is_r_formula_name <- function(name) {
    is.character(name) &&
        length(name) == 1 &&
        grepl("^[A-Za-z][A-Za-z0-9_.]*$", name)
}

workflow_glue_r_validate_r_formula_names <- function(names, label = "Column") {
    invalid <- names[!vapply(names, workflow_glue_r_is_r_formula_name, logical(1))]
    if (length(invalid) > 0) {
        stop(
            sprintf(
                "%s names must be safe for R formulas. Invalid names: %s. Names must start with a letter and contain only letters, numbers, underscores, and dots.",
                label,
                paste(invalid, collapse = ", ")
            ),
            call. = FALSE
        )
    }
    invisible(names)
}

workflow_glue_r_normalise_tsv_value <- function(value) {
    if (length(value) == 0 || all(is.na(value))) {
        return(NA_character_)
    }
    if (is.list(value)) {
        value <- unlist(value, recursive = TRUE, use.names = FALSE)
    }
    if (length(value) == 0 || all(is.na(value))) {
        return(NA_character_)
    }
    paste(as.character(value), collapse = ";")
}

workflow_glue_r_normalise_tsv_df <- function(df) {
    as.data.frame(
        lapply(df, function(column) {
            if (is.list(column)) {
                vapply(column, workflow_glue_r_normalise_tsv_value, character(1))
            } else {
                column
            }
        }),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

bambu_normalise_tsv_df <- workflow_glue_r_normalise_tsv_df

workflow_glue_r_empty_tsv <- function(columns) {
    out <- as.data.frame(matrix(nrow = 0, ncol = length(columns)))
    names(out) <- columns
    out
}

#' Extract a deduplicated id-to-name mapping from a GFF annotation file.
#'
#' Uses annotation metadata object from rtracklayer to generate a two-column
#' `data.frame` mapping feature IDs to display names.  
#' When a feature ID maps to multiple names only the first observed name is kept.
#'
#' @param annotation_path data.frame with annotation metadata
#' @param id_column Name of the GFF attribute to use as the identifier
#'   (e.g. `gene_id` or `transcript_id`).
#' @param name_column Name of the GFF attribute to use as the display name
#'   (e.g. `gene_name` or `transcript_name`).
#' @param output_id_column Column name for the identifier in the returned
#'   `data.frame` (e.g. `GENEID` or `TXNAME`).
#'
#' @return A `data.frame` with columns `output_id_column` and `name_column`.
#'   Returns an empty `data.frame` with those columns if the required
#'   attributes are absent or all values are missing.
workflow_glue_r_annotation_name_map_from_meta <- function(
    annotation_meta, 
    id_column, 
    name_column, 
    output_id_column
) {
    output_columns <- c(output_id_column, name_column)

    if (!all(c(id_column, name_column) %in% names(annotation_meta))) {
        return(workflow_glue_r_empty_tsv(output_columns))
    }

    # Keep annotation row rows only where both the id and name are present and non-empty.
    keep <- 
        !is.na(annotation_meta[[id_column]]) &
        nzchar(annotation_meta[[id_column]]) &
        !is.na(annotation_meta[[name_column]]) &
        nzchar(annotation_meta[[name_column]])
    if (!any(keep)) {
        return(workflow_glue_r_empty_tsv(output_columns))
    }

    annotation_meta <- annotation_meta[keep, , drop = FALSE]
    annotation_meta <- annotation_meta[!duplicated(annotation_meta[[id_column]]), ]
    out <- data.frame(
        annotation_meta[[id_column]],
        annotation_meta[[name_column]],
        row.names = NULL,
        stringsAsFactors = FALSE
    )
    names(out) <- output_columns
    out
}

#' Extract `GENEID->gene_name` and `TXNAME->transcript_name` mappings.
#'
#' @param annotation_path Path to a GTF or GFF file.
#'
#' @return A list with `gene` and `transcript` data.frames.
workflow_glue_r_annotation_name_maps <- function(annotation_path) {
    if (bambu_missing(annotation_path)) {
        return(list(
            gene = workflow_glue_r_empty_tsv(c("GENEID", "gene_name")),
            transcript = workflow_glue_r_empty_tsv(c("TXNAME", "transcript_name"))
        ))
    }
    annotation <- rtracklayer::import(annotation_path)
    annotation_meta <- S4Vectors::mcols(annotation)
    
    list(
        gene = workflow_glue_r_annotation_name_map_from_meta(
            annotation_meta, "gene_id", "gene_name", "GENEID"
        ),
        transcript = workflow_glue_r_annotation_name_map_from_meta(
            annotation_meta, "transcript_id", "transcript_name", "TXNAME"
        )
    )
}