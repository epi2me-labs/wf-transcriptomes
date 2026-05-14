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
