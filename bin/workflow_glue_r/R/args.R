workflow_glue_r_parse_csv_list <- function(value) {
    if (is.null(value) || length(value) == 0 || (length(value) == 1 && is.na(value))) {
        return(character(0))
    }

    values <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
    values[nzchar(values)]
}

workflow_glue_r_flag_present <- function(raw_argv, flag) {
    if (is.null(raw_argv) || length(raw_argv) == 0) {
        return(FALSE)
    }
    any(raw_argv == flag | startsWith(raw_argv, paste0(flag, "=")))
}

workflow_glue_r_arg_parser_from_spec <- function(description, arg_spec) {
    parser <- argparser::arg_parser(description)
    for (arg in arg_spec) {
        add_args <- list(
            parser = parser,
            arg = arg$flag,
            help = arg$help,
            type = arg$type
        )
        if ("default" %in% names(arg)) {
            add_args$default <- arg$default
        }
        parser <- do.call(argparser::add_argument, add_args)
    }
    parser
}

workflow_glue_r_arg_value_error <- function(arg, fallback) {
    if ("value_error" %in% names(arg)) {
        return(arg$value_error)
    }
    fallback
}

workflow_glue_r_scalar_arg <- function(value, arg) {
    if (length(value) != 1) {
        stop(
            sprintf(
                "%s must be a single %s value.",
                arg$flag,
                arg$type
            ),
            call. = FALSE
        )
    }
    value
}

workflow_glue_r_normalise_arg_value <- function(value, arg, flag_provided = FALSE) {
    value_is_na <- length(value) == 1 && is.na(value)
    value_is_absent <- is.null(value) ||
        length(value) == 0 ||
        (value_is_na && !isTRUE(flag_provided)) ||
        (is.character(value) && length(value) == 1 && !nzchar(value))

    if (value_is_absent) {
        if ("default" %in% names(arg)) {
            return(arg$default)
        }
        return(NULL)
    }

    value <- workflow_glue_r_scalar_arg(value, arg)
    if (identical(arg$type, "character")) {
        if (is.na(value)) {
            stop(
                workflow_glue_r_arg_value_error(arg, sprintf("%s must be a non-empty string.", arg$flag)),
                call. = FALSE
            )
        }
        return(as.character(value))
    }

    if (identical(arg$type, "integer")) {
        numeric_value <- suppressWarnings(as.numeric(value))
        integer_value <- suppressWarnings(as.integer(numeric_value))
        if (
            is.na(numeric_value) ||
                is.na(integer_value) ||
                !is.finite(numeric_value) ||
                numeric_value != integer_value
        ) {
            stop(
                workflow_glue_r_arg_value_error(arg, sprintf("%s must be an integer.", arg$flag)),
                call. = FALSE
            )
        }
        return(integer_value)
    }

    if (identical(arg$type, "numeric")) {
        numeric_value <- suppressWarnings(as.numeric(value))
        if (is.na(numeric_value)) {
            stop(
                workflow_glue_r_arg_value_error(arg, sprintf("%s must be numeric.", arg$flag)),
                call. = FALSE
            )
        }
        return(numeric_value)
    }

    if (identical(arg$type, "logical")) {
        logical_value <- suppressWarnings(as.logical(value))
        if (length(logical_value) != 1 || is.na(logical_value)) {
            stop(
                workflow_glue_r_arg_value_error(arg, sprintf("%s must be true or false.", arg$flag)),
                call. = FALSE
            )
        }
        return(logical_value)
    }

    value
}

workflow_glue_r_validate_arg_value <- function(value, arg) {
    if (is.null(value)) {
        return(invisible(value))
    }

    if ("choices" %in% names(arg) && !value %in% arg$choices) {
        stop(
            sprintf(
                "%s must be one of: %s",
                arg$name,
                paste(arg$choices, collapse = ", ")
            ),
            call. = FALSE
        )
    }

    if (
        ("min" %in% names(arg) && value < arg$min) ||
            ("max" %in% names(arg) && value > arg$max)
    ) {
        stop(
            workflow_glue_r_arg_value_error(
                arg,
                sprintf(
                    "%s must be between %s and %s.",
                    arg$flag,
                    if ("min" %in% names(arg)) arg$min else "-Inf",
                    if ("max" %in% names(arg)) arg$max else "Inf"
                )
            ),
            call. = FALSE
        )
    }

    invisible(value)
}

workflow_glue_r_normalise_args <- function(argv, arg_spec, raw_argv = NULL) {
    for (arg in arg_spec) {
        argv[[arg$name]] <- workflow_glue_r_normalise_arg_value(
            argv[[arg$name]],
            arg,
            flag_provided = workflow_glue_r_flag_present(raw_argv, arg$flag)
        )
        workflow_glue_r_validate_arg_value(argv[[arg$name]], arg)
    }

    required_args <- vapply(arg_spec, function(arg) {
        isTRUE(arg$required)
    }, logical(1))
    required_arg_names <- vapply(arg_spec[required_args], `[[`, character(1), "name")
    missing_args <- required_arg_names[vapply(required_arg_names, function(arg_name) {
        is.null(argv[[arg_name]])
    }, logical(1))]

    if (length(missing_args) > 0) {
        stop(
            sprintf(
                "Missing required arguments: %s",
                paste(sprintf("--%s", missing_args), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    xor_args <- vapply(arg_spec, function(arg) {
        "xor_group" %in% names(arg)
    }, logical(1))
    xor_groups <- unique(vapply(arg_spec[xor_args], `[[`, character(1), "xor_group"))
    for (xor_group in xor_groups) {
        group_args <- arg_spec[vapply(arg_spec, function(arg) {
            identical(arg$xor_group, xor_group)
        }, logical(1))]
        group_arg_names <- vapply(group_args, `[[`, character(1), "name")
        present_args <- vapply(group_arg_names, function(arg_name) {
            !is.null(argv[[arg_name]])
        }, logical(1))

        if (sum(present_args) != 1) {
            group_flags <- vapply(group_args, `[[`, character(1), "flag")
            stop(
                sprintf(
                    "Provide exactly one of %s.",
                    paste(group_flags, collapse = " or ")
                ),
                call. = FALSE
            )
        }
    }

    argv
}
