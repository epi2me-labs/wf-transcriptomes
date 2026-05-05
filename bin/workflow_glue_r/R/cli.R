workflow_glue_r_components <- function(env = globalenv()) {
    parser_suffix <- "_arg_parser"
    parser_names <- grep(
        paste0(parser_suffix, "$"),
        ls(envir = env, all.names = TRUE),
        value = TRUE
    )

    components <- list()
    for (parser_name in parser_names) {
        component <- sub(paste0(parser_suffix, "$"), "", parser_name)
        cli_name <- paste0("run_", component, "_cli")
        parser <- get(parser_name, envir = env)
        runner <- if (exists(cli_name, envir = env, mode = "function")) {
            get(cli_name, envir = env)
        } else {
            NULL
        }

        if (is.function(parser) && is.function(runner)) {
            components[[component]] <- list(
                name = component,
                parser_name = parser_name,
                runner_name = cli_name,
                parser = parser,
                runner = runner
            )
        }
    }

    components[sort(names(components))]
}

workflow_glue_r_usage <- function(components = workflow_glue_r_components()) {
    component_names <- names(components)
    lines <- c(
        "Usage: supeRglue <command> [options]",
        "",
        "Commands:",
        if (length(component_names) > 0) {
            paste0("  ", component_names)
        } else {
            "  <none found>"
        },
        "",
        "Use 'supeRglue <command> --help' for command-specific options."
    )
    paste(lines, collapse = "\n")
}

workflow_glue_r_cli <- function(argv = commandArgs(trailingOnly = TRUE), env = globalenv()) {
    components <- workflow_glue_r_components(env = env)

    if (length(argv) < 1 || argv[[1]] %in% c("-h", "--help", "help")) {
        cat(workflow_glue_r_usage(components), "\n")
        return(invisible(0L))
    }

    command <- argv[[1]]
    if (!command %in% names(components)) {
        stop(
            sprintf(
                "Unknown supeRglue command '%s'. Available commands: %s",
                command,
                paste(names(components), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    components[[command]]$runner(argv[-1])
}
