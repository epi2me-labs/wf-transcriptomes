workflow_glue_r_load <- function(pkg_dir = NULL, env = globalenv()) {
    if (is.null(pkg_dir)) {
        stop("pkg_dir must be provided when loading workflow_glue_r.", call. = FALSE)
    }

    r_dir <- file.path(pkg_dir, "R")
    if (!dir.exists(r_dir)) {
        stop(sprintf("R source directory not found: %s", r_dir), call. = FALSE)
    }

    for (path in sort(list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE))) {
        sys.source(path, envir = env)
    }

    invisible(pkg_dir)
}
